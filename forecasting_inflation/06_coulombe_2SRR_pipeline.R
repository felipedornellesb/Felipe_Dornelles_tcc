# ============================================================
# 06_coulombe_2SRR_pipeline.R   (v4 — final)
#
# Pipeline Coulombe (2024) FIEL ao hugocout, adaptado para
# a base do Medeiros (FRED-MD mensal, EUA).
#
# Estrutura do repo:
#   forecasting_inflation/
#     coulombe/         <- funções do Hugo (CVGSBHK, TVPRRcosso, etc.)
#     functions/        <- funções do Medeiros (functions.R, etc.)
#     data/data.rda     <- base carregada aqui (data frame com 'date' + vars FRED)
#     forecasts/        <- output .rda e .csv
#     results/          <- output CSVs auxiliares
#
# Ordem de execução:
#   01 → 02 → 03_call_model_felipe.R → 06_coulombe_2SRR_pipeline.R → 07_compare.R
#
# O que este script faz diferente do 03_call_model_felipe.R:
#   1. Usa make_reg_matrix() + make_last()  (Xgenerators_v190127.R)
#   2. Usa EM_sw() para imputação de NAs    (Stock & Watson)
#   3. Variável Y acumulada h-passos        (direct multi-step)
#   4. Filtro de outliers OF()              (idêntico ao Hugo)
#   5. lambda_vec = exp(linspace(-2,12,15)) (idêntico ao Hugo)
#   6. Ridge plano estimado primeiro → lambda2 para o 2SRR
#
# Modelos: m=1 Ridge | m=2 2SRR
# Horizontes: hor = c(1, 3, 6, 12)
# ============================================================

rm(list = ls())
setwd("~/tcc/Felipe_Dornelles_tcc/forecasting_inflation")
cat("Working dir:", getwd(), "\n")

# ============================================================
# 0. PACOTES
# ============================================================
pkgs <- c("pracma", "glmnet", "timeSeries", "matrixcalc", "GA", "e1071")
new  <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(new)) install.packages(new, repos = "http://cran.us.r-project.org")
invisible(lapply(pkgs, library, character.only = TRUE))

# ============================================================
# 1. BAIXA FUNÇÕES FALTANTES DE coulombe/
#    (EM_sw, ICp2, factor, TVPRR_v181111, CVKFMV)
# ============================================================
base_raw <- paste0(
  "https://raw.githubusercontent.com/hugocout/",
  "Replication-codes-for-Time-Varying-Parameters-as-Ridge-Regressions/",
  "main/Empirical/20_tools"
)

extras <- list(
  list(url = paste0(base_raw, "/EM_sw.R"),                       dest = "coulombe/EM_sw.R"),
  list(url = paste0(base_raw, "/ICp2.R"),                        dest = "coulombe/ICp2.R"),
  list(url = paste0(base_raw, "/functions/factor.R"),            dest = "coulombe/factor.R"),
  list(url = paste0(base_raw, "/functions/TVPRR_v181111.R"),     dest = "coulombe/TVPRR_v181111.R"),
  list(url = paste0(base_raw, "/functions/CVKFMV_v190214.R"),    dest = "coulombe/CVKFMV_v190214.R")
)

for (f in extras) {
  if (!file.exists(f$dest)) {
    cat(sprintf("  Baixando %s ...", basename(f$dest)))
    tryCatch({
      download.file(f$url, destfile = f$dest, quiet = TRUE, method = "libcurl")
      cat(" OK\n")
    }, error = function(e) cat(sprintf(" ERRO: %s\n", e$message)))
  }
}

# ============================================================
# 2. CARREGA TODAS AS FUNÇÕES DO COULOMBE (coulombe/ plano)
# ============================================================
cs <- function(f) {
  p <- file.path("coulombe", f)
  if (!file.exists(p))
    stop(paste0("Arquivo não encontrado: ", p,
                "\n  Baixe de github.com/hugocout/.../Empirical/20_tools/"))
  source(p)
  cat(sprintf("  [OK] %s\n", f))
}

cat("=== Carregando funções Coulombe ===\n")
cs("EM_sw.R")
cs("ICp2.R")
cs("Xgenerators_v190127.R")    # make_reg_matrix() + make_last()
cs("dualGRRmdA_v190215.R")
cs("CVGSBHK_v181127.R")
cs("zfun_v190304.R")
cs("factor.R")
cs("TVPRRcosso_v181120.R")     # TVPRR_cosso() — núcleo do 2SRR
cs("TVPRR_v181111.R")
cs("fastZrot_v181125.R")
cs("CVKFMV_v190214.R")
cat("Todas as funções carregadas.\n\n")

# ============================================================
# 3. FILTRO DE OUTLIERS (idêntico ao de Hugo, Empirical_v2.R)
# ============================================================
OF <- function(pred, y, tol = 2, go.to.pred) {
  newx          <- pred
  cond.max      <- (newx - mean(y)) > tol * (max(y) - mean(y))
  cond.min      <- (newx - mean(y)) < tol * (min(y) - mean(y))
  newx[cond.max] <- go.to.pred[cond.max]
  newx[cond.min] <- go.to.pred[cond.min]
  return(newx)
}

# ============================================================
# 4. CARREGA DATA DO MEDEIROS
#    data.rda -> data.frame: coluna 'date' + variáveis FRED-MD
# ============================================================
load("data/data.rda")

fred_raw <- as.data.frame(data)
bigt     <- nrow(fred_raw)

date_col <- grep("^date$", colnames(fred_raw), ignore.case = TRUE)
cpi_col  <- grep("^CPIAUCSL$", colnames(fred_raw), ignore.case = TRUE)

if (length(date_col) == 0 || length(cpi_col) == 0) {
  cat("Colunas disponíveis:\n"); print(colnames(fred_raw))
  stop("Nao encontrou 'date' ou 'CPIAUCSL'. Ajuste cpi_col/date_col manualmente.")
}
cpi_col  <- cpi_col[1]
date_col <- date_col[1]
dates    <- fred_raw[, date_col]
y_raw    <- fred_raw[, cpi_col]

cat(sprintf("Base: %d obs x %d vars | Periodo: %s a %s\n",
            bigt, ncol(fred_raw),
            as.character(dates[1]), as.character(dates[bigt])))
cat(sprintf("CPI: coluna '%s' (idx %d)\n\n",
            colnames(fred_raw)[cpi_col], cpi_col))

X_raw <- as.matrix(fred_raw[, -c(date_col, cpi_col)])

# ============================================================
# 5. IMPUTAÇÃO EM (Stock & Watson)
# ============================================================
cat("Imputação EM Stock & Watson (n=8, it_max=1000)...\n")
X_imp <- tryCatch({
  em_out <- EM_sw(data = as.data.frame(X_raw), n = 8, it_max = 1000)
  as.matrix(em_out$data)
}, error = function(e) {
  cat(sprintf("  EM_sw falhou (%s) — usando interpolação linear\n", e$message))
  X_r <- X_raw
  for (j in seq_len(ncol(X_r))) {
    nas <- which(is.na(X_r[, j]))
    if (length(nas) > 0 && length(nas) < nrow(X_r) - 2)
      X_r[, j] <- approx(seq_len(nrow(X_r)), X_r[, j], xout = seq_len(nrow(X_r)))$y
  }
  X_r
})
cat("Imputação concluída.\n\n")

# ============================================================
# 6. VARIÁVEL Y ACUMULADA h-PASSOS (direct multi-step)
#    y_h[t] = sum(y[t+1], ..., y[t+h])
# ============================================================
build_cumulative_y <- function(y, h) {
  n  <- length(y)
  yh <- rep(NA, n)
  for (t in seq_len(n - h)) yh[t] <- sum(y[(t + 1):(t + h)])
  yh
}

hor           <- c(1, 3, 6, 12)
forecast_vars <- sapply(hor, function(h) build_cumulative_y(y_raw, h))
colnames(forecast_vars) <- paste0("h", hor)

# ============================================================
# 7. PARÂMETROS POOS
#    nwindows = 312 (mesmo do 03_call_model_felipe.R)
# ============================================================
nwindows  <- 312
tau       <- bigt - nwindows
n_oos     <- nwindows
nf        <- 8
ly        <- 2
lf        <- 2
lambdavec <- exp(pracma::linspace(-2, 12, n = 15))

cat(sprintf("POOS: bigt=%d | tau=%d | n_oos=%d | nf=%d | ly=%d | lf=%d\n\n",
            bigt, tau, n_oos, nf, ly, lf))

dir.create("forecasts",     showWarnings = FALSE)
dir.create("results/plots", showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 8. ARRAYS DE RESULTADO
# ============================================================
fc_ridge  <- matrix(NA, nrow = bigt, ncol = length(hor))
fc_2srr   <- matrix(NA, nrow = bigt, ncol = length(hor))
lam_ridge <- matrix(NA, nrow = bigt, ncol = length(hor))
lam1_2srr <- matrix(NA, nrow = bigt, ncol = length(hor))
lam2_2srr <- matrix(NA, nrow = bigt, ncol = length(hor))

betas_2srr <- setNames(vector("list", length(hor)), paste0("h", hor))
for (hi in seq_along(hor)) betas_2srr[[hi]] <- vector("list", n_oos)

# ============================================================
# 9. LOOP POOS PRINCIPAL
# ============================================================
cat("=== INICIANDO LOOP POOS ===\n")
t0_total <- proc.time()

for (hi in seq_along(hor)) {
  h <- hor[hi]
  cat(sprintf("\n--- Horizonte h = %d ---\n", h))

  for (t in tau:(bigt - 1)) {
    idx   <- t - tau + 1
    T_end <- t - h          # information set: apenas até t-h
    if (T_end < (ly + lf + 20)) next

    # ---- PCA sobre X imputada (in-sample 1:T_end) ----
    X_is  <- X_imp[1:T_end, , drop = FALSE]
    pc_is <- prcomp(X_is, center = TRUE, scale. = TRUE)
    fac   <- pc_is$x[, 1:min(nf, ncol(pc_is$x)), drop = FALSE]

    y_h  <- as.matrix(forecast_vars[1:T_end, hi])
    na_end <- suppressWarnings(max(which(is.na(y_h))))
    si     <- if (is.finite(na_end)) na_end + 1L else 1L
    if ((T_end - si + 1) < 20) next

    y_is <- y_h[si:T_end, , drop = FALSE]
    Y_is <- as.matrix(y_raw[si:T_end])
    f_is <- fac[si:T_end, , drop = FALSE]

    # ---- make_reg_matrix ----
    reg <- tryCatch(
      make_reg_matrix(y = y_is, Y = Y_is, factors = f_is,
                      h = h, ly = ly, lf = lf),
      error = function(e) NULL
    )
    if (is.null(reg)) next

    ml  <- max(lf, ly)
    reg <- reg[(ml + 1):nrow(reg), , drop = FALSE]
    nd  <- h - 1
    if (nd > 0 && nrow(reg) > nd) reg <- reg[1:(nrow(reg) - nd), , drop = FALSE]
    reg <- reg[complete.cases(reg), , drop = FALSE]
    if (nrow(reg) < 15 || ncol(reg) < 2) next

    # ---- make_last: vetor regressores para previsão em t ----
    X_full   <- X_imp[1:t, , drop = FALSE]
    pc_full  <- prcomp(X_full, center = TRUE, scale. = TRUE)
    fac_full <- pc_full$x[, 1:min(nf, ncol(pc_full$x)), drop = FALSE]
    Y_full   <- as.matrix(y_raw[si:t])
    f_full   <- fac_full[si:t, , drop = FALSE]

    last <- tryCatch(
      make_last(y = y_is, Y = Y_full, factors = f_full,
                h = h, ly = ly, lf = lf),
      error = function(e) NULL
    )
    if (is.null(last) || any(is.na(last))) next
    last <- as.numeric(last)

    # ============================================================
    # m=1: RIDGE PLANO
    # ============================================================
    CV <- tryCatch(
      cv.glmnet(x = reg[, -1, drop = FALSE],
                y = reg[, 1],
                family = "gaussian", alpha = 0),
      error = function(e) NULL
    )
    if (is.null(CV)) next

    mdl_r   <- glmnet(x = reg[, -1, drop = FALSE],
                      y = reg[, 1],
                      family = "gaussian", alpha = 0,
                      lambda = CV$lambda.min)
    p_ridge <- as.numeric(predict(mdl_r, newx = matrix(last, nrow = 1)))

    y_ref           <- y_raw[T_end]
    fc_ridge[t, hi] <- p_ridge - y_ref
    lam_ridge[t, hi] <- CV$lambda.min

    # ============================================================
    # m=2: 2SRR (Coulombe — TVPRR_cosso)
    # ============================================================
    aa <- tryCatch(
      TVPRR_cosso(
        y         = reg[, 1],
        X         = reg[, -1, drop = FALSE],
        lambdavec = lambdavec,
        sweigths  = 1,
        type      = 2,
        alpha     = 0.01,
        silent    = 1,
        kfold     = 5,
        lambda2   = CV$lambda.min,
        tol       = 1e-6,
        maxit     = 10,
        oosX      = last
      ),
      error = function(e) NULL
    )

    if (!is.null(aa)) {
      p_raw <- if (!is.null(aa$fcast)) {
        as.numeric(aa$fcast)
      } else {
        bm <- aa$grrats$betas_grr
        bl <- if (is.matrix(bm)) bm[nrow(bm), ] else as.numeric(bm)
        sum(c(1, last) * bl)
      }

      p_filt <- OF(pred = p_raw, y = reg[, 1], go.to.pred = p_ridge)

      fc_2srr[t, hi]    <- p_filt - y_ref
      lam1_2srr[t, hi]  <- CV$lambda.min
      lam2_2srr[t, hi]  <- if (!is.null(aa$grrats$lambdas)) aa$grrats$lambdas[1] else NA

      betas_2srr[[hi]][[idx]] <- list(
        t     = t,
        date  = dates[t],
        betas = aa$grrats$betas_grr
      )
    }

    if (idx %% 24 == 0) {
      el <- (proc.time() - t0_total)["elapsed"]
      cat(sprintf("  h=%d | t=%d/%d (%.0f%%) | %.1f min\n",
                  h, t, bigt - 1, 100 * idx / n_oos, el / 60))
    }
  }
  cat(sprintf("  h=%d concluído.\n", h))
}

cat(sprintf("\nPOOS completo: %.1f min\n\n",
            (proc.time() - t0_total)["elapsed"] / 60))

# ============================================================
# 10. SALVA RESULTADOS
# ============================================================
save(fc_ridge,  file = "forecasts/coulombe_ridge.rda")
save(fc_2srr,   file = "forecasts/coulombe_2SRR.rda")
save(lam_ridge, lam1_2srr, lam2_2srr,
     file = "forecasts/coulombe_lambdas.rda")
save(betas_2srr, file = "forecasts/coulombe_betas_2SRR.rda")

oos_idx <- (tau + 1):bigt
for (hi in seq_along(hor)) {
  h <- hor[hi]
  df_out <- data.frame(
    date      = dates[oos_idx],
    realized  = forecast_vars[oos_idx, hi],
    fc_ridge  = fc_ridge[oos_idx, hi],
    fc_2srr   = fc_2srr[oos_idx, hi],
    lam_ridge = lam_ridge[oos_idx, hi],
    lam2_2srr = lam2_2srr[oos_idx, hi]
  )
  fname <- sprintf("forecasts/coulombe_fc_h%02d.csv", h)
  write.csv(df_out, file = fname, row.names = FALSE)
  cat(sprintf("Salvo: %s  (valid ridge=%d | 2SRR=%d)\n",
              fname,
              sum(!is.na(df_out$fc_ridge)),
              sum(!is.na(df_out$fc_2srr))))
}

# CSV de betas para h=1
hi1      <- which(hor == 1)
valid_b  <- Filter(Negate(is.null), betas_2srr[[hi1]])
if (length(valid_b) > 0) {
  beta_rows <- lapply(valid_b, function(b) {
    bm   <- b$betas
    bvec <- if (is.matrix(bm)) bm[nrow(bm), ] else as.numeric(bm)
    c(date = as.character(b$date), t = b$t, bvec)
  })
  df_betas_h1 <- as.data.frame(do.call(rbind, beta_rows))
  write.csv(df_betas_h1, "results/coulombe_betas_h1.csv", row.names = FALSE)
  cat(sprintf("Betas h=1 salvos: results/coulombe_betas_h1.csv (%d linhas)\n",
              nrow(df_betas_h1)))
}

cat("\n=== 06_coulombe_2SRR_pipeline.R v4 — COMPLETO ===\n")
