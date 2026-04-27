# ============================================================
# 06_coulombe_2SRR_pipeline.R   (v5 — PhD-grade)
#
# Pipeline Coulombe (2024) alinhado 1:1 com o script original
# hugocout/forecasting_table16to17.R  (mod=2, 2SRR).
#
# CORREÇÕES v5 em relação à v4:
#   [C1] UMA única PCA por iteração t (usando data[1:(t-h),:])
#        → elimina inconsistência de espaço fatorial entre
#          reg e last. No Hugo: factors = data[1:(t-h),:] puro.
#   [C2] `last` extraído de `train` ANTES do corte temporal,
#        exatamente como Hugo:
#          last  = train[nrow(train), ]
#          train = train[1:(nrow(train)-h), ]
#   [C3] Remoção do bloco `nd <- h-1` que removia linhas
#        corretas já alinhadas pelo make_reg_matrix().
#   [C4] forecast[t] = OF(aa$fcast) - last[1]
#        → subtrai last[1] = y_{t-h} acumulado, tal como Hugo.
#        Para h=1 isso equivale à diferença de nível;
#        para h>1 é a diferença acumulada — semanticamente
#        correta porque forecast_vars já é cumulativa.
#   [C5] closeAllConnections() antes do loop resolve o limite
#        de 128 conexões no Windows (problema operacional,
#        zero impacto nos resultados econométricos).
#   [C6] nf=2 fatores EM (como Hugo mod=2 que usa EM_sw n=2)
#        em vez de PCA convencional com nf=8.
#   [C7] Parâmetros do Hugo preservados exatamente:
#        ly=2, lf=2, lambdavec=exp(linspace(-2,12,15)),
#        alpha=0.01, kfold=5, tol=1e-6, maxit=10.
#
# Estrutura esperada:
#   forecasting_inflation/
#     coulombe/         <- funções do Hugo
#     data/data.rda     <- base do Medeiros
#     forecasts/        <- saída .rda e .csv
#     results/          <- CSVs auxiliares
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
# 2. CARREGA FUNÇÕES DO COULOMBE
# ============================================================
cs <- function(f) {
  p <- file.path("coulombe", f)
  if (!file.exists(p))
    stop(paste0("Arquivo nao encontrado: ", p))
  source(p)
  cat(sprintf("  [OK] %s\n", f))
}

cat("=== Carregando funcoes Coulombe ===\n")
cs("EM_sw.R")
cs("ICp2.R")
cs("Xgenerators_v190127.R")    # make_reg_matrix()
cs("dualGRRmdA_v190215.R")
cs("CVGSBHK_v181127.R")
cs("zfun_v190304.R")
cs("factor.R")
cs("TVPRRcosso_v181120.R")     # TVPRR_cosso() — nucleo do 2SRR
cs("TVPRR_v181111.R")
cs("fastZrot_v181125.R")
cs("CVKFMV_v190214.R")
cat("Todas as funcoes carregadas.\n\n")

# ============================================================
# 3. FILTRO DE OUTLIERS — identico ao Hugo
#    OF(pred, y=train[subset,1], tol=2, go.to.pred=pred.lin)
# ============================================================
OF <- function(pred, y, tol = 2, go.to.pred) {
  newx           <- pred
  cond.max       <- (newx - mean(y)) > tol * (max(y) - mean(y))
  cond.min       <- (newx - mean(y)) < tol * (min(y) - mean(y))
  newx[cond.max] <- go.to.pred[cond.max]
  newx[cond.min] <- go.to.pred[cond.min]
  return(newx)
}

# ============================================================
# 4. CARREGA BASE DO MEDEIROS
# ============================================================
load("data/data.rda")

fred_raw <- as.data.frame(data)
bigt     <- nrow(fred_raw)

date_col <- grep("^date$",     colnames(fred_raw), ignore.case = TRUE)[1]
cpi_col  <- grep("^CPIAUCSL$", colnames(fred_raw), ignore.case = TRUE)[1]

if (is.na(date_col) || is.na(cpi_col)) {
  cat("Colunas disponiveis:\n"); print(colnames(fred_raw))
  stop("Nao encontrou 'date' ou 'CPIAUCSL'. Ajuste manualmente.")
}

dates <- fred_raw[, date_col]
y_raw <- fred_raw[, cpi_col]          # nivel de CPIAUCSL transformado

# Matriz de preditores (sem date e sem CPI — o CPI entra via Y em make_reg_matrix)
X_raw <- as.matrix(fred_raw[, -c(date_col, cpi_col)])

cat(sprintf("Base: %d obs x %d preditores | %s a %s\n",
            bigt, ncol(X_raw),
            as.character(dates[1]), as.character(dates[bigt])))
cat(sprintf("Variavel alvo: '%s'\n\n", colnames(fred_raw)[cpi_col]))

# ============================================================
# 5. VARIAVEL Y ACUMULADA h-PASSOS (direct multi-step forecast)
#
#    forecast_vars[t, h] = sum(y[t+1], ..., y[t+h])
#
#    Equivalente a newQ_targets.csv do Hugo:
#    a variavel dependente ja carrega a acumulacao de h passos.
#    make_reg_matrix alinha y[t] com X[t-h] internamente.
# ============================================================
build_cumulative_y <- function(y, h) {
  n  <- length(y)
  yh <- rep(NA_real_, n)
  for (t in seq_len(n - h)) yh[t] <- sum(y[(t + 1L):(t + h)])
  yh
}

hor           <- c(1, 3, 6, 12)
forecast_vars <- sapply(hor, build_cumulative_y, y = y_raw)
colnames(forecast_vars) <- paste0("h", hor)

# ============================================================
# 6. PARAMETROS POOS
#
#    Hugo: tau=158, bigt=222 → ~71% in-sample (trimestral)
#    Medeiros: nwindows=312 (mantido do 03_call_model_felipe.R)
#    ly=2, lf=2, nf=2 (EM_sw n=2, como mod=2 do Hugo)
#    lambdavec = exp(linspace(-2,12,15)) — identico ao Hugo
# ============================================================
nwindows  <- 312
tau       <- bigt - nwindows
n_oos     <- nwindows
ly        <- 2
lf        <- 2
nf_em     <- 2     # [C6] fatores via EM_sw, nao PCA — como Hugo mod=2
lambdavec <- exp(pracma::linspace(-2, 12, n = 15))
silenceplz <- 1

cat(sprintf("POOS: bigt=%d | tau=%d | n_oos=%d | ly=%d | lf=%d | nf_em=%d\n\n",
            bigt, tau, n_oos, ly, lf, nf_em))

dir.create("forecasts",     showWarnings = FALSE)
dir.create("results",       showWarnings = FALSE)

# ============================================================
# 7. ARRAYS DE RESULTADO
# ============================================================
fc_ridge  <- matrix(NA_real_, nrow = bigt, ncol = length(hor))
fc_2srr   <- matrix(NA_real_, nrow = bigt, ncol = length(hor))
lam_ridge <- matrix(NA_real_, nrow = bigt, ncol = length(hor))
lam1_2srr <- matrix(NA_real_, nrow = bigt, ncol = length(hor))
lam2_2srr <- matrix(NA_real_, nrow = bigt, ncol = length(hor))

betas_2srr <- setNames(vector("list", length(hor)), paste0("h", hor))
for (hi in seq_along(hor)) betas_2srr[[hi]] <- vector("list", n_oos)

# ============================================================
# 8. LOOP POOS PRINCIPAL
#
#    Logica identica ao Hugo (forecasting_table16to17.R, mod=2):
#
#    for t in tau:bigt
#      data    = fred2[1:(t-h), ]           <- info set t-h
#      y       = forecast_vars[1:(t-h), v]  <- y acumulada
#      Y       = forecast_vars[1:(t-h), 1]  <- y nivel (lags)
#      factors = EM_sw(data, n=2)$factors   <- 2 fatores EM
#      start   = sum(is.na(y)) + 1
#      train   = make_reg_matrix(y,Y,factors,h,ly,lf)
#      last    = train[nrow(train), ]        <- ultima linha
#      train   = train[1:(nrow(train)-h), ]  <- remove h linhas finais
#      train   = train[(maxlag+1):nrow, ]    <- remove lags iniciais
#      m=1: Ridge plain
#      forecast[t] = predict(mdl, last[-1]) - last[1]
#      m=2: TVPRR_cosso (2SRR)
#      forecast[t] = OF(aa$fcast) - last[1]
# ============================================================

# [C5] Fecha conexoes orfas antes do loop (fix Windows)
closeAllConnections()

cat("=== INICIANDO LOOP POOS ===\n")
t0_total <- proc.time()

for (hi in seq_along(hor)) {
  h <- hor[hi]
  cat(sprintf("\n--- Horizonte h = %d ---\n", h))

  # Posicao do target acumulado (h=1 usa col 1 de forecast_vars, etc.)
  # Hugo usa hor_pos = c(1,2,3,3) mas aqui cada h tem sua propria coluna
  y_target <- forecast_vars[, hi]   # y acumulada para este h

  for (t in tau:(bigt - 1)) {
    idx <- t - tau + 1L

    # -----------------------------------------------------------
    # INFO SET: apenas observacoes ate t-h
    # Identico ao Hugo: data = fred2[1:(t-h), ]
    # -----------------------------------------------------------
    T_end <- t - h
    if (T_end < (ly + lf + 10L)) next

    # [C1] UMA unica PCA/EM por iteracao, sobre X[1:(t-h)]
    #      Hugo mod=2: factors = EM_sw(data[,-vars[v]], n=2)$factors
    #      Aqui: EM_sw sobre X_raw[1:T_end] (todos os preditores menos CPI)
    data_th <- X_raw[1:T_end, , drop = FALSE]

    factors_em <- tryCatch(
      EM_sw(data = as.data.frame(data_th), n = nf_em, it_max = 1000)$factors,
      error = function(e) NULL
    )
    if (is.null(factors_em)) next
    factors_em <- as.matrix(factors_em)

    # Variaveis dependentes alinhadas ao mesmo info set
    y_h <- as.matrix(y_target[1:T_end])
    Y_h <- as.matrix(y_raw[1:T_end])       # nivel de y (para lags em make_reg_matrix)

    # Remove NAs iniciais (identico ao Hugo: start = sum(is.na(y)) + 1)
    start <- sum(is.na(y_h)) + 1L
    if (start > T_end || (T_end - start + 1L) < 20L) next

    y_is  <- y_h[start:T_end, , drop = FALSE]
    Y_is  <- Y_h[start:T_end, , drop = FALSE]
    f_is  <- factors_em[start:T_end, , drop = FALSE]

    # -----------------------------------------------------------
    # make_reg_matrix — constroi [y_h | const | lags_Y | lags_f]
    # [C2] last extraido ANTES do corte de h linhas finais
    # -----------------------------------------------------------
    train <- tryCatch(
      make_reg_matrix(y = y_is, Y = Y_is, factors = f_is,
                      h = h, ly = ly, lf = lf),
      error = function(e) NULL
    )
    if (is.null(train) || nrow(train) < (h + ly + lf + 5L)) next

    # [C2] last = ultima linha de train (inclui obs de previsao)
    last  <- as.numeric(train[nrow(train), ])

    # [C2] Remove h linhas finais (como Hugo: train[1:(nrow-h), ])
    # [C3] NAO remove nd adicional — make_reg_matrix ja alinha corretamente
    train <- train[1:(nrow(train) - h), , drop = FALSE]

    # Remove lags iniciais (como Hugo: train[(maxlag+1):nrow, ])
    maxlag <- max(lf, ly)
    if (nrow(train) <= maxlag) next
    train <- train[(maxlag + 1L):nrow(train), , drop = FALSE]

    # Remove linhas com NA
    train <- train[complete.cases(train), , drop = FALSE]
    if (nrow(train) < 15L || ncol(train) < 2L) next

    subset <- seq_len(nrow(train))

    # -----------------------------------------------------------
    # m=1: RIDGE PLANO (Plain Ridge)
    # Identico ao Hugo: CV = cv.glmnet(...); predict - last[1]
    # -----------------------------------------------------------
    CV <- tryCatch(
      cv.glmnet(x      = train[subset, -1, drop = FALSE],
                y      = train[subset,  1],
                family = "gaussian",
                alpha  = 0),
      error = function(e) NULL
    )
    if (is.null(CV)) next

    mdl_r   <- glmnet(x      = train[subset, -1, drop = FALSE],
                      y      = train[subset,  1],
                      family = "gaussian",
                      alpha  = 0,
                      lambda = CV$lambda.min)
    pred_lin <- as.numeric(predict(mdl_r, newx = matrix(last[-1], nrow = 1)))

    # [C4] forecast = predict - last[1], identico ao Hugo
    fc_ridge[t, hi]   <- pred_lin - last[1]
    lam_ridge[t, hi]  <- CV$lambda.min

    # -----------------------------------------------------------
    # m=2: 2SRR — TVPRR_cosso
    # Identico ao Hugo: aa <- TVPRR_cosso(...); OF(aa$fcast) - last[1]
    # -----------------------------------------------------------
    aa <- tryCatch(
      TVPRR_cosso(
        y         = train[subset, 1],
        X         = train[subset, -1, drop = FALSE],
        lambdavec = lambdavec,
        sweigths  = 1,
        type      = 2,
        alpha     = 0.01,
        silent    = silenceplz,
        kfold     = 5,
        lambda2   = CV$lambda.min,
        tol       = 1e-6,
        maxit     = 10,
        oosX      = last[-1]
      ),
      error = function(e) NULL
    )

    if (!is.null(aa)) {
      # Previsao bruta do 2SRR
      p_raw <- if (!is.null(aa$fcast)) {
        as.numeric(aa$fcast)
      } else {
        # Fallback: produto interno do ultimo beta com last
        bm   <- aa$grrats$betas_grr
        bl   <- if (is.matrix(bm)) bm[nrow(bm), ] else as.numeric(bm)
        sum(c(1, last[-1]) * bl)
      }

      # Filtro de outliers (identico ao Hugo)
      p_filt <- OF(pred      = p_raw,
                   y         = train[subset, 1],
                   go.to.pred = pred_lin)

      # [C4] forecast = OF(aa$fcast) - last[1], identico ao Hugo
      fc_2srr[t, hi]   <- p_filt - last[1]
      lam1_2srr[t, hi] <- CV$lambda.min
      lam2_2srr[t, hi] <- if (!is.null(aa$grrats$lambdas)) aa$grrats$lambdas[1] else NA_real_

      # Salva betas TVP para analise posterior
      betas_2srr[[hi]][[idx]] <- list(
        t     = t,
        date  = dates[t],
        betas = aa$grrats$betas_grr
      )
    }

    if (idx %% 24L == 0L) {
      el <- (proc.time() - t0_total)["elapsed"]
      cat(sprintf("  h=%d | t=%d/%d (%.0f%%) | %.1f min\n",
                  h, t, bigt - 1L, 100 * idx / n_oos, el / 60))
    }
  }  # fim loop t

  cat(sprintf("  h=%d concluido.\n", h))
}  # fim loop hi

cat(sprintf("\nPOOS completo: %.1f min\n\n",
            (proc.time() - t0_total)["elapsed"] / 60))

# ============================================================
# 9. SALVA RESULTADOS
# ============================================================
save(fc_ridge,                    file = "forecasts/coulombe_ridge.rda")
save(fc_2srr,                     file = "forecasts/coulombe_2SRR.rda")
save(lam_ridge, lam1_2srr, lam2_2srr, file = "forecasts/coulombe_lambdas.rda")
save(betas_2srr,                  file = "forecasts/coulombe_betas_2SRR.rda")

oos_idx <- (tau + 1L):bigt

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
  cat(sprintf("Salvo: %s  (ridge=%d validos | 2SRR=%d validos)\n",
              fname,
              sum(!is.na(df_out$fc_ridge)),
              sum(!is.na(df_out$fc_2srr))))
}

# CSV de betas para h=1 (para 07_compare.R)
hi1     <- which(hor == 1L)
valid_b <- Filter(Negate(is.null), betas_2srr[[hi1]])
if (length(valid_b) > 0L) {
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

cat("\n=== 06_coulombe_2SRR_pipeline.R v5 — COMPLETO ===\n")
