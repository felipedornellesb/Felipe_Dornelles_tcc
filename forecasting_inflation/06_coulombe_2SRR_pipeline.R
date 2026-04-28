# ============================================================
# 06_coulombe_2SRR_pipeline.R   (v8.0)
#
# CORRECOES vs v7.0:
#   [C1] fGarch adicionado aos pacotes (garchFit necessario
#        para TVPRR_cosso type=2)
#   [C2] timeSeries adicionado (dependencia de fGarch)
#   [C3] Erro capturado com mensagem detalhada para debug
#   [C4] Exporta betas para TODOS os horizontes (nao so h=1)
#   [C5] Teste pre-loop de TVPRR_cosso com dados sinteticos
#   [C6] factor() com fallback seguro
# ============================================================

rm(list = ls())
gc()

setwd("~/tcc/Felipe_Dornelles_tcc/forecasting_inflation")
cat("============================================================\n")
cat("  2SRR Pipeline v8.0\n")
cat("  Working dir:", getwd(), "\n")
cat("============================================================\n\n")

# ============================================================
# 0. PACOTES — INCLUINDO fGarch
# ============================================================
pkgs <- c("pracma", "glmnet", "matrixcalc", "GA", "e1071",
          "fGarch", "timeSeries")
new  <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if (length(new)) {
  cat("Instalando pacotes faltantes:", paste(new, collapse=", "), "\n")
  install.packages(new, repos = "https://cran.r-project.org")
}
invisible(lapply(pkgs, library, character.only = TRUE))
cat("[OK] Pacotes carregados (incluindo fGarch para garchFit)\n")

# Verifica critico
if (!exists("garchFit", mode = "function")) {
  stop("garchFit nao encontrada mesmo apos library(fGarch). Verifique instalacao.")
}
cat("[OK] garchFit disponivel\n\n")

# ============================================================
# 1. factor() PCA
# ============================================================
base_factor_backup <- base::factor
factor <- function(X, n_fac = NULL, ...) {
  if (!is.null(n_fac)) {
    X    <- as.matrix(X)
    Tobs <- nrow(X)
    S    <- (1 / Tobs) * t(X) %*% X
    eig  <- eigen(S, symmetric = TRUE)
    nf   <- min(n_fac, ncol(X))
    lam  <- eig$vectors[, 1:nf, drop = FALSE]
    fac  <- X %*% lam
    fit  <- fac %*% t(lam)
    mse  <- mean((X - fit)^2, na.rm = TRUE)
    return(list(factors = fac, lambda = lam, mse = mse))
  }
  base_factor_backup(X, ...)
}
cat("[OK] factor() PCA definida\n")

# ============================================================
# 2. DOWNLOAD E SOURCE DAS FUNCOES COULOMBE
# ============================================================
dir.create("coulombe", showWarnings = FALSE, recursive = TRUE)

base_url_tools <- paste0(
  "https://raw.githubusercontent.com/hugocout/",
  "Replication-codes-for-Time-Varying-Parameters-as-Ridge-Regressions/",
  "main/Empirical/20_tools"
)
base_url_func <- paste0(base_url_tools, "/functions")

download_list <- list(
  list(url = paste0(base_url_tools, "/EM_sw.R"),           dest = "coulombe/EM_sw.R"),
  list(url = paste0(base_url_tools, "/ICp2.R"),            dest = "coulombe/ICp2.R"),
  list(url = paste0(base_url_func, "/factor.R"),           dest = "coulombe/factor.R"),
  list(url = paste0(base_url_func, "/TVPRR_v181111.R"),    dest = "coulombe/TVPRR_v181111.R"),
  list(url = paste0(base_url_func, "/CVKFMV_v190214.R"),   dest = "coulombe/CVKFMV_v190214.R"),
  list(url = paste0(base_url_func, "/Xgenerators_v190127.R"), dest = "coulombe/Xgenerators_v190127.R"),
  list(url = paste0(base_url_func, "/dualGRRmdA_v190215.R"),  dest = "coulombe/dualGRRmdA_v190215.R"),
  list(url = paste0(base_url_func, "/CVGSBHK_v181127.R"),     dest = "coulombe/CVGSBHK_v181127.R"),
  list(url = paste0(base_url_func, "/zfun_v190304.R"),         dest = "coulombe/zfun_v190304.R"),
  list(url = paste0(base_url_func, "/TVPRRcosso_v181120.R"),   dest = "coulombe/TVPRRcosso_v181120.R"),
  list(url = paste0(base_url_func, "/fastZrot_v181125.R"),     dest = "coulombe/fastZrot_v181125.R")
)

cat("\n=== Download dos scripts Coulombe ===\n")
for (f in download_list) {
  if (!file.exists(f$dest)) {
    cat(sprintf("  Baixando %-40s", basename(f$dest)))
    tryCatch({
      download.file(f$url, destfile = f$dest, quiet = TRUE, method = "libcurl")
      cat(" [OK]\n")
    }, error = function(e) cat(sprintf(" [ERRO: %s]\n", e$message)))
  }
}

# Source na ordem correta
cs <- function(f) {
  p <- file.path("coulombe", f)
  if (!file.exists(p)) { warning(paste0("Nao encontrado: ", p)); return(invisible(NULL)) }
  source(p, local = FALSE)
  cat(sprintf("  [OK] %s\n", f))
}

cat("\n=== Carregando funcoes Coulombe ===\n")
cs("EM_sw.R"); cs("ICp2.R"); cs("Xgenerators_v190127.R")
cs("dualGRRmdA_v190215.R"); cs("CVGSBHK_v181127.R"); cs("zfun_v190304.R")
cs("factor.R"); cs("TVPRRcosso_v181120.R"); cs("TVPRR_v181111.R")
cs("fastZrot_v181125.R"); cs("CVKFMV_v190214.R")

# Verifica todas as funcoes criticas
critical_fns <- c("make_reg_matrix", "TVPRR_cosso", "Zfun",
                  "cvgs.bhk2015", "dualGRR", "hush", "garchFit")
cat("\n=== Verificacao de funcoes ===\n")
all_ok <- TRUE
for (fn in critical_fns) {
  ok <- exists(fn, mode = "function")
  cat(sprintf("  %-25s %s\n", fn, ifelse(ok, "OK", "*** FALTA ***")))
  if (!ok) all_ok <- FALSE
}
if (!all_ok) stop("Funcoes criticas faltando. Abortando.")
cat("\n")

# ============================================================
# 3. TESTE PRE-LOOP DE TVPRR_cosso [C5]
# ============================================================
cat("=== Teste pre-loop de TVPRR_cosso ===\n")
set.seed(123)
n_t <- 150; k_t <- 5
X_t <- matrix(rnorm(n_t * k_t), n_t, k_t)
b_t <- c(1, -0.5, 0.3, 0.1, -0.2)
y_t <- as.numeric(X_t %*% b_t + rnorm(n_t, 0, 0.5))

test_ok <- tryCatch({
  r_test <- TVPRR_cosso(
    X = X_t, y = y_t, type = 2,
    lambdavec = exp(seq(-2, 6, length.out = 8)),
    lambda2 = 0.1, silent = 1, kfold = 3,
    tol = 1e-3, maxit = 3
  )
  cat("  TVPRR_cosso funciona! Estrutura:\n")
  cat(sprintf("    grrats$betas_grr: %s\n",
              paste(dim(r_test$grrats$betas_grr), collapse = " x ")))
  TRUE
}, error = function(e) {
  cat(sprintf("  TVPRR_cosso FALHOU: %s\n", e$message))
  FALSE
})

if (!test_ok) {
  stop("TVPRR_cosso falha mesmo com dados sinteticos. Verifique dependencias.")
}
cat("[OK] Teste pre-loop passou\n\n")

# ============================================================
# 4. HELPERS
# ============================================================
OF <- function(pred, y, tol = 2, go.to.pred) {
  newx <- pred
  cm   <- (newx - mean(y)) > tol * (max(y) - mean(y))
  cmi  <- (newx - mean(y)) < tol * (min(y) - mean(y))
  newx[cm]  <- go.to.pred[cm]
  newx[cmi] <- go.to.pred[cmi]
  newx
}

rm_const <- function(X) {
  keep <- apply(X, 2, function(col) {
    v <- var(col, na.rm = TRUE)
    !is.na(v) && v > .Machine$double.eps
  })
  X[, keep, drop = FALSE]
}

# ============================================================
# 5. CARREGA BASE DO MEDEIROS
# ============================================================
cat("=== Carregando base de dados ===\n")
load("data/data.rda")
fred_raw <- as.data.frame(data)
bigt     <- nrow(fred_raw)

date_col <- grep("^date$", colnames(fred_raw), ignore.case = TRUE)[1]
cpi_col  <- grep("^CPIAUCSL$", colnames(fred_raw), ignore.case = TRUE)[1]
if (is.na(date_col) || is.na(cpi_col)) stop("Colunas date/CPIAUCSL nao encontradas.")

dates <- fred_raw[, date_col]
y_raw <- as.numeric(fred_raw[, cpi_col])
X_raw <- as.matrix(fred_raw[, -c(date_col, cpi_col)])
cat(sprintf("  Base: %d obs x %d preditores | %s a %s\n\n",
            bigt, ncol(X_raw), as.character(dates[1]), as.character(dates[bigt])))

# ============================================================
# 6. ALINHAMENTO tau
# ============================================================
if (file.exists("forecasts/yout.rda")) {
  load("forecasts/yout.rda")
  n_oos <- nrow(yout)
  tau   <- bigt - n_oos
  cat(sprintf("  yout.rda -> n_oos=%d, tau=%d\n", n_oos, tau))
} else {
  n_oos <- 312; tau <- bigt - n_oos
  cat(sprintf("  Fallback -> n_oos=%d, tau=%d\n", n_oos, tau))
}
if (tau < 50) stop("tau muito pequeno.")
cat(sprintf("  OOS: %s a %s\n\n", as.character(dates[tau+1]), as.character(dates[bigt])))

# ============================================================
# 7. IMPUTACAO EM (1x)
# ============================================================
cat("=== Imputacao EM ===\n")
X_imp <- tryCatch({
  em_out <- EM_sw(data = as.data.frame(X_raw), n = 8, it_max = 1000)
  cat(sprintf("  EM_sw OK | NAs restantes: %d\n", sum(is.na(em_out$data))))
  as.matrix(em_out$data)
}, error = function(e) {
  cat(sprintf("  EM_sw falhou: %s -> interpolacao\n", e$message))
  Xr <- X_raw
  for (j in seq_len(ncol(Xr))) {
    nas <- which(is.na(Xr[, j]))
    if (length(nas) > 0 && length(nas) < nrow(Xr) - 2)
      Xr[, j] <- approx(seq_len(nrow(Xr)), Xr[, j], xout = seq_len(nrow(Xr)), rule = 2)$y
  }
  Xr[is.na(Xr)] <- 0
  Xr
})
cat("\n")

# ============================================================
# 8. VARIAVEL Y ACUMULADA
# ============================================================
build_cumulative_y <- function(y, h) {
  n <- length(y); yh <- rep(NA_real_, n)
  for (t in h:n) yh[t] <- sum(y[(t - h + 1):t])
  yh
}

hor <- c(1, 3, 6, 12)
forecast_vars <- sapply(hor, build_cumulative_y, y = y_raw)
colnames(forecast_vars) <- paste0("h", hor)

# ============================================================
# 9. PARAMETROS
# ============================================================
nf  <- 8; ly <- 2; lf <- 2
lambdavec <- exp(pracma::linspace(-2, 12, n = 15))
silent <- 1

cat(sprintf("Parametros: nf=%d ly=%d lf=%d lambdas=%d\n", nf, ly, lf, length(lambdavec)))
cat(sprintf("bigt=%d tau=%d n_oos=%d\n\n", bigt, tau, n_oos))

dir.create("forecasts",   showWarnings = FALSE)
dir.create("results",     showWarnings = FALSE)
dir.create("checkpoints", showWarnings = FALSE)

# ============================================================
# 10. ARRAYS
# ============================================================
fc_ridge  <- matrix(NA_real_, nrow = bigt, ncol = length(hor))
fc_2srr   <- matrix(NA_real_, nrow = bigt, ncol = length(hor))
lam_ridge <- matrix(NA_real_, nrow = bigt, ncol = length(hor))
lam2_2srr <- matrix(NA_real_, nrow = bigt, ncol = length(hor))

betas_2srr  <- setNames(vector("list", length(hor)), paste0("h", hor))
betas_ridge <- setNames(vector("list", length(hor)), paste0("h", hor))
for (hi in seq_along(hor)) {
  betas_2srr[[hi]]  <- vector("list", n_oos)
  betas_ridge[[hi]] <- vector("list", n_oos)
}

fail_counts <- list(pca_fail=0, reg_fail=0, cv_fail=0,
                    tvp_fail=0, dim_fail=0, skip_small=0)
# Guarda mensagens de erro para debug
tvp_errors <- character(0)

# ============================================================
# 11. LOOP POOS
# ============================================================
closeAllConnections()
cat("============================================================\n")
cat("  INICIANDO LOOP POOS\n")
cat("============================================================\n\n")

t0_total <- proc.time()

for (hi in seq_along(hor)) {
  h <- hor[hi]
  cat(sprintf("\n######## HORIZONTE h = %d ########\n", h))

  n_ridge_ok <- 0; n_2srr_ok <- 0

  for (t in tau:(bigt - 1)) {
    idx   <- t - tau + 1L
    T_end <- t
    min_obs <- ly + lf + h + 30L
    if (T_end < min_obs) { fail_counts$skip_small <- fail_counts$skip_small + 1; next }

    # PCA
    X_is <- rm_const(X_imp[1:T_end, , drop = FALSE])
    if (ncol(X_is) < nf) { fail_counts$pca_fail <- fail_counts$pca_fail + 1; next }

    pc <- tryCatch(prcomp(X_is, center = TRUE, scale. = TRUE), error = function(e) NULL)
    if (is.null(pc)) { fail_counts$pca_fail <- fail_counts$pca_fail + 1; next }

    n_fac_use <- min(nf, ncol(pc$x))
    fac <- pc$x[, 1:n_fac_use, drop = FALSE]

    # Y
    y_h <- forecast_vars[1:T_end, hi]
    Y_h <- y_raw[1:T_end]
    first_valid <- which(!is.na(y_h))[1]
    if (is.na(first_valid) || first_valid >= T_end) next

    si <- first_valid
    if ((T_end - si + 1L) < min_obs) { fail_counts$skip_small <- fail_counts$skip_small + 1; next }

    y_is <- as.matrix(y_h[si:T_end])
    Y_is <- as.matrix(Y_h[si:T_end])
    f_is <- fac[si:T_end, , drop = FALSE]

    # make_reg_matrix
    reg <- tryCatch(
      make_reg_matrix(y = y_is, Y = Y_is, factors = f_is, h = h, ly = ly, lf = lf),
      error = function(e) NULL
    )
    if (is.null(reg) || nrow(reg) < (ly + lf + h + 10L)) {
      fail_counts$reg_fail <- fail_counts$reg_fail + 1; next
    }

    last <- as.numeric(reg[nrow(reg), ])
    reg  <- reg[1:(nrow(reg) - 1L), , drop = FALSE]
    ml   <- max(ly, lf)
    if (nrow(reg) <= ml + 5L) { fail_counts$reg_fail <- fail_counts$reg_fail + 1; next }
    reg <- reg[(ml + 1L):nrow(reg), , drop = FALSE]
    reg <- reg[complete.cases(reg), , drop = FALSE]
    if (nrow(reg) < 20L || ncol(reg) < 2L) { fail_counts$reg_fail <- fail_counts$reg_fail + 1; next }

    yy   <- reg[, 1]
    XX   <- reg[, -1, drop = FALSE]
    xnew <- matrix(last[-1], nrow = 1)

    if (ncol(xnew) != ncol(XX)) {
      fail_counts$dim_fail <- fail_counts$dim_fail + 1
      if (ncol(xnew) > ncol(XX)) xnew <- xnew[, 1:ncol(XX), drop = FALSE]
      else xnew <- cbind(xnew, matrix(0, 1, ncol(XX) - ncol(xnew)))
    }

    # --- RIDGE ---
    CV <- tryCatch(
      cv.glmnet(x = XX, y = yy, family = "gaussian",
                alpha = 0, nfolds = min(10, nrow(XX))),
      error = function(e) NULL
    )
    if (is.null(CV)) { fail_counts$cv_fail <- fail_counts$cv_fail + 1; next }

    mdl_r    <- glmnet(x = XX, y = yy, family = "gaussian",
                       alpha = 0, lambda = CV$lambda.min)
    pred_lin <- as.numeric(predict(mdl_r, newx = xnew))
    fc_ridge[t, hi]  <- pred_lin
    lam_ridge[t, hi] <- CV$lambda.min
    n_ridge_ok <- n_ridge_ok + 1

    rc <- as.numeric(coef(mdl_r))
    betas_ridge[[hi]][[idx]] <- list(t=t, date=dates[t], beta0=rc[1], betas=rc[-1])

    # --- 2SRR ---
    aa <- tryCatch({
      TVPRR_cosso(
        X         = XX,          # X PRIMEIRO (conforme assinatura)
        y         = yy,          # y SEGUNDO
        type      = 2,
        lambdavec = lambdavec,
        lambda2   = CV$lambda.min,
        sweigths  = 1,
        oosX      = as.numeric(xnew),
        kfold     = 5,
        silent    = silent,
        alpha     = 0.01,
        tol       = 1e-6,
        maxit     = 10
      )
    }, error = function(e) {
      # Guarda a mensagem de erro para analise
      if (length(tvp_errors) < 10) {
        tvp_errors <<- c(tvp_errors,
                         sprintf("h=%d t=%d: %s", h, t, e$message))
      }
      fail_counts$tvp_fail <<- fail_counts$tvp_fail + 1
      NULL
    })

    if (!is.null(aa)) {
      # Extrair previsao
      if (!is.null(aa$fcast) && length(aa$fcast) > 0 && !all(is.na(aa$fcast))) {
        p_raw <- as.numeric(aa$fcast[length(aa$fcast)])
      } else if (!is.null(aa$grrats) && !is.null(aa$grrats$betas_grr)) {
        bm <- aa$grrats$betas_grr
        if (is.array(bm) && length(dim(bm)) == 3) {
          # betas_grr e 3D: [1, K+1, T]
          bl <- bm[1, , dim(bm)[3]]
        } else if (is.matrix(bm)) {
          bl <- bm[nrow(bm), ]
        } else {
          bl <- as.numeric(bm)
        }
        if (length(bl) >= 1 + ncol(xnew)) {
          p_raw <- as.numeric(bl[1] + sum(bl[-1] * as.numeric(xnew)))
        } else {
          p_raw <- pred_lin
        }
      } else {
        p_raw <- pred_lin
      }

      p_filt <- OF(pred = p_raw, y = yy, go.to.pred = pred_lin)
      fc_2srr[t, hi]   <- p_filt
      lam2_2srr[t, hi] <- if (!is.null(aa$grrats$lambdas)) aa$grrats$lambdas[1] else NA_real_

      betas_2srr[[hi]][[idx]] <- list(
        t = t, date = dates[t],
        betas = if (!is.null(aa$grrats$betas_grr)) aa$grrats$betas_grr else NULL
      )
      n_2srr_ok <- n_2srr_ok + 1
    }

    # Progresso
    if (idx %% 12L == 0L) {
      el <- (proc.time() - t0_total)["elapsed"]
      rem <- el / idx * (n_oos - idx)
      oos_sf <- (tau + 1):t
      real_sf <- forecast_vars[oos_sf, hi]
      rmse_r <- sqrt(mean((fc_ridge[oos_sf, hi] - real_sf)^2, na.rm=T))
      rmse_2 <- sqrt(mean((fc_2srr[oos_sf, hi] - real_sf)^2, na.rm=T))
      ratio  <- ifelse(is.finite(rmse_r) && rmse_r > 0, rmse_2/rmse_r, NA)
      cat(sprintf("  h=%2d | %s | %3d/%d (%3.0f%%) | R:%.3f 2S:%.3f rat:%.3f | %.1fm (~%.0fm)\n",
                  h, as.character(dates[t]), idx, n_oos, 100*idx/n_oos,
                  rmse_r, rmse_2, ifelse(is.na(ratio),99,ratio), el/60, rem/60))
    }

    # Checkpoint
    if (idx %% 50L == 0L) {
      cp <- list(hi=hi, t=t, idx=idx, h=h,
                 fc_ridge=fc_ridge, fc_2srr=fc_2srr,
                 lam_ridge=lam_ridge, lam2_2srr=lam2_2srr,
                 betas_2srr=betas_2srr, betas_ridge=betas_ridge,
                 fail_counts=fail_counts, tvp_errors=tvp_errors)
      save(cp, file=sprintf("checkpoints/cp_h%d_t%d.rda", h, t))
    }
  }

  cat(sprintf("\n  h=%d CONCLUIDO | Ridge:%d/%d | 2SRR:%d/%d | %.1fmin\n",
              h, n_ridge_ok, n_oos, n_2srr_ok, n_oos,
              (proc.time()-t0_total)["elapsed"]/60))
}

el_total <- (proc.time() - t0_total)["elapsed"]
cat(sprintf("\nPOOS COMPLETO: %.1f min\n\n", el_total/60))

# Diagnostico
cat("=== Diagnostico ===\n")
for (nm in names(fail_counts)) cat(sprintf("  %-15s: %d\n", nm, fail_counts[[nm]]))
if (length(tvp_errors) > 0) {
  cat("\nPrimeiros erros do TVPRR_cosso:\n")
  for (e in tvp_errors) cat(sprintf("  %s\n", e))
}
cat("\n")

# ============================================================
# 12. SALVA RESULTADOS
# ============================================================
save(fc_ridge, fc_2srr, lam_ridge, lam2_2srr,
     file = "forecasts/coulombe_forecasts.rda")
save(betas_2srr,  file = "forecasts/coulombe_betas_2SRR.rda")
save(betas_ridge, file = "forecasts/coulombe_betas_ridge.rda")
cat("RDAs salvos\n")

# CSVs
oos_idx <- (tau + 1L):bigt
for (hi in seq_along(hor)) {
  h    <- hor[hi]
  real <- forecast_vars[oos_idx, hi]
  fr   <- fc_ridge[oos_idx, hi]
  f2   <- fc_2srr[oos_idx, hi]

  df_out <- data.frame(
    date=dates[oos_idx], realized=real,
    fc_ridge=fr, fc_2srr=f2,
    err_ridge=fr-real, err_2srr=f2-real,
    lam_ridge=lam_ridge[oos_idx,hi], lam2_2srr=lam2_2srr[oos_idx,hi])

  rmse_r <- sqrt(mean((fr-real)^2, na.rm=T))
  rmse_2 <- sqrt(mean((f2-real)^2, na.rm=T))
  fname <- sprintf("forecasts/coulombe_fc_h%02d.csv", h)
  write.csv(df_out, file=fname, row.names=FALSE)
  cat(sprintf("  h=%2d | Ridge RMSE=%.4f (%d) | 2SRR RMSE=%.4f (%d) | ratio=%.4f\n",
              h, rmse_r, sum(!is.na(fr)), rmse_2, sum(!is.na(f2)), rmse_2/rmse_r))
}

# ============================================================
# 13. EXPORTA BETAS PARA TODOS OS HORIZONTES [C4]
# ============================================================
cat("\n=== Exportando betas (todos os horizontes) ===\n")

for (hi in seq_along(hor)) {
  h <- hor[hi]

  # Ridge betas
  valid_r <- Filter(Negate(is.null), betas_ridge[[hi]])
  if (length(valid_r) > 0) {
    df_r <- do.call(rbind, lapply(valid_r, function(b) {
      c(date = as.character(b$date), t = b$t, beta0 = b$beta0,
        setNames(b$betas, paste0("b", seq_along(b$betas))))
    }))
    fname_r <- sprintf("results/betas_ridge_h%02d.csv", h)
    write.csv(as.data.frame(df_r), fname_r, row.names = FALSE)
    cat(sprintf("  Ridge h=%2d: %d janelas -> %s\n", h, nrow(df_r), fname_r))
  }

  # 2SRR betas
  valid_b <- Filter(Negate(is.null), betas_2srr[[hi]])
  if (length(valid_b) > 0) {
    df_b <- do.call(rbind, lapply(valid_b, function(b) {
      bm <- b$betas
      if (is.array(bm) && length(dim(bm)) == 3) {
        bvec <- bm[1, , dim(bm)[3]]
      } else if (is.matrix(bm)) {
        bvec <- bm[nrow(bm), ]
      } else {
        bvec <- as.numeric(bm)
      }
      c(date = as.character(b$date), t = b$t,
        setNames(bvec, paste0("b", seq_along(bvec))))
    }))
    fname_b <- sprintf("results/betas_2srr_h%02d.csv", h)
    write.csv(as.data.frame(df_b), fname_b, row.names = FALSE)
    cat(sprintf("  2SRR  h=%2d: %d janelas -> %s\n", h, nrow(df_b), fname_b))
  }
}

cat("\n============================================================\n")
cat("  06_coulombe_2SRR_pipeline.R v8.0 --- COMPLETO\n")
cat("============================================================\n")
