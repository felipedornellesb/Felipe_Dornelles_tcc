# ============================================================
# 08_tvp_comparison.R
#
# Compara 3 especificacoes TVP (2SRR):
#   Caso 1: TVP-AR     — so lags de y (univariado)
#   Caso 2: TVP-Factor — so fatores PCA
#   Caso 3: TVP-FAVAR  — fatores + lags de y (JA RODOU no 06)
#
# + Ridge constante como benchmark em cada caso
# ============================================================

rm(list = ls())
gc()

setwd("~/TCC/tcc/forecasting_inflation")
cat("  08_tvp_comparison.R — Comparacao de especificacoes TVP\n")

# ============================================================
# 0. PACOTES
# ============================================================
pkgs <- c("pracma", "glmnet", "matrixcalc", "GA", "e1071",
          "fGarch", "timeSeries")
invisible(lapply(pkgs, library, character.only = TRUE))
cat("[OK] Pacotes carregados\n")

# ============================================================
# 1. FUNCOES COULOMBE
# ============================================================
base_factor_backup <- base::factor
factor <- function(X, n_fac = NULL, ...) {
  if (!is.null(n_fac)) {
    X <- as.matrix(X)
    Tobs <- nrow(X)
    S <- (1 / Tobs) * t(X) %*% X
    eig <- eigen(S, symmetric = TRUE)
    nf <- min(n_fac, ncol(X))
    lam <- eig$vectors[, 1:nf, drop = FALSE]
    fac <- X %*% lam
    fit <- fac %*% t(lam)
    mse <- mean((X - fit)^2, na.rm = TRUE)
    return(list(factors = fac, lambda = lam, mse = mse))
  }
  base_factor_backup(X, ...)
}

cs <- function(f) {
  p <- file.path("coulombe", f)
  if (!file.exists(p)) { warning(paste0("Nao encontrado: ", p)); return(invisible(NULL)) }
  source(p, local = FALSE)
}

cat("Carregando funcoes Coulombe...\n")
cs("EM_sw.R"); cs("ICp2.R"); cs("Xgenerators_v190127.R")
cs("dualGRRmdA_v190215.R"); cs("CVGSBHK_v181127.R"); cs("zfun_v190304.R")
cs("factor.R"); cs("TVPRRcosso_v181120.R"); cs("TVPRR_v181111.R")
cs("fastZrot_v181125.R"); cs("CVKFMV_v190214.R")
cat("[OK] Funcoes carregadas\n\n")

# ============================================================
# 2. HELPERS
# ============================================================
OF <- function(pred, y, tol = 2, go.to.pred) {
  newx <- pred
  cm  <- (newx - mean(y)) > tol * (max(y) - mean(y))
  cmi <- (newx - mean(y)) < tol * (min(y) - mean(y))
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

# Funcao para construir matriz de regressao para TVP-AR (so lags de y)
make_ar_matrix <- function(y_cum, y_level, h, ly) {
  n <- length(y_cum)
  if (n < ly + h + 5) return(NULL)

  # Colunas: [y_cum_t | y_level_{t-1} | y_level_{t-2} | ...]
  max_lag <- ly
  start <- max_lag + 1
  end_row <- n

  rows <- list()
  for (t in start:end_row) {
    row <- y_cum[t]  # dependente
    for (l in 1:ly) {
      row <- c(row, y_level[t - l])  # lags de y em nivel
    }
    rows[[length(rows) + 1]] <- row
  }

  mat <- do.call(rbind, rows)
  colnames(mat) <- c("y", paste0("y_lag", 1:ly))
  mat
}

# Funcao para construir matriz com so fatores (sem lags de y)
make_factor_matrix <- function(y_cum, factors, h, lf) {
  n <- nrow(factors)
  if (length(y_cum) != n) {
    min_n <- min(length(y_cum), n)
    y_cum <- y_cum[1:min_n]
    factors <- factors[1:min_n, , drop = FALSE]
    n <- min_n
  }

  nf <- ncol(factors)
  max_lag <- lf
  start <- max_lag + 1

  rows <- list()
  for (t in start:n) {
    row <- y_cum[t]  # dependente
    for (l in 1:lf) {
      row <- c(row, factors[t - l, ])  # lags dos fatores
    }
    rows[[length(rows) + 1]] <- row
  }

  mat <- do.call(rbind, rows)
  col_names <- "y"
  for (l in 1:lf) {
    col_names <- c(col_names, paste0("f", 1:nf, "_lag", l))
  }
  colnames(mat) <- col_names
  mat
}

# ============================================================
# 3. CARREGA DADOS
# ============================================================
cat("=== Carregando dados ===\n")
load("data/data.rda")
fred_raw <- as.data.frame(data)
bigt <- nrow(fred_raw)

date_col <- grep("^date$", colnames(fred_raw), ignore.case = TRUE)[1]
cpi_col  <- grep("^CPIAUCSL$", colnames(fred_raw), ignore.case = TRUE)[1]
dates <- fred_raw[, date_col]
y_raw <- as.numeric(fred_raw[, cpi_col])
X_raw <- as.matrix(fred_raw[, -c(date_col, cpi_col)])

load("forecasts/yout.rda")
n_oos <- nrow(yout)
tau   <- bigt - n_oos

cat(sprintf("  Base: %d obs | tau=%d | n_oos=%d\n", bigt, tau, n_oos))

# Imputacao EM (1x)
cat("  Imputacao EM...\n")
X_imp <- tryCatch({
  em_out <- EM_sw(data = as.data.frame(X_raw), n = 8, it_max = 1000)
  as.matrix(em_out$data)
}, error = function(e) {
  Xr <- X_raw
  for (j in seq_len(ncol(Xr))) {
    nas <- which(is.na(Xr[, j]))
    if (length(nas) > 0 && length(nas) < nrow(Xr) - 2)
      Xr[, j] <- approx(seq_len(nrow(Xr)), Xr[, j],
                         xout = seq_len(nrow(Xr)), rule = 2)$y
  }
  Xr[is.na(Xr)] <- 0; Xr
})

# Y acumulado
build_cumulative_y <- function(y, h) {
  n <- length(y); yh <- rep(NA_real_, n)
  for (t in h:n) yh[t] <- sum(y[(t - h + 1):t])
  yh
}

hor <- c(1, 3, 6, 12)
forecast_vars <- sapply(hor, build_cumulative_y, y = y_raw)
cat("[OK] Dados prontos\n\n")

# ============================================================
# 4. PARAMETROS
# ============================================================
nf <- 8; ly <- 2; lf <- 2
lambdavec <- exp(pracma::linspace(-2, 12, n = 15))
silent <- 1

# Definicao dos casos
cases <- list(
  list(name = "TVP_AR",     univar = TRUE,  factonly = FALSE, nofact = TRUE),
  list(name = "TVP_Factor", univar = FALSE, factonly = TRUE,  nofact = FALSE)
  # Caso 3 (TVP_FAVAR) ja rodou no script 06 — carregaremos os resultados
)

cat(sprintf("Casos a rodar: %d (+ Caso 3 ja rodado)\n", length(cases)))
cat(sprintf("Horizontes: %s\n\n", paste(hor, collapse = ", ")))

dir.create("forecasts", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("checkpoints", showWarnings = FALSE)

# ============================================================
# 5. LOOP PRINCIPAL POR CASO
# ============================================================

for (ci in seq_along(cases)) {
  caso <- cases[[ci]]
  case_name <- caso$name

  cat(sprintf("  CASO: %s\n", case_name))
  cat(sprintf("  univar=%s | factonly=%s | nofact=%s\n",
              caso$univar, caso$factonly, caso$nofact))

  # Arrays de resultado para este caso
  fc_ridge_c <- matrix(NA_real_, nrow = bigt, ncol = length(hor))
  fc_2srr_c  <- matrix(NA_real_, nrow = bigt, ncol = length(hor))
  betas_2srr_c  <- setNames(vector("list", length(hor)), paste0("h", hor))
  betas_ridge_c <- setNames(vector("list", length(hor)), paste0("h", hor))
  for (hi in seq_along(hor)) {
    betas_2srr_c[[hi]]  <- vector("list", n_oos)
    betas_ridge_c[[hi]] <- vector("list", n_oos)
  }

  fail_c <- list(skip=0, reg=0, cv=0, tvp=0)
  tvp_err_c <- character(0)

  t0 <- proc.time()

  for (hi in seq_along(hor)) {
    h <- hor[hi]
    cat(sprintf("\n  ### %s h=%d ###\n", case_name, h))

    n_ok_r <- 0; n_ok_2 <- 0

    for (t_idx in tau:(bigt - 1)) {
      idx   <- t_idx - tau + 1L
      T_end <- t_idx
      min_obs <- ly + lf + h + 30L
      if (T_end < min_obs) { fail_c$skip <- fail_c$skip + 1; next }

      y_h <- forecast_vars[1:T_end, hi]
      Y_h <- y_raw[1:T_end]
      first_valid <- which(!is.na(y_h))[1]
      if (is.na(first_valid) || first_valid >= T_end) next
      si <- first_valid
      if ((T_end - si + 1L) < min_obs) { fail_c$skip <- fail_c$skip + 1; next }

      y_cum_is <- y_h[si:T_end]
      y_lev_is <- Y_h[si:T_end]

      # --------------------------------------------------------
      # Construir matriz de regressao conforme o caso
      # --------------------------------------------------------
      reg <- NULL

      if (caso$univar && caso$nofact) {
        # CASO 1: TVP-AR — so lags de y
        reg <- tryCatch(
          make_ar_matrix(y_cum = y_cum_is, y_level = y_lev_is,
                         h = h, ly = ly),
          error = function(e) NULL)

      } else if (!caso$univar && caso$factonly) {
        # CASO 2: TVP-Factor — so fatores PCA
        X_is <- rm_const(X_imp[si:T_end, , drop = FALSE])
        if (ncol(X_is) < nf) { fail_c$skip <- fail_c$skip + 1; next }
        pc <- tryCatch(prcomp(X_is, center = TRUE, scale. = TRUE),
                       error = function(e) NULL)
        if (is.null(pc)) { fail_c$skip <- fail_c$skip + 1; next }
        n_fac_use <- min(nf, ncol(pc$x))
        fac_is <- pc$x[, 1:n_fac_use, drop = FALSE]

        reg <- tryCatch(
          make_factor_matrix(y_cum = y_cum_is, factors = fac_is,
                             h = h, lf = lf),
          error = function(e) NULL)
      }

      if (is.null(reg)) { fail_c$reg <- fail_c$reg + 1; next }

      # Limpar reg
      reg <- reg[complete.cases(reg), , drop = FALSE]
      if (nrow(reg) < 20 || ncol(reg) < 2) { fail_c$reg <- fail_c$reg + 1; next }

      # Ultima linha para previsao
      last <- as.numeric(reg[nrow(reg), ])
      reg  <- reg[1:(nrow(reg) - 1L), , drop = FALSE]
      if (nrow(reg) < 20) { fail_c$reg <- fail_c$reg + 1; next }

      yy   <- reg[, 1]
      XX   <- reg[, -1, drop = FALSE]
      xnew <- matrix(last[-1], nrow = 1)

      if (ncol(xnew) != ncol(XX)) {
        if (ncol(xnew) > ncol(XX)) xnew <- xnew[, 1:ncol(XX), drop = FALSE]
        else xnew <- cbind(xnew, matrix(0, 1, ncol(XX) - ncol(xnew)))
      }

      # --- RIDGE ---
      CV <- tryCatch(
        cv.glmnet(x = XX, y = yy, family = "gaussian",
                  alpha = 0, nfolds = min(10, nrow(XX))),
        error = function(e) NULL)
      if (is.null(CV)) { fail_c$cv <- fail_c$cv + 1; next }

      mdl_r <- glmnet(x = XX, y = yy, family = "gaussian",
                      alpha = 0, lambda = CV$lambda.min)
      pred_lin <- as.numeric(predict(mdl_r, newx = xnew))
      fc_ridge_c[t_idx, hi] <- pred_lin
      n_ok_r <- n_ok_r + 1

      rc <- as.numeric(coef(mdl_r))
      betas_ridge_c[[hi]][[idx]] <- list(t=t_idx, date=dates[t_idx],
                                         beta0=rc[1], betas=rc[-1])

      # --- 2SRR ---
      aa <- tryCatch({
        TVPRR_cosso(
          X = XX, y = yy, type = 2,
          lambdavec = lambdavec,
          lambda2 = CV$lambda.min,
          sweigths = 1,
          oosX = as.numeric(xnew),
          kfold = 5, silent = silent,
          alpha = 0.01, tol = 1e-6, maxit = 10
        )
      }, error = function(e) {
        if (length(tvp_err_c) < 5)
          tvp_err_c <<- c(tvp_err_c, sprintf("h=%d t=%d: %s", h, t_idx, e$message))
        fail_c$tvp <<- fail_c$tvp + 1
        NULL
      })

      if (!is.null(aa)) {
        if (!is.null(aa$fcast) && length(aa$fcast) > 0 && !all(is.na(aa$fcast))) {
          p_raw <- as.numeric(aa$fcast[length(aa$fcast)])
        } else if (!is.null(aa$grrats) && !is.null(aa$grrats$betas_grr)) {
          bm <- aa$grrats$betas_grr
          if (is.array(bm) && length(dim(bm)) == 3) bl <- bm[1,,dim(bm)[3]]
          else if (is.matrix(bm)) bl <- bm[nrow(bm),]
          else bl <- as.numeric(bm)
          if (length(bl) >= 1 + ncol(xnew))
            p_raw <- as.numeric(bl[1] + sum(bl[-1] * as.numeric(xnew)))
          else p_raw <- pred_lin
        } else {
          p_raw <- pred_lin
        }

        p_filt <- OF(pred = p_raw, y = yy, go.to.pred = pred_lin)
        fc_2srr_c[t_idx, hi] <- p_filt
        n_ok_2 <- n_ok_2 + 1

        betas_2srr_c[[hi]][[idx]] <- list(
          t = t_idx, date = dates[t_idx],
          betas = if (!is.null(aa$grrats$betas_grr)) aa$grrats$betas_grr else NULL)
      }

      # Progresso
      if (idx %% 24L == 0L) {
        el <- (proc.time() - t0)["elapsed"]
        rem <- el / idx * (n_oos - idx)
        oos_sf <- (tau + 1):t_idx
        real_sf <- forecast_vars[oos_sf, hi]
        rmse_r <- sqrt(mean((fc_ridge_c[oos_sf, hi] - real_sf)^2, na.rm=T))
        rmse_2 <- sqrt(mean((fc_2srr_c[oos_sf, hi] - real_sf)^2, na.rm=T))
        ratio  <- ifelse(is.finite(rmse_r) && rmse_r > 0, rmse_2/rmse_r, NA)
        cat(sprintf("    %s h=%d | %3d/%d (%2.0f%%) | R:%.3f 2S:%.3f r:%.3f | %.0fm\n",
                    case_name, h, idx, n_oos, 100*idx/n_oos,
                    rmse_r, rmse_2, ifelse(is.na(ratio),99,ratio), el/60))
      }

      # Checkpoint
      if (idx %% 100L == 0L) {
        cp <- list(case=case_name, hi=hi, t=t_idx, idx=idx,
                   fc_ridge=fc_ridge_c, fc_2srr=fc_2srr_c,
                   fail=fail_c, tvp_err=tvp_err_c)
        save(cp, file=sprintf("checkpoints/cp_%s_h%d_t%d.rda", case_name, h, t_idx))
      }
    }

    cat(sprintf("    %s h=%d DONE | Ridge:%d/%d | 2SRR:%d/%d\n",
                case_name, h, n_ok_r, n_oos, n_ok_2, n_oos))
  }

  el_case <- (proc.time() - t0)["elapsed"]
  cat(sprintf("\n  %s COMPLETO: %.1f min\n", case_name, el_case/60))
  cat("  Falhas:", paste(names(fail_c), fail_c, sep="=", collapse=" | "), "\n")
  if (length(tvp_err_c) > 0) {
    cat("  Erros TVP:\n")
    for (e in tvp_err_c) cat(sprintf("    %s\n", e))
  }

  # Salva resultados deste caso
  save(fc_ridge_c, fc_2srr_c, betas_2srr_c, betas_ridge_c,
       file = sprintf("forecasts/tvp_%s_forecasts.rda", case_name))

  # CSVs
  oos_idx <- (tau + 1L):bigt
  for (hi in seq_along(hor)) {
    h <- hor[hi]
    real <- forecast_vars[oos_idx, hi]
    fr <- fc_ridge_c[oos_idx, hi]
    f2 <- fc_2srr_c[oos_idx, hi]

    df_out <- data.frame(
      date=dates[oos_idx], realized=real,
      fc_ridge=fr, fc_2srr=f2,
      err_ridge=fr-real, err_2srr=f2-real)

    fname <- sprintf("forecasts/tvp_%s_h%02d.csv", case_name, h)
    write.csv(df_out, file=fname, row.names=FALSE)

    rmse_r <- sqrt(mean((fr-real)^2, na.rm=T))
    rmse_2 <- sqrt(mean((f2-real)^2, na.rm=T))
    cat(sprintf("  %s h=%2d | Ridge=%.4f | 2SRR=%.4f | ratio=%.4f\n",
                case_name, h, rmse_r, rmse_2, rmse_2/rmse_r))
  }
  cat("\n")
}

cat("  08_tvp_comparison.R --- COMPLETO\n")
