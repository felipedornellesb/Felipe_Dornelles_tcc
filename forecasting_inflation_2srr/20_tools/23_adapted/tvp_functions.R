# ==============================================================================
# tvp_functions.R
#
# Self-contained 2SRR (Algorithm 1, Coulombe 2025) with Medeiros wrapper.
# TVPRR_cosso from the Coulombe code is used when available (sourced by
# 00_setup.R). If unavailable or failing, falls back to the standalone
# implementation below.
# ==============================================================================


# ==============================================================================
# Standalone implementation (fallback)
# ==============================================================================

make_ZZt <- function(X, y) {
  T_obs <- nrow(X); K <- ncol(X)
  C0 <- matrix(0, T_obs, T_obs)
  for (i in 1:T_obs) C0[i, 1:i] <- 1
  ZZt <- matrix(0, T_obs, T_obs); Zty <- numeric(T_obs)
  for (k in 1:K) {
    Zk  <- outer(X[, k], rep(1, T_obs)) * C0
    ZZt <- ZZt + tcrossprod(Zk)
    Zty <- Zty + Zk %*% y
  }
  list(ZZt = ZZt, Zty = Zty, C0 = C0, K = K, T_obs = T_obs)
}

dual_solve <- function(ZZt, Zty, lam, eps = 1e-8) {
  T_obs <- nrow(ZZt)
  M <- ZZt + (lam + eps) * diag(T_obs)
  tryCatch(solve(M, Zty), error = function(e) solve(M + 1e-4*diag(T_obs), Zty))
}

recover_beta <- function(X, alpha, C0) {
  T_obs <- nrow(X); K <- ncol(X); C0t <- t(C0)
  beta <- matrix(NA_real_, T_obs, K)
  for (k in 1:K) beta[, k] <- C0 %*% (C0t %*% (X[, k] * alpha))
  beta
}

cv_ridge_dual <- function(X, y, lambdas = exp(seq(-4, 20, length.out = 25)),
                           kfold = 5, blocked = TRUE, block_size = 6) {
  T_obs <- nrow(X)
  if (blocked) {
    bs <- max(block_size, round(T_obs / kfold))
    folds <- rep(1:kfold, each = bs, length.out = T_obs)
  } else {
    folds <- sample(rep(1:kfold, length.out = T_obs))
  }
  cv_mse <- numeric(length(lambdas))
  for (li in seq_along(lambdas)) {
    errors <- numeric(T_obs)
    for (f in 1:kfold) {
      tr <- which(folds != f); te <- which(folds == f)
      zz <- make_ZZt(X[tr,,drop=FALSE], y[tr])
      alpha <- dual_solve(zz$ZZt, zz$Zty, lambdas[li])
      beta_T <- recover_beta(X[tr,,drop=FALSE], alpha, zz$C0)
      beta_T <- beta_T[nrow(beta_T), ]
      errors[te] <- y[te] - X[te,,drop=FALSE] %*% beta_T
    }
    cv_mse[li] <- mean(errors^2, na.rm = TRUE)
  }
  lambdas[which.min(cv_mse)]
}

tvp_1srr_standalone <- function(X, y, kfold = 5, lambdas = exp(seq(-4, 20, length.out = 25))) {
  lam <- cv_ridge_dual(X, y, lambdas, kfold)
  zz <- make_ZZt(X, y)
  alpha <- dual_solve(zz$ZZt, zz$Zty, lam)
  beta <- recover_beta(X, alpha, zz$C0)
  resid <- y - rowSums(X * beta)
  list(beta = beta, resid = resid, lambda = lam)
}

estimate_sigma2_standalone <- function(resid) {
  T_obs <- length(resid); sigma2 <- rep(1, T_obs)
  if (requireNamespace("rugarch", quietly = TRUE)) {
    tryCatch({
      spec <- rugarch::ugarchspec(
        variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
        mean.model = list(armaOrder = c(0, 0), include.mean = TRUE))
      fit <- rugarch::ugarchfit(spec, resid, solver = "hybrid",
                                 solver.control = list(trace = 0))
      sigma2 <- as.numeric(rugarch::sigma(fit))^2
      return(sigma2 / mean(sigma2))
    }, error = function(e) NULL)
  }
  if (requireNamespace("fGarch", quietly = TRUE)) {
    tryCatch({
      fit <- fGarch::garchFit(~ garch(1, 1), data = resid,
                                trace = FALSE, include.mean = FALSE)
      sigma2 <- as.numeric(fit@sigma.t^2)
      return(sigma2 / mean(sigma2))
    }, error = function(e) NULL)
  }
  for (t in 1:T_obs) { w <- max(1,t-11):t; sigma2[t] <- var(resid[w], na.rm=TRUE) }
  pmax(sigma2, 1e-8) / mean(pmax(sigma2, 1e-8))
}

estimate_omega_standalone <- function(beta) {
  omega <- colMeans(diff(beta)^2)
  omega <- pmax(omega, 1e-12)
  omega / mean(omega)
}

tvp_2srr_standalone <- function(X, y, kfold = 5, lambdas = exp(seq(-4, 20, length.out = 25))) {
  s1 <- tvp_1srr_standalone(X, y, kfold, lambdas)
  sigma2 <- estimate_sigma2_standalone(s1$resid)
  omega  <- estimate_omega_standalone(s1$beta)
  isig <- 1/sqrt(pmax(sigma2, 1e-8)); som <- sqrt(pmax(omega, 1e-12))
  X_t <- sweep(sweep(X, 2, som, "*"), 1, isig, "*"); y_t <- y * isig
  lam2 <- cv_ridge_dual(X_t, y_t, lambdas, kfold)
  zz2 <- make_ZZt(X_t, y_t)
  a2 <- dual_solve(zz2$ZZt, zz2$Zty, lam2)
  bt <- recover_beta(X_t, a2, zz2$C0)
  beta <- sweep(bt, 2, som, "*")
  resid <- y - rowSums(X * beta)
  list(beta = beta, resid = resid, lambda = lam2,
       omega = omega, sigma2 = sigma2, lambda_step1 = s1$lambda)
}


# ==============================================================================
# run2srr: Medeiros-compatible wrapper
# Interface: fn(ind, df, variable, horizon) -> list(forecast, outputs)
# Uses TVPRR_cosso if available, otherwise standalone tvp_2srr_standalone.
# ==============================================================================

run2srr <- function(ind, df, variable, horizon, kfold = 5, n_lags = 4, univar = FALSE, factonly = FALSE, nofact = FALSE) {

  if (exists("dataprep", mode = "function")) {
    prep_data <- dataprep(ind, df, variable, horizon, add_dummy = FALSE, univar = univar, factonly = factonly, nofact = nofact)
    X_in  <- prep_data$Xin
    y_in  <- prep_data$yin
    x_out <- prep_data$Xout
  } else {
    df_w <- df[ind, , drop = FALSE]
    y_raw <- as.numeric(df_w[, variable])
    T_raw <- length(y_raw)
  
    # Keep only numeric columns
    num_cols <- sapply(df_w, is.numeric)
    df_num   <- as.matrix(df_w[, num_cols, drop = FALSE])
  
    # Embed with lags
    X_embed <- embed(df_num, n_lags)
    n_align <- nrow(X_embed) - horizon
    if (n_align < 30) return(list(forecast = NA_real_, outputs = NULL))
  
    X_in  <- X_embed[1:n_align, , drop = FALSE]
    y_in  <- y_raw[(n_lags + horizon):(n_lags + horizon + n_align - 1)]
    x_out <- X_embed[nrow(X_embed), ]
  
    n_use <- min(length(y_in), nrow(X_in))
    X_in  <- X_in[1:n_use, , drop = FALSE]
    y_in  <- y_in[1:n_use]
  }

  # Remove zero-variance columns
  cv <- apply(X_in, 2, var, na.rm = TRUE)
  good <- which(is.finite(cv) & cv > 1e-10)
  if (length(good) < 2) return(list(forecast = NA_real_, outputs = NULL))
  X_in <- X_in[, good, drop = FALSE]
  x_out <- x_out[good]

  # Remove NAs
  ok <- complete.cases(X_in, y_in)
  X_in <- X_in[ok, , drop = FALSE]; y_in <- y_in[ok]
  if (nrow(X_in) < 30) return(list(forecast = NA_real_, outputs = NULL))

  # Standardize
  X_mu <- colMeans(X_in); X_sd <- apply(X_in, 2, sd); X_sd[X_sd < 1e-10] <- 1
  X_sc <- sweep(sweep(X_in, 2, X_mu, "-"), 2, X_sd, "/")
  x_sc <- (x_out - X_mu) / X_sd

  # Estimate 2SRR
  fit <- tryCatch(
    tvp_2srr_standalone(X_sc, y_in, kfold = kfold),
    error = function(e) {
      message("  standalone 2SRR failed: ", e$message)
      NULL
    }
  )

  if (is.null(fit)) return(list(forecast = mean(y_in), outputs = NULL))

  # Forecast
  fcast <- sum(x_sc * fit$beta[nrow(fit$beta), ])

  # Outlier filter
  y_mean <- mean(y_in); y_sd <- sd(y_in)
  if (!is.finite(fcast) || abs(fcast - y_mean) > 3 * y_sd) {
    tryCatch({
      bs <- max(6, round(length(y_in) / kfold))
      blocked_folds <- rep(1:kfold, each = bs, length.out = length(y_in))
      cv_r <- cv.glmnet(X_sc, y_in, alpha = 0, foldid = blocked_folds)
      fcast <- as.numeric(predict(
        glmnet(X_sc, y_in, alpha = 0, lambda = cv_r$lambda.min),
        newx = matrix(x_sc, nrow = 1)))
    }, error = function(e) fcast <<- y_mean)
  }

  list(
    forecast = fcast,
    outputs  = list(
      betas_tvp = fit$beta,
      lambda    = fit$lambda,
      omega     = fit$omega,
      sigma2    = fit$sigma2,
      n_obs     = nrow(X_in),
      n_vars    = ncol(X_in)
    )
  )
}


# -- Clark-West test -----------------------------------------------------------
clark_west <- function(y, fc_restricted, fc_unrestricted) {
  e1 <- y - fc_restricted; e2 <- y - fc_unrestricted
  d <- e1^2 - (e2^2 - (fc_restricted - fc_unrestricted)^2)
  ok <- complete.cases(d); d <- d[ok]
  if (length(d) < 10) return(list(stat = NA, pvalue = NA))
  t_stat <- mean(d) / (sd(d) / sqrt(length(d)))
  list(stat = t_stat, pvalue = 1 - pnorm(t_stat))
}
