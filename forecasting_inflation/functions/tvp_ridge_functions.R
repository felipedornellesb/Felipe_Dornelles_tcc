# ==============================================================
# v2: tvp_ridge_functions.R
#
# Implementação do 2SRR (Coulombe 2024) em R.
# Padrão Medeiros (ForecastingInflation) — compatível com
#   rolling_window(), dataprep(), accumulate_model().
#
# Baseado em: "Time-Varying Parameters as Ridge Regressions"
#   International Journal of Forecasting, 2025.
# Autor: Felipe Dornelles
# ==============================================================
# 1. make_ZZt — constrói ZZ' e Z'y (chamado UMA vez com lambda ótimo)
# ==============================================================

make_ZZt <- function(X, y) {
  T_obs <- nrow(X)
  K     <- ncol(X)

  C0 <- matrix(0, T_obs, T_obs)
  for (i in seq_len(T_obs)) C0[i, seq_len(i)] <- 1L

  ZZt <- matrix(0, T_obs, T_obs)
  Zty <- numeric(T_obs)

  for (k in seq_len(K)) {
    Zk  <- C0 * X[, k]
    ZZt <- ZZt + crossprod(Zk)
    Zty <- Zty + t(Zk) %*% y
  }

  list(ZZt = ZZt, Zty = as.numeric(Zty), C0 = C0)
}


# ==============================================================
# 2. dual_solve — resolve (ZZ' + λI)^{-1} Z'y
# ==============================================================

dual_solve <- function(ZZt, Zty, lam, eps = 1e-8) {
  T_obs <- nrow(ZZt)
  M     <- ZZt + (lam + eps) * diag(T_obs)
  tryCatch(
    solve(M, Zty),
    error = function(e) solve(M + 1e-4 * diag(T_obs), Zty)
  )
}


# ==============================================================
# 3. recover_beta — paths β(t) a partir da solução dual
# ==============================================================

recover_beta <- function(X, alpha, C0) {
  T_obs <- nrow(X)
  K     <- ncol(X)
  beta  <- matrix(NA_real_, T_obs, K)
  C0t   <- t(C0)

  for (k in seq_len(K)) {
    Zkt_alpha <- C0t %*% (X[, k] * alpha)
    beta[, k] <- C0 %*% Zkt_alpha
  }
  beta
}


# ==============================================================
# 4. tvp_2srr_fit — estima o modelo dado X, y, lambda
#   lambda selecionado externamente via glmnet (rápido)
# ==============================================================

tvp_2srr_fit <- function(X, y, lam) {

  # Passo 1: Ridge TVP homogêneo
  zz1    <- make_ZZt(X, y)
  alpha1 <- dual_solve(zz1$ZZt, zz1$Zty, lam)
  beta1  <- recover_beta(X, alpha1, zz1$C0)
  resid1 <- y - rowSums(X * beta1)

  # Variância dos resíduos (janela móvel 12 meses — sem GARCH)
  T_obs      <- length(resid1)
  var_global <- var(resid1, na.rm = TRUE)
  if (!is.finite(var_global) || var_global <= 0) var_global <- 1
  sigma2 <- rep(var_global, T_obs)
  win    <- 12L
  for (t in seq_len(T_obs)) {
    idx <- max(1L, t - win + 1L):t
    r   <- resid1[idx]
    n   <- length(r)
    if (n >= 2L) {
      v <- sum((r - mean(r))^2) / (n - 1L)
      if (is.finite(v) && v > 0) sigma2[t] <- v
    }
  }
  sigma2 <- pmax(sigma2, 1e-8)
  sigma2 <- sigma2 / mean(sigma2)

  # Variância dos parâmetros
  dbeta <- diff(beta1)
  omega <- colMeans(dbeta^2)
  omega <- pmax(omega, 1e-12)
  omega <- omega / mean(omega)

  # Passo 2: GLS ponderado
  inv_sig    <- 1 / sqrt(pmax(sigma2, 1e-8))
  sqrt_omega <- sqrt(pmax(omega, 1e-12))

  X_tilde <- sweep(sweep(X,       2, sqrt_omega, "*"), 1, inv_sig, "*")
  y_tilde <- y * inv_sig

  zz2        <- make_ZZt(X_tilde, y_tilde)
  alpha2     <- dual_solve(zz2$ZZt, zz2$Zty, lam)
  beta_tilde <- recover_beta(X_tilde, alpha2, zz2$C0)
  beta_orig  <- sweep(beta_tilde, 2, sqrt_omega, "/")

  list(
    beta         = beta_orig,
    beta_step1   = beta1,
    omega        = omega,
    sigma2       = sigma2,
    lambda       = lam
  )
}


# ==============================================================
# 5. run2srr — wrapper padrão Medeiros
#
# Usa glmnet::cv.glmnet para selecionar lambda (igual ao runlasso
# do Medeiros) — rápido. make_ZZt só é chamado 2x (passo 1 e 2).
# ==============================================================

# Original do Coulombe:
#run2srr <- function(ind, df, variable, horizon,
#                    K_pca_max = 40L,    # mais fatores PCA
#                    var_expl  = 0.90,   # 90% variância explicada
#                    kfold     = 5) {    # 5-fold CV (padrão do paper)

run2srr <- function(ind, df, variable, horizon,
                    K_pca_max = 20L,
                    var_expl  = 0.85,
                    kfold     = 3) {

  result <- tryCatch({

    # 1. dataprep idêntico ao Medeiros
    prep <- dataprep(ind, df, variable, horizon, nofact = TRUE)
    Xin  <- prep$Xin
    yin  <- prep$yin
    Xout <- as.numeric(prep$Xout)

    if (!is.matrix(Xin)) Xin <- as.matrix(Xin)
    if (nrow(Xin) < 20 || ncol(Xin) == 0)
      return(list(forecast = NA_real_, outputs = NULL))

    # 2. Separa dummy (última coluna — add_dummy=TRUE no dataprep)
    nc        <- ncol(Xin)
    dummy_col <- Xin[, nc, drop = FALSE]
    dummy_out <- Xout[nc]
    Xin       <- Xin[, -nc, drop = FALSE]
    Xout      <- Xout[-nc]

    # 3. Remove colunas com variância zero
    col_var <- apply(Xin, 2, var, na.rm = TRUE)
    ok      <- col_var > 1e-10
    Xin     <- Xin[, ok, drop = FALSE]
    Xout    <- Xout[ok]
    if (ncol(Xin) < 2)
      return(list(forecast = NA_real_, outputs = NULL))

    # 4. Padroniza
    X_mu  <- colMeans(Xin)
    X_sd  <- apply(Xin, 2, sd)
    X_sd[X_sd < 1e-10] <- 1
    Xin_s  <- sweep(sweep(Xin, 2, X_mu, "-"), 2, X_sd, "/")
    Xout_s <- (Xout - X_mu) / X_sd

    # 5. PCA — reduz dimensão antes de make_ZZt
    pca      <- prcomp(Xin_s, center = FALSE, scale. = FALSE)
    cum_var  <- cumsum(pca$sdev^2) / sum(pca$sdev^2)
    K_pca    <- max(5L, min(
      K_pca_max,
      which(cum_var >= var_expl)[1L],
      nrow(Xin_s) - kfold - 1L
    ))

    Xin_pca  <- pca$x[, 1:K_pca, drop = FALSE]
    Xout_pca <- as.numeric(matrix(Xout_s, nrow = 1) %*%
                             pca$rotation[, 1:K_pca])

    # Reincorpora dummy
    Xin_f  <- cbind(Xin_pca, dummy_col)
    Xout_f <- c(Xout_pca, dummy_out)

    # 6. Seleciona lambda via glmnet (rápido — igual ao Medeiros)
    cv_fit <- glmnet::cv.glmnet(Xin_f, yin,
                                 alpha  = 0,
                                 nfolds = kfold)
    lam_opt <- cv_fit$lambda.min

    # 7. Estima 2SRR com lambda ótimo (make_ZZt chamado só 2x)
    fit    <- tvp_2srr_fit(Xin_f, yin, lam_opt)
    beta_T <- fit$beta[nrow(fit$beta), ]
    fcast  <- sum(Xout_f * beta_T)

    # Fallback: se previsão implausível, usa predição Ridge do glmnet
    y_mu <- mean(yin)
    y_sd <- sd(yin)
    if (!is.finite(fcast) || abs(fcast - y_mu) > 5 * y_sd) {
      fcast <- as.numeric(predict(
        glmnet::glmnet(Xin_f, yin, alpha = 0, lambda = lam_opt),
        newx = matrix(Xout_f, nrow = 1)
      ))
    }

    outputs <- list(
      betas_time_varying = fit$beta,
      betas_step1        = fit$beta_step1,
      pca_rotation       = pca$rotation[, 1:K_pca],
      pca_center         = X_mu,
      pca_scale          = X_sd,
      K_pca              = K_pca,
      lambda             = lam_opt,
      omega              = fit$omega,
      sigma2             = fit$sigma2
    )

    list(forecast = fcast, outputs = outputs)

  }, error = function(e) {
    message(sprintf("  [run2srr h=%d] ERRO: %s", horizon, conditionMessage(e)))
    list(forecast = NA_real_, outputs = NULL)
  })

  result
}


# ==============================================================
# 6. extract_betas_over_time
# ==============================================================

extract_betas_over_time <- function(model_list, df, nwindows) {
  n_windows  <- length(model_list$outputs)
  dates_all  <- as.Date(rownames(df))
  win_dates  <- tail(dates_all, n_windows)
  betas_list <- vector("list", n_windows)

  for (i in seq_len(n_windows)) {
    out <- model_list$outputs[[i]]
    if (is.null(out) || is.null(out$betas_time_varying)) next
    b      <- out$betas_time_varying
    last_b <- if (is.matrix(b)) as.numeric(b[nrow(b), ]) else as.numeric(b)
    betas_list[[i]] <- data.frame(
      window_date = win_dates[i],
      var_idx     = seq_along(last_b),
      beta        = last_b
    )
  }
  do.call(rbind, betas_list)
}