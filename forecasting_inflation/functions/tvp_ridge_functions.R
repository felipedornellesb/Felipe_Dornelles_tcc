# ==============================================================
# tvp_ridge_functions.R
#
# Implementação do 2SRR (Coulombe 2024) em R moderno.
# Baseado diretamente nas equações do artigo:
#   "Time-Varying Parameters as Ridge Regressions"
#   International Journal of Forecasting, 2025.
#
# Funções principais:
#   make_ZZt()                — ZZ' e Z'y sem montar Z (Eq. 4 do artigo)
#   dual_solve()              — solução dual (Eq. 9/11), numericamente estável
#   cv_ridge_dual()           — k-fold CV para selecionar lambda
#   tvp_1srr()                — Step 1: Ridge homogêneo
#   tvp_2srr()                — Step 2: 2SRR com Omega e Sigma estimados
#   run2srr()                 — wrapper no padrão rolling_window() do Medeiros
#   extract_betas_over_time() — consolida betas de todas as janelas
#   plot_betas_over_time()    — plot ggplot dos betas no tempo
#
# Dependências: glmnet (obrigatório); rugarch (opcional, fallback automático)
# ==============================================================


# ==============================================================
# 1. make_ZZt — ZZ' e Z'y sem montar Z explicitamente
# ==============================================================
# FIX 1: for (k in 1:K) → for (k in seq_len(K))
#         Em R, 1:0 itera k=1 e k=0, corrompendo ZZt/Zty quando K=0.
#         seq_len(0) retorna integer(0), ou seja, loop não executa.
#         Guard adicional: retorna matrizes zero se K==0.
# ==============================================================

make_ZZt <- function(X, y) {
  T_obs <- nrow(X)
  K     <- ncol(X)

  # Guard: K == 0 → retorna estrutura vazia (evita loop 1:0)
  if (K == 0L) {
    return(list(
      ZZt   = matrix(0, T_obs, T_obs),
      Zty   = numeric(T_obs),
      C0    = matrix(0, T_obs, T_obs),
      K     = 0L,
      T_obs = T_obs
    ))
  }

  # C0: lower triangular de 1s (random-walk cumulativa)
  C0  <- matrix(0, T_obs, T_obs)
  for (i in seq_len(T_obs)) C0[i, seq_len(i)] <- 1

  ZZt <- matrix(0, T_obs, T_obs)
  Zty <- numeric(T_obs)

  # FIX 1 aplicado: seq_len(K) em vez de 1:K
  for (k in seq_len(K)) {
    xk  <- X[, k]
    Zk  <- outer(xk, rep(1, T_obs)) * C0   # T×T: Zk[t,s] = X[t,k]*C0[t,s]
    ZZt <- ZZt + tcrossprod(Zk)
    Zty <- Zty + Zk %*% y
  }

  list(ZZt = ZZt, Zty = Zty, C0 = C0, K = K, T_obs = T_obs)
}


# ==============================================================
# 2. dual_solve — solução (ZZ' + λI)^{-1} Z'y
# ==============================================================

dual_solve <- function(ZZt, Zty, lam, eps = 1e-8) {
  T_obs <- nrow(ZZt)
  M     <- ZZt + (lam + eps) * diag(T_obs)
  alpha <- tryCatch(
    solve(M, Zty),
    error = function(e) {
      # Aumenta regularização se ainda singular
      tryCatch(
        solve(M + 1e-4 * diag(T_obs), Zty),
        error = function(e2) rep(0, T_obs)   # último recurso: zeros
      )
    }
  )
  alpha
}


# ==============================================================
# 3. recover_beta — paths β(t) a partir da solução dual
# ==============================================================
# FIX aplicado aqui também: seq_len(K) no loop interno

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


#' OOS forecast: x_out · beta(T)
oos_forecast_from_beta <- function(x_out, beta) {
  sum(x_out * beta[nrow(beta), ])
}


# ==============================================================
# 4. cv_ridge_dual — K-fold CV para selecionar lambda
# ==============================================================

cv_ridge_dual <- function(X, y,
                           lambdas = exp(seq(-4, 20, length.out = 30)),
                           kfold   = 5) {
  T_obs  <- nrow(X)
  folds  <- cut(seq_len(T_obs), breaks = kfold, labels = FALSE)
  folds  <- sample(folds)

  cv_mse <- numeric(length(lambdas))

  for (li in seq_along(lambdas)) {
    lam    <- lambdas[li]
    errors <- numeric(T_obs)

    for (f in seq_len(kfold)) {
      test_idx  <- which(folds == f)
      train_idx <- which(folds != f)

      X_tr <- X[train_idx, , drop = FALSE]
      y_tr <- y[train_idx]
      X_te <- X[test_idx,  , drop = FALSE]
      y_te <- y[test_idx]

      zz     <- make_ZZt(X_tr, y_tr)
      alpha  <- dual_solve(zz$ZZt, zz$Zty, lam)
      beta   <- recover_beta(X_tr, alpha, zz$C0)
      beta_T <- beta[nrow(beta), ]

      y_hat            <- X_te %*% beta_T
      errors[test_idx] <- y_te - y_hat
    }
    cv_mse[li] <- mean(errors^2, na.rm = TRUE)
  }

  lambdas[which.min(cv_mse)]
}


# ==============================================================
# 5. tvp_1srr — Step 1: Ridge TVP homogêneo (Eq. 9)
# ==============================================================

tvp_1srr <- function(X, y, kfold = 5,
                      lambdas = exp(seq(-4, 20, length.out = 30))) {
  lam_opt <- cv_ridge_dual(X, y, lambdas = lambdas, kfold = kfold)
  zz      <- make_ZZt(X, y)
  alpha   <- dual_solve(zz$ZZt, zz$Zty, lam_opt)
  beta    <- recover_beta(X, alpha, zz$C0)
  resid   <- y - rowSums(X * beta)

  list(beta = beta, resid = resid, lambda = lam_opt,
       alpha = alpha, ZZt = zz$ZZt, C0 = zz$C0)
}


# ==============================================================
# 6. estimate_sigma2 — variâncias dos resíduos (Passo 2 do artigo)
# ==============================================================
# FIX 2: fallback robusto sem rugarch.
#   Problema original: var(x) com window de 1 elemento retorna NA,
#   propagando NAs para toda a série sigma2, quebrando o GLS.
#   Solução: usar soma de quadrados manual (sempre >= 0, nunca NA)
#   com janela expansiva mínima de 12 obs.
# ==============================================================

estimate_sigma2 <- function(resid) {
  T_obs  <- length(resid)
  sigma2 <- rep(1, T_obs)
  fit_ok <- FALSE

  # Tenta GARCH(1,1) via rugarch se disponível
  if (requireNamespace("rugarch", quietly = TRUE)) {
    tryCatch({
      spec   <- rugarch::ugarchspec(
        variance.model     = list(model = "sGARCH", garchOrder = c(1, 1)),
        mean.model         = list(armaOrder = c(0, 0), include.mean = TRUE),
        distribution.model = "norm"
      )
      fit    <- rugarch::ugarchfit(spec = spec, data = resid,
                                   solver = "hybrid",
                                   solver.control = list(trace = 0))
      s2     <- as.numeric(rugarch::sigma(fit))^2
      if (all(is.finite(s2)) && all(s2 > 0)) {
        sigma2 <- s2
        fit_ok <- TRUE
      }
    }, error = function(e) NULL)
  }

  # FIX 2: fallback robusto — variância móvel via soma de quadrados
  # Usa window mínima de 12 obs (expansiva nas primeiras linhas).
  # sum(x^2)/n - mean(x)^2  ≥ 0 sempre e nunca produz NA para n ≥ 1.
  if (!fit_ok) {
    win <- 12L
    for (t in seq_len(T_obs)) {
      idx <- max(1L, t - win + 1L):t
      r   <- resid[idx]
      n   <- length(r)
      if (n < 2L) {
        # Para as primeiras obs: usa variância global como prior
        sigma2[t] <- var(resid, na.rm = TRUE)
      } else {
        sigma2[t] <- sum((r - mean(r))^2) / (n - 1L)
      }
    }
    # Garante positividade estrita
    sigma2 <- pmax(sigma2, 1e-8)
  }

  # Normaliza: mean(sigma2) = 1  (conforme artigo, Seção 3.2)
  mu_s2 <- mean(sigma2, na.rm = TRUE)
  if (!is.finite(mu_s2) || mu_s2 <= 0) return(rep(1, T_obs))
  sigma2 / mu_s2
}


# ==============================================================
# 7. estimate_omega — variâncias dos parâmetros (Passo 3 do artigo)
# ==============================================================

estimate_omega <- function(beta) {
  dbeta <- diff(beta)            # (T-1) × K
  omega <- colMeans(dbeta^2)
  omega <- pmax(omega, 1e-12)
  omega / mean(omega)
}


# ==============================================================
# 8. tvp_2srr — Step 2: GLS ponderado (Eq. 11 do artigo)
# ==============================================================
# X̃[t,k] = X[t,k] * sqrt(omega_k) / sigma_t
# ỹ[t]    = y[t]   / sigma_t
# Roda ridge homogêneo em X̃, ỹ e desescala betas.

tvp_2srr <- function(X, y, kfold = 5,
                      lambdas = exp(seq(-4, 20, length.out = 30))) {

  step1 <- tvp_1srr(X, y, kfold = kfold, lambdas = lambdas)

  sigma2     <- estimate_sigma2(step1$resid)
  omega      <- estimate_omega(step1$beta)
  inv_sigma  <- 1 / sqrt(pmax(sigma2, 1e-8))
  sqrt_omega <- sqrt(pmax(omega, 1e-12))

  X_tilde <- sweep(X, 2, sqrt_omega, "*")
  X_tilde <- sweep(X_tilde, 1, inv_sigma, "*")
  y_tilde <- y * inv_sigma

  lam_opt2   <- cv_ridge_dual(X_tilde, y_tilde,
                               lambdas = lambdas, kfold = kfold)
  zz2        <- make_ZZt(X_tilde, y_tilde)
  alpha2     <- dual_solve(zz2$ZZt, zz2$Zty, lam_opt2)
  beta_tilde <- recover_beta(X_tilde, alpha2, zz2$C0)
  beta_orig  <- sweep(beta_tilde, 2, sqrt_omega, "*")
  resid2     <- y - rowSums(X * beta_orig)

  list(
    beta         = beta_orig,
    resid        = resid2,
    lambda       = lam_opt2,
    omega        = omega,
    sigma2       = sigma2,
    beta_step1   = step1$beta,
    lambda_step1 = step1$lambda
  )
}


# ==============================================================
# 9. run2srr — wrapper padrão Medeiros
# ==============================================================
# FIX 3: guard explícito is.matrix(Xin) + tryCatch total envolvendo
#   todo o corpo da função, evitando que erros internos disparem o
#   modo interativo Browse[1]> no rolling_window() do Medeiros.
# ==============================================================

run2srr <- function(ind, df, variable, horizon,
                    kfold   = 5,
                    lambdas = exp(seq(-4, 20, length.out = 25))) {

  # FIX 3: tryCatch externo — qualquer erro retorna NA sem travar o loop
  result <- tryCatch({

    # 1. dataprep idêntico ao Medeiros (nofact=TRUE: sem PCA, com dummy)
    prep <- dataprep(ind, df, variable, horizon, nofact = TRUE)
    Xin  <- prep$Xin
    yin  <- prep$yin
    Xout <- as.numeric(prep$Xout)

    # FIX 3: guards explícitos de tipo e dimensão
    if (!is.matrix(Xin))
      Xin <- as.matrix(Xin)
    if (nrow(Xin) == 0 || ncol(Xin) == 0)
      return(list(forecast = NA_real_, outputs = NULL))
    if (length(yin) < 20 || ncol(Xin) < 2)
      return(list(forecast = NA_real_, outputs = NULL))

    # 2. Remove colunas com variância zero
    col_var   <- apply(Xin, 2, var, na.rm = TRUE)
    good_cols <- which(col_var > 1e-10)
    if (length(good_cols) < 2)
      return(list(forecast = NA_real_, outputs = NULL))
    Xin  <- Xin[, good_cols, drop = FALSE]
    Xout <- Xout[good_cols]

    # 3. Padroniza X (estabilidade numérica do dual ridge)
    X_means <- colMeans(Xin)
    X_sds   <- apply(Xin, 2, sd)
    X_sds[X_sds < 1e-10] <- 1
    Xin_sc  <- sweep(sweep(Xin, 2, X_means, "-"), 2, X_sds, "/")
    Xout_sc <- (Xout - X_means) / X_sds

    # 4. Estima 2SRR
    fit <- tvp_2srr(Xin_sc, yin, kfold = kfold, lambdas = lambdas)

    # 5. Previsão OOS
    beta_T <- fit$beta[nrow(fit$beta), ]
    fcast  <- sum(Xout_sc * beta_T)

    # Fallback Ridge se previsão implausível
    y_mu <- mean(yin)
    y_sd <- sd(yin)
    if (!is.finite(fcast) || abs(fcast - y_mu) > 5 * y_sd) {
      fcast <- tryCatch({
        cv_r <- glmnet::cv.glmnet(Xin_sc, yin, alpha = 0, nfolds = kfold)
        as.numeric(predict(
          glmnet::glmnet(Xin_sc, yin, alpha = 0, lambda = cv_r$lambda.min),
          newx = matrix(Xout_sc, nrow = 1)
        ))
      }, error = function(e) y_mu)
    }

    # 6. Despadroniza betas para escala original
    beta_orig <- sweep(fit$beta, 2, X_sds[good_cols], "/")

    # 7. Outputs no formato Medeiros
    outputs <- list(
      betas_time_varying = beta_orig,
      lambda             = fit$lambda,
      omega              = fit$omega,
      sigma2             = fit$sigma2,
      n_obs              = nrow(Xin),
      n_vars             = ncol(Xin)
    )

    list(forecast = fcast, outputs = outputs)

  }, error = function(e) {
    # Qualquer erro inesperado: loga e retorna NA sem interromper o loop
    message(sprintf("  [run2srr h=%d] ERRO: %s", horizon, conditionMessage(e)))
    list(forecast = NA_real_, outputs = NULL)
  })

  result
}


# ==============================================================
# 10. extract_betas_over_time — consolida betas de todas as janelas
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


# ==============================================================
# 11. plot_betas_over_time — top-N betas com maior variância
# ==============================================================

plot_betas_over_time <- function(df_betas, var_names = NULL,
                                  top_n = 10, variable = "CPIAUCSL",
                                  save_path = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Instale ggplot2: install.packages('ggplot2')")

  beta_var <- tapply(df_betas$beta, df_betas$var_idx, var, na.rm = TRUE)
  top_vars <- order(beta_var, decreasing = TRUE)[seq_len(min(top_n,
                                                              length(beta_var)))]
  df_p <- df_betas[df_betas$var_idx %in% top_vars, ]
  df_p$var_label <- if (!is.null(var_names) &&
                          max(df_p$var_idx) <= length(var_names)) {
    var_names[df_p$var_idx]
  } else {
    paste0("X", df_p$var_idx)
  }

  p <- ggplot2::ggplot(
    df_p,
    ggplot2::aes(x = window_date, y = beta,
                 colour = var_label, group = var_label)
  ) +
    ggplot2::geom_line(linewidth = 0.7, alpha = 0.85) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        colour = "grey60", linewidth = 0.4) +
    ggplot2::scale_x_date(date_labels = "%Y", date_breaks = "2 years") +
    ggplot2::labs(
      title    = sprintf("Betas Time-Varying — 2SRR | %s", variable),
      subtitle = sprintf("Top %d variáveis por variância dos betas", top_n),
      x        = NULL,
      y        = "Coeficiente \u03b2(t)",
      colour   = "Variável"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      legend.position    = "right",
      panel.grid.minor   = ggplot2::element_blank()
    )

  if (!is.null(save_path))
    ggplot2::ggsave(save_path, p, width = 11, height = 5, dpi = 150)

  p
}
