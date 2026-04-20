# ==============================================================
# felipe_functions.R
#
# Implementação do 2SRR (Coulombe 2024) em R.
# Baseado em: "Time-Varying Parameters as Ridge Regressions"
#   International Journal of Forecasting, 2025.
#
# CHANGELOG v3 (final):
#   FIX-1: seq_len(K) em todos os loops (1:0 é perigoso)
#   FIX-2: estimate_sigma2 robusto sem rugarch
#   FIX-3: tryCatch externo em run2srr()
#   FIX-4: cv_ridge_dual() usa folds TEMPORAIS (sem embaralhar)
#   FIX-5: make_ZZt() usa C0 * xk sem outer() (memória eficiente)
#   FIX-6: dummy preservada em run2srr() — não removida pelo good_cols
#   FIX-7: plot_betas_over_time() usa "\u03b2" correto (1 barra)
#
# Dependências: glmnet (obrigatório); rugarch (opcional)
# ==============================================================


# ==============================================================
# 1. make_ZZt — ZZ' e Z'y sem montar Z explicitamente
# ==============================================================
# FIX-1: seq_len(K) em vez de 1:K
# FIX-5: usa C0 * xk (broadcast column-wise) em vez de outer()
#   C0 * xk equivale a sweep(C0, 1, xk, "*")
#   É numericamente idêntico mas evita alocar T×T temporário K vezes.

make_ZZt <- function(X, y) {
  T_obs <- nrow(X)
  K     <- ncol(X)

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
  for (i in seq_len(T_obs)) C0[i, seq_len(i)] <- 1L

  ZZt <- matrix(0, T_obs, T_obs)
  Zty <- numeric(T_obs)

  # FIX-5: C0 * xk multiplica cada LINHA t de C0 por X[t,k]
  # (R recicla xk como vetor coluna quando aplicado sobre matriz por linhas)
  # Equivalente a sweep(C0, 1, xk, "*") mas sem overhead de sweep()
  for (k in seq_len(K)) {
    Zk  <- C0 * X[, k]           # T×T: Zk[t,s] = C0[t,s] * X[t,k]
    ZZt <- ZZt + crossprod(Zk)   # t(Zk) %*% Zk — equivale a tcrossprod(t(Zk))
    Zty <- Zty + t(Zk) %*% y
  }

  list(ZZt = ZZt, Zty = as.numeric(Zty), C0 = C0, K = K, T_obs = T_obs)
}


# ==============================================================
# 2. dual_solve — (ZZ' + λI)^{-1} Z'y
# ==============================================================

dual_solve <- function(ZZt, Zty, lam, eps = 1e-8) {
  T_obs <- nrow(ZZt)
  M     <- ZZt + (lam + eps) * diag(T_obs)
  alpha <- tryCatch(
    solve(M, Zty),
    error = function(e) tryCatch(
      solve(M + 1e-4 * diag(T_obs), Zty),
      error = function(e2) rep(0, T_obs)
    )
  )
  alpha
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


# (auxiliar, exportado para uso externo se necessário)
oos_forecast_from_beta <- function(x_out, beta) {
  sum(x_out * beta[nrow(beta), ])
}


# ==============================================================
# 4. cv_ridge_dual — K-fold CV TEMPORAL para selecionar lambda
# ==============================================================
# FIX-4: folds temporais em vez de sample(folds).
#   Em séries temporais, o fold de teste deve sempre ser POSTERIOR
#   ao fold de treino. Usa "blocked time-series CV":
#     - divide a série em kfold blocos contíguos por ordem temporal
#     - para fold f: treina em 1:(inicio_f - 1), testa em bloco_f
#   Isso evita vazamento de dados futuros para o treino.
#   Ref: Bergmeir & Benítez (2012), Hyndman & Athanasopoulos (2021).
# ==============================================================

cv_ridge_dual <- function(X, y,
                           lambdas = exp(seq(-4, 20, length.out = 30)),
                           kfold   = 5) {
  T_obs <- nrow(X)

  # FIX-4: folds temporais contíguos (NÃO embaralhados)
  # Bloco mínimo: pelo menos 10 obs de treino antes do primeiro fold
  min_train <- max(10L, floor(T_obs * 0.3))
  usable    <- T_obs - min_train          # obs disponíveis para testar
  fold_size <- max(1L, floor(usable / kfold))

  cv_mse <- numeric(length(lambdas))

  for (li in seq_along(lambdas)) {
    lam    <- lambdas[li]
    sq_err <- numeric(kfold)
    n_err  <- integer(kfold)

    for (f in seq_len(kfold)) {
      test_end   <- min_train + f * fold_size
      test_start <- min_train + (f - 1L) * fold_size + 1L
      if (test_start > T_obs) break
      test_end   <- min(test_end, T_obs)

      train_idx <- seq_len(test_start - 1L)
      test_idx  <- test_start:test_end

      X_tr <- X[train_idx, , drop = FALSE]
      y_tr <- y[train_idx]
      X_te <- X[test_idx,  , drop = FALSE]
      y_te <- y[test_idx]

      if (length(train_idx) < 5L || length(test_idx) < 1L) next

      zz     <- make_ZZt(X_tr, y_tr)
      alpha  <- dual_solve(zz$ZZt, zz$Zty, lam)
      beta   <- recover_beta(X_tr, alpha, zz$C0)
      beta_T <- beta[nrow(beta), ]     # betas do fim do treino

      y_hat      <- X_te %*% beta_T
      sq_err[f]  <- sum((y_te - y_hat)^2, na.rm = TRUE)
      n_err[f]   <- length(y_te)
    }

    total_n <- sum(n_err)
    cv_mse[li] <- if (total_n > 0) sum(sq_err) / total_n else Inf
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

  list(beta    = beta,
       resid   = resid,
       lambda  = lam_opt,
       alpha   = alpha,
       ZZt     = zz$ZZt,
       C0      = zz$C0)
}


# ==============================================================
# 6. estimate_sigma2 — variâncias dos resíduos (Passo 2)
# ==============================================================
# FIX-2: fallback robusto via soma de quadrados manual.
#   var() com n=1 retorna NA; sum((r - mean(r))^2)/(n-1) com n<2
#   é tratado explicitamente com prior = var global.

estimate_sigma2 <- function(resid) {
  T_obs   <- length(resid)
  sigma2  <- rep(1, T_obs)
  fit_ok  <- FALSE

  if (requireNamespace("rugarch", quietly = TRUE)) {
    tryCatch({
      spec <- rugarch::ugarchspec(
        variance.model     = list(model = "sGARCH", garchOrder = c(1, 1)),
        mean.model         = list(armaOrder = c(0, 0), include.mean = TRUE),
        distribution.model = "norm"
      )
      fit  <- rugarch::ugarchfit(spec  = spec, data = resid,
                                  solver = "hybrid",
                                  solver.control = list(trace = 0))
      s2   <- as.numeric(rugarch::sigma(fit))^2
      if (all(is.finite(s2)) && all(s2 > 0)) {
        sigma2 <- s2
        fit_ok <- TRUE
      }
    }, error = function(e) NULL)
  }

  # FIX-2: fallback — variância móvel expansiva (janela = 12)
  if (!fit_ok) {
    var_global <- var(resid, na.rm = TRUE)
    if (!is.finite(var_global) || var_global <= 0) var_global <- 1
    win <- 12L
    for (t in seq_len(T_obs)) {
      idx <- max(1L, t - win + 1L):t
      r   <- resid[idx]
      n   <- length(r)
      sigma2[t] <- if (n < 2L) var_global
                   else sum((r - mean(r))^2) / (n - 1L)
    }
    sigma2 <- pmax(sigma2, 1e-8)
  }

  mu_s2 <- mean(sigma2, na.rm = TRUE)
  if (!is.finite(mu_s2) || mu_s2 <= 0) return(rep(1, T_obs))
  sigma2 / mu_s2
}


# ==============================================================
# 7. estimate_omega — variâncias dos parâmetros (Passo 3)
# ==============================================================

estimate_omega <- function(beta) {
  dbeta <- diff(beta)
  omega <- colMeans(dbeta^2)
  omega <- pmax(omega, 1e-12)
  omega / mean(omega)
}


# ==============================================================
# 8. tvp_2srr — Step 2: GLS ponderado (Eq. 11)
# ==============================================================

tvp_2srr <- function(X, y, kfold = 5,
                      lambdas = exp(seq(-4, 20, length.out = 30))) {

  step1      <- tvp_1srr(X, y, kfold = kfold, lambdas = lambdas)
  sigma2     <- estimate_sigma2(step1$resid)
  omega      <- estimate_omega(step1$beta)
  inv_sigma  <- 1 / sqrt(pmax(sigma2, 1e-8))
  sqrt_omega <- sqrt(pmax(omega, 1e-12))

  X_tilde <- sweep(X,       2, sqrt_omega, "*")
  X_tilde <- sweep(X_tilde, 1, inv_sigma,  "*")
  y_tilde <- y * inv_sigma

  lam_opt2   <- cv_ridge_dual(X_tilde, y_tilde,
                               lambdas = lambdas, kfold = kfold)
  zz2        <- make_ZZt(X_tilde, y_tilde)
  alpha2     <- dual_solve(zz2$ZZt, zz2$Zty, lam_opt2)
  beta_tilde <- recover_beta(X_tilde, alpha2, zz2$C0)
  beta_orig  <- sweep(beta_tilde, 2, sqrt_omega, "*")
  resid2     <- y - rowSums(X * beta_orig)

  list(
    beta          = beta_orig,
    resid         = resid2,
    lambda        = lam_opt2,
    omega         = omega,
    sigma2        = sigma2,
    beta_step1    = step1$beta,
    lambda_step1  = step1$lambda
  )
}


# ==============================================================
# 9. run2srr — wrapper padrão Medeiros
# ==============================================================
# FIX-3: tryCatch externo envolve tudo
# FIX-6: dummy é preservada explicitamente — não pode cair no
#   filtro good_cols porque sua variância pode ser 0 (sem crise
#   na janela). A dummy é separada de Xin antes do filtro e
#   reincorporada depois, garantindo que Xout também bata.

run2srr <- function(ind, df, variable, horizon,
                    kfold   = 5,
                    lambdas = exp(seq(-4, 20, length.out = 25))) {

  result <- tryCatch({

    # 1. dataprep idêntico ao Medeiros (nofact=TRUE, com dummy)
    prep <- dataprep(ind, df, variable, horizon, nofact = TRUE)
    Xin  <- prep$Xin
    yin  <- prep$yin
    Xout <- as.numeric(prep$Xout)

    if (!is.matrix(Xin)) Xin <- as.matrix(Xin)
    if (nrow(Xin) == 0 || ncol(Xin) == 0 || length(yin) < 20)
      return(list(forecast = NA_real_, outputs = NULL))

    # FIX-6: separa a dummy (última coluna) antes de filtrar por variância
    # dataprep() sempre adiciona dummy como última coluna (add_dummy=TRUE)
    n_cols     <- ncol(Xin)
    dummy_col  <- Xin[, n_cols, drop = FALSE]   # salva dummy
    dummy_out  <- Xout[n_cols]                   # escalar OOS da dummy
    Xin_main   <- Xin[, -n_cols, drop = FALSE]  # remove dummy do filtro
    Xout_main  <- Xout[-n_cols]

    # 2. Filtra variância zero apenas nas colunas principais
    col_var   <- apply(Xin_main, 2, var, na.rm = TRUE)
    good_cols <- which(col_var > 1e-10)
    if (length(good_cols) < 2)
      return(list(forecast = NA_real_, outputs = NULL))

    Xin_main  <- Xin_main[, good_cols, drop = FALSE]
    Xout_main <- Xout_main[good_cols]

    # Reincorpora dummy
    Xin_f  <- cbind(Xin_main, dummy_col)
    Xout_f <- c(Xout_main, dummy_out)

    # 3. Padroniza X (exceto dummy — já está em 0/1)
    n_main  <- length(good_cols)
    X_means <- colMeans(Xin_f[, seq_len(n_main), drop = FALSE])
    X_sds   <- apply(Xin_f[, seq_len(n_main), drop = FALSE], 2, sd)
    X_sds[X_sds < 1e-10] <- 1

    Xin_sc  <- Xin_f
    Xout_sc <- Xout_f
    Xin_sc[, seq_len(n_main)]  <- sweep(
      sweep(Xin_f[, seq_len(n_main), drop = FALSE], 2, X_means, "-"),
      2, X_sds, "/"
    )
    Xout_sc[seq_len(n_main)] <- (Xout_f[seq_len(n_main)] - X_means) / X_sds

    # 4. Estima 2SRR
    fit <- tvp_2srr(Xin_sc, yin, kfold = kfold, lambdas = lambdas)

    # 5. Previsão OOS com betas do último período
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
    # (apenas colunas principais; dummy permanece em escala 0/1)
    beta_orig <- fit$beta
    beta_orig[, seq_len(n_main)] <- sweep(
      fit$beta[, seq_len(n_main), drop = FALSE],
      2, X_sds, "/"
    )

    outputs <- list(
      betas_time_varying = beta_orig,
      lambda             = fit$lambda,
      omega              = fit$omega,
      sigma2             = fit$sigma2,
      n_obs              = nrow(Xin_f),
      n_vars             = ncol(Xin_f)
    )

    list(forecast = fcast, outputs = outputs)

  }, error = function(e) {
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
# FIX-7: "\u03b2" com UMA barra — R interpreta \u em tempo de parse.
#   Com DUAS barras (\\u03b2) o R trata como string literal e imprime
#   o texto "\u03b2" no eixo em vez do caractere β.

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
      title    = sprintf("Betas Time-Varying \u2014 2SRR | %s", variable),
      subtitle = sprintf("Top %d vari\u00e1veis por vari\u00e2ncia dos betas", top_n),
      x        = NULL,
      y        = "Coeficiente \u03b2(t)",   # FIX-7: 1 barra
      colour   = "Vari\u00e1vel"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      legend.position  = "right",
      panel.grid.minor = ggplot2::element_blank()
    )

  if (!is.null(save_path))
    ggplot2::ggsave(save_path, p, width = 11, height = 5, dpi = 150)

  p
}