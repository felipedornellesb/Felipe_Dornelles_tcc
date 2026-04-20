# ==============================================================
# tvp_ridge_functions.R
#
# Implementação do 2SRR (Coulombe 2024) em R moderno.
# Baseado diretamente nas equações do artigo:
#   "Time-Varying Parameters as Ridge Regressions"
#   International Journal of Forecasting, 2025.
#
# Funções principais:
#   make_Z_matrix()   — constrói a matriz Z = W*C (Eq. 4 do artigo)
#   dual_ridge()      — solução dual (Eq. 9/11), numericamente estável
#   cv_ridge_dual()   — k-fold CV para selecionar lambda
#   tvp_1srr()        — Step 1: Ridge homogêneo
#   tvp_2srr()        — Step 2: 2SRR com Omega e Sigma estimados
#   run2srr()         — wrapper no padrão rolling_window() do Medeiros
#   extract_betas_over_time() — consolida betas de todas as janelas
#   plot_betas_over_time()    — plot ggplot dos betas no tempo
#
# Dependências: glmnet, Matrix, rugarch (ou fGarch como fallback)
# ==============================================================


# ==============================================================
# 1. Construção da matriz Z (reparametrização do artigo, Seção 2.2)
# ==============================================================
# Z = W * C
# W = diagblock(X1, ..., XK)  — T x KT
# C = I_K ⊗ C0               — KT x KT
# C0 = lower triangular de 1s (cumulativa, random walk)
# Resultado: Z é T x KT, estrutura em blocos (Eq. Z no artigo)
#
# Para o dual, só precisamos de Z Z' (T x T), que é computável
# como sum_k (X_k * C0) * (X_k * C0)' de forma eficiente.
# ==============================================================

#' Calcula Z %*% t(Z) e Z %*% y diretamente (sem montar Z explicitamente)
#' @param X  matrix T x K — regressores no tempo
#' @param y  vector T     — variável dependente
#' @return list(ZZt = T×T, Zty = T×1, C0 = T×T)
make_ZZt <- function(X, y) {
  T_obs <- nrow(X)
  K     <- ncol(X)

  # Matriz cumulativa C0 (lower triangular de 1s) — random walk
  C0  <- matrix(0, T_obs, T_obs)
  for (i in 1:T_obs) C0[i, 1:i] <- 1

  # ZZ' = sum_{k=1}^{K} (X_k * C0) (X_k * C0)'
  # Onde X_k * C0: cada linha t de Z_k é X[t,k] * C0[t,]
  # (X_k * C0)' * (X_k * C0) contribui para ZZ' como outer product escalado
  ZZt <- matrix(0, T_obs, T_obs)
  Zty <- numeric(T_obs)

  for (k in 1:K) {
    xk   <- X[, k]               # vetor T
    Zk   <- outer(xk, rep(1, T_obs)) * C0  # T x T: Zk[t,s] = X[t,k]*C0[t,s]
    ZZt  <- ZZt + tcrossprod(Zk)            # T x T
    Zty  <- Zty + Zk %*% y
  }

  list(ZZt = ZZt, Zty = Zty, C0 = C0, K = K, T_obs = T_obs)
}


# ==============================================================
# 2. Solução dual da ridge (Eq. 9 e 11 do artigo)
# ==============================================================
# Solução: θ = C * Z' * (Z*Z' + λ*I_T)^{-1} * y
# Para previsão fora da amostra usamos apenas α = (ZZ'+λI)^{-1} * Zy
# e depois fcast = Zout * C * α equivalente via:
#   fcast = Zout' * α  onde Zout é o vetor de regressores OOS expandido
#
# Versão com matrizes de peso Omega (param variances) e Sigma (resid var):
# θ = C * Omega^{1/2} * Z̃' * (Z̃*Z̃' + λ*I)^{-1} * ỹ
# onde Z̃ = Sigma^{-1/2} * Z * Omega^{1/2}, ỹ = Sigma^{-1/2} * y
# ==============================================================

#' Solução dual ridge (numericamente estável via solve com ridge numérico)
#' @param ZZt   T×T matrix Z*Z'
#' @param Zty   T×1 vector Z'y
#' @param lam   escalar lambda (penalidade)
#' @param eps   regularização numérica adicional (default 1e-8)
#' @return alpha — vetor T (dual solution)
dual_solve <- function(ZZt, Zty, lam, eps = 1e-8) {
  T_obs <- nrow(ZZt)
  # (ZZ' + λ*I + eps*I)^{-1} * Z'y
  # eps protege contra singularidade numérica
  M     <- ZZt + (lam + eps) * diag(T_obs)
  alpha <- tryCatch(
    solve(M, Zty),
    error = function(e) {
      # fallback: aumenta regularização se ainda singular
      solve(M + 1e-4 * diag(T_obs), Zty)
    }
  )
  alpha
}


#' Recover beta paths (T x K matrix) from dual solution alpha
#' β = C * Z' * α  = C * (sum_k diag(X_k) * C0)' * alpha
#' Para cada k: β_k = C0 * diag(X_k)' * alpha = cumsum(X_k * (C0' * alpha))
#' @param X     T×K matrix
#' @param alpha T×1 dual solution
#' @param C0    T×T lower triangular cumulative matrix
#' @return beta T×K matrix of time-varying coefficients
recover_beta <- function(X, alpha, C0) {
  T_obs <- nrow(X)
  K     <- ncol(X)
  beta  <- matrix(NA_real_, T_obs, K)
  C0t   <- t(C0)  # T×T

  for (k in 1:K) {
    # Z_k' * alpha = (diag(X[,k]) * C0)' * alpha = C0' * (X[,k] * alpha)
    Zkt_alpha  <- C0t %*% (X[, k] * alpha)  # T×1
    # beta_k = C0 * Zkt_alpha (random walk: cumsum)
    beta[, k]  <- C0 %*% Zkt_alpha
  }
  beta
}


#' OOS forecast via dual solution
#' fcast = x_out * beta_T = x_out * (últimos coeficientes)
#' @param x_out  K×1 vetor de regressores OOS
#' @param beta   T×K matrix (apenas a última linha é usada)
oos_forecast_from_beta <- function(x_out, beta) {
  sum(x_out * beta[nrow(beta), ])
}


# ==============================================================
# 3. K-fold Cross Validation para selecionar lambda
# ==============================================================

#' K-fold CV para o TVP ridge (dual)
#' @param X       T×K matrix de regressores
#' @param y       T×1 vetor dependente
#' @param lambdas vetor de lambdas candidatos
#' @param kfold   número de folds (default 5)
#' @return lambda ótimo (scalar)
cv_ridge_dual <- function(X, y, lambdas = exp(seq(-4, 20, length.out = 30)),
                           kfold = 5) {
  T_obs  <- nrow(X)
  # divide em folds aleatórios (sem respeitar ordem temporal — OK para k-fold
  # em séries com lags suficientes, cf. Bergmeir et al. 2018 citado no artigo)
  folds  <- cut(seq_len(T_obs), breaks = kfold, labels = FALSE)
  folds  <- sample(folds)  # embaralha

  cv_mse <- numeric(length(lambdas))

  for (li in seq_along(lambdas)) {
    lam    <- lambdas[li]
    errors <- numeric(T_obs)

    for (f in 1:kfold) {
      test_idx  <- which(folds == f)
      train_idx <- which(folds != f)

      X_tr  <- X[train_idx, , drop = FALSE]
      y_tr  <- y[train_idx]
      X_te  <- X[test_idx,  , drop = FALSE]
      y_te  <- y[test_idx]

      # Computa ZZt e Zty no conjunto de treino
      zz    <- make_ZZt(X_tr, y_tr)
      alpha <- dual_solve(zz$ZZt, zz$Zty, lam)
      beta  <- recover_beta(X_tr, alpha, zz$C0)
      beta_T <- beta[nrow(beta), ]  # coefs no fim do treino

      # Previsão para cada obs de teste (usa beta fixo do treino)
      y_hat <- X_te %*% beta_T
      errors[test_idx] <- y_te - y_hat
    }
    cv_mse[li] <- mean(errors^2, na.rm = TRUE)
  }

  lambdas[which.min(cv_mse)]
}


# ==============================================================
# 4. Step 1 do 2SRR: TVP ridge homogêneo (Eq. 9)
# ==============================================================

tvp_1srr <- function(X, y, kfold = 5,
                      lambdas = exp(seq(-4, 20, length.out = 30))) {
  # Seleciona lambda por CV
  lam_opt <- cv_ridge_dual(X, y, lambdas = lambdas, kfold = kfold)

  # Estimação final com lambda ótimo
  zz      <- make_ZZt(X, y)
  alpha   <- dual_solve(zz$ZZt, zz$Zty, lam_opt)
  beta    <- recover_beta(X, alpha, zz$C0)
  resid   <- y - rowSums(X * beta)

  list(beta = beta, resid = resid, lambda = lam_opt,
       alpha = alpha, ZZt = zz$ZZt, C0 = zz$C0)
}


# ==============================================================
# 5. Estimação de Sigma e Omega (passos 2 e 3 do Algorithm 1)
# ==============================================================

#' Estima variâncias de resíduos via GARCH(1,1) — passo 2 do artigo
#' Usa rugarch como backend moderno; fallback para variância móvel
#' @param resid  vetor T de resíduos
#' @return sigma2 vetor T de variâncias condicionais (normalizadas com média=1)
estimate_sigma2 <- function(resid) {
  T_obs  <- length(resid)
  sigma2 <- rep(1, T_obs)

  # Tenta GARCH(1,1) via rugarch (moderno, substitui fGarch)
  fit_ok <- FALSE
  if (requireNamespace("rugarch", quietly = TRUE)) {
    tryCatch({
      spec <- rugarch::ugarchspec(
        variance.model  = list(model = "sGARCH", garchOrder = c(1, 1)),
        mean.model      = list(armaOrder = c(0, 0), include.mean = TRUE),
        distribution.model = "norm"
      )
      fit  <- rugarch::ugarchfit(spec = spec, data = resid,
                                  solver = "hybrid", solver.control = list(trace = 0))
      sigma2  <- as.numeric(rugarch::sigma(fit))^2
      fit_ok  <- TRUE
    }, error = function(e) NULL)
  }

  # Fallback: variância móvel (janela = 12) se GARCH falhar
  if (!fit_ok) {
    for (t in 1:T_obs) {
      window   <- max(1, t - 11):t
      sigma2[t] <- var(resid[window], na.rm = TRUE)
    }
    sigma2 <- pmax(sigma2, 1e-8)
  }

  # Normaliza: mean(sigma2) = 1 (conforme artigo)
  sigma2 / mean(sigma2)
}


#' Estima variâncias dos parâmetros omega_k — passo 3 do artigo
#' omega_k = (1/(T-1)) * sum_t (Δbeta_k,t)^2
#' @param beta  T×K matrix de betas
#' @return omega vetor K (normalizado com média = 1)
estimate_omega <- function(beta) {
  T_obs  <- nrow(beta)
  K      <- ncol(beta)
  # Primeira diferença dos betas
  dbeta  <- diff(beta)  # (T-1) x K
  omega  <- colMeans(dbeta^2)
  omega  <- pmax(omega, 1e-12)
  # Normaliza: mean(omega) = 1
  omega / mean(omega)
}


# ==============================================================
# 6. Step 2 do 2SRR: solução GLS ponderada (Eq. 11 do artigo)
# ==============================================================
# θ̃ = C Ω^{1/2} Z̃' (Z̃Z̃' + λI)^{-1} ỹ
# Z̃ = Σ^{-1/2} Z Ω^{1/2},  ỹ = Σ^{-1/2} y
# O que equivale a reescalar: X̃[t,k] = X[t,k] * sqrt(omega[k]) / sigma[t]
#                              ỹ[t]    = y[t] / sigma[t]
# e rodar o ridge homogêneo em X̃, ỹ

tvp_2srr <- function(X, y, kfold = 5,
                      lambdas = exp(seq(-4, 20, length.out = 30))) {

  # --- Passo 1: ridge homogêneo ---
  step1 <- tvp_1srr(X, y, kfold = kfold, lambdas = lambdas)

  # --- Passo 2: estima sigma2 e omega ---
  sigma2 <- estimate_sigma2(step1$resid)
  omega  <- estimate_omega(step1$beta)

  # --- GLS rescaling ---
  inv_sigma <- 1 / sqrt(pmax(sigma2, 1e-8))
  sqrt_omega <- sqrt(pmax(omega, 1e-12))

  # X̃[t,k] = X[t,k] * sqrt(omega_k) / sigma_t
  X_tilde <- sweep(X, 2, sqrt_omega, "*")        # escala colunas por omega
  X_tilde <- sweep(X_tilde, 1, inv_sigma, "*")   # escala linhas por 1/sigma
  y_tilde <- y * inv_sigma

  # --- Passo 2 CV + estimação final ---
  lam_opt2 <- cv_ridge_dual(X_tilde, y_tilde, lambdas = lambdas, kfold = kfold)
  zz2      <- make_ZZt(X_tilde, y_tilde)
  alpha2   <- dual_solve(zz2$ZZt, zz2$Zty, lam_opt2)

  # Recupera betas no espaço original:
  # β_k(t) = sqrt(omega_k) * beta_tilde_k(t) (desescala)
  beta_tilde <- recover_beta(X_tilde, alpha2, zz2$C0)
  beta_orig  <- sweep(beta_tilde, 2, sqrt_omega, "*")

  resid2 <- y - rowSums(X * beta_orig)

  list(
    beta   = beta_orig,      # T×K betas time-varying no espaço original
    resid  = resid2,
    lambda = lam_opt2,
    omega  = omega,
    sigma2 = sigma2,
    # step 1 info
    beta_step1  = step1$beta,
    lambda_step1 = step1$lambda
  )
}


# ==============================================================
# 7. Previsão OOS com 2SRR
# ==============================================================
#' Gera previsão h-passos à frente com 2SRR
#' @param X_in   T×K — regressores in-sample
#' @param y_in   T   — variável alvo in-sample
#' @param x_out  K   — regressores OOS (um vetor)
#' @param kfold  número de folds para CV
#' @param lambdas vetor de lambdas candidatos
#' @return list(forecast, beta, lambda, omega, sigma2)
forecast_2srr <- function(X_in, y_in, x_out,
                           kfold   = 5,
                           lambdas = exp(seq(-4, 20, length.out = 25))) {

  # Valida inputs
  stopifnot(is.matrix(X_in), length(y_in) == nrow(X_in),
            length(x_out) == ncol(X_in))
  stopifnot(all(is.finite(X_in)), all(is.finite(y_in)),
            all(is.finite(x_out)))

  # Remove colunas com variância zero
  col_var <- apply(X_in, 2, var, na.rm = TRUE)
  good_cols <- col_var > 1e-10
  if (sum(good_cols) < 2) {
    warning("Menos de 2 colunas com variância > 0. Retornando NA.")
    return(list(forecast = NA_real_, beta = NULL, lambda = NA,
                omega = NULL, sigma2 = NULL))
  }
  X_in  <- X_in[, good_cols, drop = FALSE]
  x_out <- x_out[good_cols]

  # Padroniza X (mean=0, sd=1) para estabilidade numérica
  X_means <- colMeans(X_in)
  X_sds   <- apply(X_in, 2, sd)
  X_sds[X_sds < 1e-10] <- 1
  X_in_sc  <- sweep(sweep(X_in, 2, X_means, "-"), 2, X_sds, "/")
  x_out_sc <- (x_out - X_means) / X_sds

  # Roda 2SRR
  fit <- tvp_2srr(X_in_sc, y_in, kfold = kfold, lambdas = lambdas)

  # Previsão OOS: usa beta do último período in-sample
  fcast <- sum(x_out_sc * fit$beta[nrow(fit$beta), ])

  # Sanity check: fallback para ridge simples se previsão implausível
  y_mean <- mean(y_in)
  y_sd   <- sd(y_in)
  if (!is.finite(fcast) || abs(fcast - y_mean) > 5 * y_sd) {
    tryCatch({
      cv_r    <- glmnet::cv.glmnet(X_in_sc, y_in, alpha = 0,
                                    nfolds = kfold)
      fcast_r <- as.numeric(predict(
        glmnet::glmnet(X_in_sc, y_in, alpha = 0,
                       lambda = cv_r$lambda.min),
        newx = matrix(x_out_sc, nrow = 1)
      ))
      fcast <- fcast_r
      warning("2SRR instável — usando Ridge como fallback.")
    }, error = function(e) fcast <<- y_mean)
  }

  # Devolve betas no espaço padronizado (suficiente para análise)
  list(
    forecast = fcast,
    beta     = fit$beta,        # T×K time-varying betas
    lambda   = fit$lambda,
    omega    = fit$omega,
    sigma2   = fit$sigma2
  )
}


# ==============================================================
# 8. Wrapper padrão Medeiros: run2srr()
# run2srr() — versão alinhada ao padrão Medeiros
# Usa dataprep() diretamente.
# ==============================================================

run2srr <- function(ind, df, variable, horizon,
                    kfold   = 5,
                    lambdas = exp(seq(-4, 20, length.out = 25))) {

  # --- 1. Mesmo dataprep do Medeiros (nofact=TRUE: sem PCA, com dummy) ---
  prep <- dataprep(ind, df, variable, horizon, nofact = TRUE)
  Xin  <- prep$Xin    # matrix (T-h) x K  — idêntico ao Ridge
  yin  <- prep$yin    # vector (T-h)       — idêntico ao Ridge
  Xout <- as.numeric(prep$Xout)  # vetor K — ponto de previsão OOS

  if (length(yin) < 20 || ncol(Xin) < 2)
    return(list(forecast = NA_real_, outputs = NULL))

  # --- 2. Remove colunas com variância zero (robustez) ---
  col_var   <- apply(Xin, 2, var, na.rm = TRUE)
  good_cols <- col_var > 1e-10
  if (sum(good_cols) < 2)
    return(list(forecast = NA_real_, outputs = NULL))
  Xin  <- Xin[, good_cols, drop = FALSE]
  Xout <- Xout[good_cols]

  # --- 3. Padroniza X (para estabilidade numérica do dual ridge) ---
  X_means <- colMeans(Xin)
  X_sds   <- apply(Xin, 2, sd)
  X_sds[X_sds < 1e-10] <- 1
  Xin_sc   <- sweep(sweep(Xin, 2, X_means, "-"), 2, X_sds, "/")
  Xout_sc  <- (Xout - X_means) / X_sds

  # --- 4. Estima 2SRR (dual ridge, dois passos) ---
  fit <- tvp_2srr(Xin_sc, yin, kfold = kfold, lambdas = lambdas)

  # --- 5. Previsão OOS com betas do último período ---
  beta_T <- fit$beta[nrow(fit$beta), ]   # coefs finais (escala padronizada)
  fcast  <- sum(Xout_sc * beta_T)

  # Fallback Ridge simples se previsão for implausível
  y_sd <- sd(yin)
  y_mu <- mean(yin)
  if (!is.finite(fcast) || abs(fcast - y_mu) > 5 * y_sd) {
    tryCatch({
      cv_r  <- glmnet::cv.glmnet(Xin_sc, yin, alpha = 0, nfolds = kfold)
      fcast <- as.numeric(predict(
        glmnet::glmnet(Xin_sc, yin, alpha = 0, lambda = cv_r$lambda.min),
        newx = matrix(Xout_sc, nrow = 1)
      ))
    }, error = function(e) fcast <<- y_mu)
  }

  # --- 6. Despadroniza betas para escala original (comparável ao Ridge) ---
  # beta_orig[k] = beta_sc[k] / sd(Xin[,k])
  # (intercepto implícito absorvido na média)
  beta_orig <- sweep(fit$beta, 2, X_sds[good_cols], "/")

  # --- 7. Outputs no mesmo formato que os outros modelos do Medeiros ---
  outputs <- list(
    betas_time_varying = beta_orig,   # T×K — betas no tempo (escala original)
    lambda             = fit$lambda,
    omega              = fit$omega,   # variâncias dos parâmetros (2SRR step 2)
    sigma2             = fit$sigma2,  # variâncias dos resíduos (GARCH step 2)
    n_obs              = nrow(Xin),
    n_vars             = ncol(Xin)
  )

  list(forecast = fcast, outputs = outputs)
}

# ==============================================================
# 9. Utilitários de análise dos betas no tempo
# ==============================================================

#' Consolida betas de todas as janelas em data.frame para análise
#' @param model_list resultado de rolling_window(run2srr, ...)
#' @param df         data matrix original (com rownames = datas)
#' @param nwindows   número de janelas OOS
extract_betas_over_time <- function(model_list, df, nwindows) {
  n_windows  <- length(model_list$outputs)
  dates_all  <- as.Date(rownames(df))
  win_dates  <- tail(dates_all, n_windows)
  betas_list <- vector("list", n_windows)

  for (i in seq_len(n_windows)) {
    out <- model_list$outputs[[i]]
    if (is.null(out) || is.null(out$betas_time_varying)) next
    b   <- out$betas_time_varying
    last_b <- if (is.matrix(b)) as.numeric(b[nrow(b), ]) else as.numeric(b)
    betas_list[[i]] <- data.frame(
      window_date = win_dates[i],
      var_idx     = seq_along(last_b),
      beta        = last_b
    )
  }
  do.call(rbind, betas_list)
}


#' Plot dos top-N betas com maior variância no tempo
#' @param df_betas  data.frame de extract_betas_over_time()
#' @param var_names vetor de nomes das variáveis (opcional)
#' @param top_n     quantas variáveis plotar
#' @param save_path caminho para salvar PNG (opcional)
plot_betas_over_time <- function(df_betas, var_names = NULL,
                                  top_n = 10, variable = "CPIAUCSL",
                                  save_path = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Instale ggplot2")

  beta_var  <- tapply(df_betas$beta, df_betas$var_idx, var, na.rm = TRUE)
  top_vars  <- order(beta_var, decreasing = TRUE)[seq_len(min(top_n,
                                                               length(beta_var)))]
  df_p <- df_betas[df_betas$var_idx %in% top_vars, ]
  df_p$var_label <- if (!is.null(var_names) &&
                         max(df_p$var_idx) <= length(var_names)) {
    var_names[df_p$var_idx]
  } else paste0("X", df_p$var_idx)

  p <- ggplot2::ggplot(df_p,
         ggplot2::aes(x = window_date, y = beta,
                      colour = var_label, group = var_label)) +
    ggplot2::geom_line(linewidth = 0.7, alpha = 0.85) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        colour = "grey60", linewidth = 0.4) +
    ggplot2::scale_x_date(date_labels = "%Y", date_breaks = "2 years") +
    ggplot2::labs(
      title    = sprintf("Betas Time-Varying — 2SRR | %s", variable),
      subtitle = sprintf("Top %d variáveis por variância dos betas", top_n),
      x = NULL, y = "Coeficiente β(t)", colour = "Variável"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "right",
                   panel.grid.minor = ggplot2::element_blank())

  if (!is.null(save_path))
    ggplot2::ggsave(save_path, p, width = 11, height = 5, dpi = 150)
  p
}