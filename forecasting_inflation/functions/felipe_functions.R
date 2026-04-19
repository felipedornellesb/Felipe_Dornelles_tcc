# ============================================================
# felipe_functions.R
#
# Funções no padrão Medeiros (ForecastingInflation) para rodar
# o modelo 2SRR de Coulombe (2019).
#
# Interface obrigatória (idêntica a functions.R do Medeiros):
#   fn(ind, df, variable, horizon, ...)
#   retorna list(forecast, outputs)
#
# DEPENDÊNCIAS — fazer source ANTES de usar este arquivo:
#   source("functions/functions.R")          # dataprep() do Medeiros
#   source("coulombe/TVPRRcosso_v181120.R")
#   source("coulombe/zfun_v190304.R")
#   source("coulombe/fastZrot_v181125.R")
#   source("coulombe/CVGSBHK_v181127.R")
# ============================================================


# ============================================================
# run2srr()
#
# Roda o modelo 2SRR (Two-Step Ridge Regression, Coulombe 2019)
# usando a infraestrutura do Medeiros.
#
# Parâmetros:
#   ind         : índice de linhas (janela de treino), igual ao Medeiros
#   df          : matriz de dados com rownames = datas, igual ao Medeiros
#   variable    : variável alvo (ex: "CPIAUCSL")
#   horizon     : horizonte de previsão (1 a 12)
#   lambda_vec  : grid de lambdas para o COSSO/2SRR
#   kfold       : número de folds para cross-validation
#   alpha_2srr  : parâmetro alpha do COSSO (padrão Coulombe = 0.01)
#   silent      : 1 = sem output no console durante estimação
#
# Retorna:
#   list(
#     forecast  : previsão escalar para o horizonte h
#     outputs   : list com betas_time_varying (matrix T x K),
#                 lambda selecionado por CV, e metadados
#   )
# ============================================================

run2srr <- function(ind, df, variable, horizon,
                    lambda_vec  = exp(seq(-6, 20, length.out = 15)),
                    kfold       = 5,
                    alpha_2srr  = 0.01,
                    silent      = 1) {

  # --- 1. Prepara dados usando dataprep() do Medeiros (nofact = TRUE) ---
  # nofact = TRUE: sem PCA, usa todas as variáveis diretamente
  # O 2SRR lida com alta dimensão nativamente via regularização
  prep  <- dataprep(ind, df, variable, horizon, nofact = TRUE)
  Xin   <- prep$Xin
  yin   <- prep$yin
  Xout  <- as.numeric(prep$Xout)

  # --- 2. Ridge via glmnet como fallback e para obter lambda2 inicial ---
  # lambda2 é o parâmetro de Ridge do segundo estágio do 2SRR
  pred_ridge <- NA_real_
  lv2        <- 1.0

  tryCatch({
    cv_r       <- glmnet::cv.glmnet(x = Xin, y = yin,
                                     alpha = 0, nfolds = kfold,
                                     family = "gaussian")
    mdl_r      <- glmnet::glmnet(x = Xin, y = yin,
                                  alpha = 0, lambda = cv_r$lambda.min,
                                  family = "gaussian")
    pred_ridge <- as.numeric(
      predict(mdl_r, newx = matrix(Xout, nrow = 1))
    )
    lv2        <- cv_r$lambda.min
  }, error = function(e) {
    message(sprintf("[run2srr | Ridge fallback] Erro: %s", e$message))
  })

  # --- 3. Estimação 2SRR via TVPRR_cosso (type = 2) ---
  forecast   <- NA_real_
  betas_tvp  <- NULL
  lambda_sel <- NA_real_

  tryCatch({
    fit <- TVPRR_cosso(
      y         = yin,
      X         = Xin,
      lambdavec = lambda_vec,
      sweigths  = 1,
      type      = 2,          # type = 2 → 2SRR (Coulombe)
      alpha     = alpha_2srr,
      silent    = silent,
      kfold     = kfold,
      lambda2   = lv2,
      tol       = 1e-6,
      maxit     = 10,
      oosX      = Xout
    )

    forecast   <- as.numeric(fit$fcast)
    lambda_sel <- fit$lambda         # lambda ótimo selecionado por CV
    betas_tvp  <- fit$beta           # matrix T x K — betas time-varying
                                     # cada linha = um período no tempo
                                     # cada coluna = uma variável

    # Outlier filter: se previsão absurda, usa Ridge como fallback
    y_mean <- mean(yin, na.rm = TRUE)
    y_sd   <- sd(yin,   na.rm = TRUE)
    if (!is.finite(forecast) ||
        abs(forecast - y_mean) > 6 * y_sd) {
      if (is.finite(pred_ridge)) {
        message("[run2srr] Previsão fora dos limites — usando Ridge como fallback.")
        forecast <- pred_ridge
      }
    }

  }, error = function(e) {
    message(sprintf("[run2srr | 2SRR] Erro: %s", e$message))
    if (is.finite(pred_ridge)) forecast <<- pred_ridge
  })

  # --- 4. Monta outputs ---
  outputs <- list(
    betas_time_varying = betas_tvp,   # PRINCIPAL: betas no tempo
    lambda             = lambda_sel,  # lambda CV ótimo
    ridge_fallback     = pred_ridge,  # benchmark interno
    n_obs              = length(yin),
    n_vars             = ncol(Xin)
  )

  return(list(forecast = forecast, outputs = outputs))
}


# ============================================================
# extract_betas_over_time()
#
# Consolida os betas time-varying de TODAS as janelas rolling
# em um único data.frame pronto para análise e visualização.
#
# O campo outputs[[i]]$betas_time_varying de cada janela é uma
# matrix T x K. Extraímos a ÚLTIMA linha (= betas usados para
# gerar a previsão daquela janela) e empilhamos no tempo.
#
# Uso:
#   model_list  <- rolling_window(run2srr, data, nwindows+h-1, h, "CPIAUCSL")
#   df_betas    <- extract_betas_over_time(model_list, data, nwindows)
#
# Retorna data.frame com colunas:
#   window_date  : data da previsão (última observação da janela)
#   var_idx      : índice da variável (coluna de Xin)
#   beta         : coeficiente time-varying naquela janela
# ============================================================

extract_betas_over_time <- function(model_list, df, nwindows) {

  n_windows  <- length(model_list$outputs)
  dates_all  <- as.Date(rownames(df))

  # data associada a cada previsão (fim da janela de treino)
  window_dates <- dates_all[(nrow(df) - n_windows + 1):nrow(df)]

  betas_list <- vector("list", n_windows)

  for (i in seq_len(n_windows)) {
    out <- model_list$outputs[[i]]
    if (is.null(out) || is.null(out$betas_time_varying)) next

    b <- out$betas_time_varying

    if (is.matrix(b) && nrow(b) > 0) {
      last_betas <- as.numeric(b[nrow(b), ])   # última linha = betas "atuais"
    } else if (is.numeric(b)) {
      last_betas <- b
    } else {
      next
    }

    betas_list[[i]] <- data.frame(
      window_date = window_dates[i],
      var_idx     = seq_along(last_betas),
      beta        = last_betas,
      stringsAsFactors = FALSE
    )
  }

  df_betas <- do.call(rbind, betas_list)
  return(df_betas)
}


# ============================================================
# plot_betas_over_time()
#
# Gera gráfico ggplot dos betas time-varying mais relevantes.
# Seleciona automaticamente os top_n por maior variância no tempo.
#
# Parâmetros:
#   df_betas   : output de extract_betas_over_time()
#   var_names  : vetor de nomes das variáveis (colnames de Xin)
#                Se NULL, usa "X1", "X2", ...
#   top_n      : quantas variáveis plotar (default = 10)
#   variable   : nome da variável alvo (para o título)
#   save_path  : caminho para salvar o PNG (NULL = não salva)
#
# Uso:
#   plot_betas_over_time(df_betas, var_names = colnames(data),
#                        top_n = 10, save_path = "results/betas.png")
# ============================================================

plot_betas_over_time <- function(df_betas,
                                  var_names  = NULL,
                                  top_n      = 10,
                                  variable   = "CPIAUCSL",
                                  save_path  = NULL) {

  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Instale ggplot2: install.packages('ggplot2')")

  # Seleciona top_n variáveis com maior variância dos betas no tempo
  beta_var <- tapply(df_betas$beta, df_betas$var_idx, var, na.rm = TRUE)
  top_vars <- order(beta_var, decreasing = TRUE)[
    seq_len(min(top_n, length(beta_var)))
  ]

  df_plot <- df_betas[df_betas$var_idx %in% top_vars, ]

  # Labels das variáveis
  if (!is.null(var_names)) {
    df_plot$var_label <- ifelse(
      df_plot$var_idx <= length(var_names),
      var_names[df_plot$var_idx],
      paste0("X", df_plot$var_idx)
    )
  } else {
    df_plot$var_label <- paste0("X", df_plot$var_idx)
  }

  p <- ggplot2::ggplot(
    df_plot,
    ggplot2::aes(x     = window_date,
                 y     = beta,
                 colour = var_label,
                 group  = var_label)
  ) +
    ggplot2::geom_line(linewidth = 0.7, alpha = 0.85) +
    ggplot2::geom_hline(yintercept = 0,
                        linetype = "dashed",
                        colour   = "grey50",
                        linewidth = 0.4) +
    ggplot2::scale_x_date(date_labels = "%Y",
                          date_breaks = "2 years") +
    ggplot2::labs(
      title    = sprintf("Betas Time-Varying — 2SRR | %s", variable),
      subtitle = sprintf("Top %d variáveis por variância dos coeficientes", top_n),
      x        = NULL,
      y        = expression(beta(t)),
      colour   = "Variável"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      legend.position  = "right",
      panel.grid.minor = ggplot2::element_blank(),
      plot.title       = ggplot2::element_text(face = "bold")
    )

  if (!is.null(save_path)) {
    ggplot2::ggsave(save_path, p, width = 12, height = 5, dpi = 150)
    message(sprintf("Gráfico salvo em: %s", save_path))
  }

  return(p)
}