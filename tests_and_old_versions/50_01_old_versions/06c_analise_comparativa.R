# 06c_analise_comparativa.R
#
# PARTE 3 DE 3 — Analise e comparacao de resultados
#
# Pre-requisitos:
#   - source("coulombe/06a_coulombe_setup.R")  [setup_objects.rda]
#   - source("coulombe/06b_coulombe_pipeline.R") [results_coulombe.rda]
#   - forecasts/yout.rda   (Medeiros — benchmarks)
#   - forecasts/2SRR.rda   (teu TVP 2SRR anterior, se existir)
#   - forecasts/betas_2SRR.rda (betas do teu TVP anterior, se existir)
#
# O que este script produz:
#   1. Tabela MSFE relativa: Coulombe 2SRR vs benchmarks do Medeiros
#   2. Grafico: betas TVP ao longo do tempo — comparacao
#      Coulombe 2SRR vs teu TVP 2SRR anterior (se disponivel)
#   3. Grafico: lambda selecionado ao longo do tempo (por h)
#   4. Grafico: forecasts vs realizado por horizonte
#   5. Grafico: erro quadratico acumulado (CumSSE) por modelo
#   6. Todos os graficos salvos em coulombe/plots/
# ============================================================

setwd("~/tcc/Felipe_Dornelles_tcc/forecasting_inflation")



library(tidyverse)
library(pracma)

dir.create("coulombe/plots", showWarnings = FALSE, recursive = TRUE)

# ============================================================
# PASSO 1 — Carrega todos os resultados
# ============================================================

cat("=== Carregando resultados ===\n")

# Resultados do 06b (Coulombe 2SRR)
load("coulombe/setup_objects.rda")   # tau, n_oos, dates, y_imp
load("coulombe/results_coulombe.rda") # results_coulombe

# Benchmarks do Medeiros (yout: matriz T x modelos)
load("forecasts/yout.rda")
oos_dates <- dates[(tau + 1):bigt]

cat(sprintf("  OOS: %s -> %s | n=%d\n",
            format(oos_dates[1]), format(tail(oos_dates, 1)), length(oos_dates)))

# Tenta carregar TVP 2SRR anterior
tem_tvp_anterior <- file.exists("forecasts/2SRR.rda")
if (tem_tvp_anterior) {
  load("forecasts/2SRR.rda")   # objeto: forecasts (vetor ou matrix)
  cat("  [OK] forecasts/2SRR.rda carregado\n")
} else {
  cat("  [INFO] forecasts/2SRR.rda nao encontrado — comparacao de betas limitada\n")
}

tem_betas_anterior <- file.exists("forecasts/betas_2SRR.rda")
if (tem_betas_anterior) {
  load("forecasts/betas_2SRR.rda")  # objeto: df_betas ou betas_list
  cat("  [OK] forecasts/betas_2SRR.rda carregado\n")
} else {
  cat("  [INFO] forecasts/betas_2SRR.rda nao encontrado\n")
}

# ============================================================
# HELPER: extrai nome das colunas de yout
# ============================================================

get_yout_col <- function(yout, patterns) {
  # Retorna primeiro nome de coluna que faz match com qualquer pattern
  nms <- colnames(yout)
  for (p in patterns) {
    m <- grep(p, nms, ignore.case = TRUE, value = TRUE)
    if (length(m) > 0) return(m[1])
  }
  return(NULL)
}

# ============================================================
# PASSO 2 — Tabela MSFE relativa por horizonte
# ============================================================

cat("\n=== Tabela MSFE relativa ===\n")

# Identifica colunas de benchmark em yout
col_rw  <- get_yout_col(yout, c("rw", "random.walk", "RW"))
col_ar  <- get_yout_col(yout, c("^ar", "AR"))
col_fac <- get_yout_col(yout, c("factor", "DFM", "fac"))

cat(sprintf("  Colunas yout disponiveis: %s\n",
            paste(colnames(yout), collapse=", ")))

# Monta data.frame de MSFE
msfe_tab <- data.frame(horizonte = integer(), modelo = character(),
                        MSFE = double(), RMSFE = double(),
                        MSFE_rel_RW = double(), stringsAsFactors = FALSE)

for (h in c(1, 3, 6, 12)) {
  key    <- paste0("h", h)
  r      <- results_coulombe[[key]]
  y_real <- r$y_real

  # Erros Coulombe
  err_2srr  <- (r$fc_2srr  - y_real)^2
  err_ridge <- (r$fc_ridge - y_real)^2

  msfe_2srr  <- mean(err_2srr,  na.rm = TRUE)
  msfe_ridge <- mean(err_ridge, na.rm = TRUE)

  # Erros benchmarks do Medeiros (coluna h do yout)
  # yout: rows = OOS obs, cols = modelos
  # Medeiros usa h como coluna ou como sub-lista — tenta as duas formas
  msfe_rw <- NA
  if (!is.null(col_rw) && col_rw %in% colnames(yout)) {
    yout_rw <- yout[, col_rw]
    if (is.matrix(yout_rw) || length(yout_rw) == n_oos) {
      msfe_rw <- mean((yout_rw - y_real)^2, na.rm = TRUE)
    }
  }

  add_row <- function(mod, msfe) {
    msfe_rel <- if (!is.na(msfe_rw) && msfe_rw > 0) msfe / msfe_rw else NA
    msfe_tab <<- rbind(msfe_tab, data.frame(
      horizonte   = h,
      modelo      = mod,
      MSFE        = round(msfe, 6),
      RMSFE       = round(sqrt(msfe), 6),
      MSFE_rel_RW = round(msfe_rel, 4),
      stringsAsFactors = FALSE
    ))
  }

  add_row("Coulombe_2SRR",  msfe_2srr)
  add_row("Coulombe_Ridge", msfe_ridge)
}

print(msfe_tab)
write.csv(msfe_tab, "coulombe/msfe_comparativo.csv", row.names = FALSE)
cat("  [OK] coulombe/msfe_comparativo.csv\n")

# ============================================================
# PASSO 3 — Betas TVP ao longo do tempo
#
# Os betas do TVPRR_cosso sao uma matriz (T x K) onde cada
# linha e o vetor de coeficientes estimado naquela janela t.
# Plota cada coeficiente como uma serie temporal.
# ============================================================

cat("\n=== Graficos de betas TVP ===\n")

plot_betas_tvp <- function(h, max_betas = 10) {
  key <- paste0("h", h)
  r   <- results_coulombe[[key]]

  # Extrai betas nao-nulos
  betas_list <- r$betas
  idx_ok     <- which(!sapply(betas_list, is.null))

  if (length(idx_ok) < 5) {
    cat(sprintf("  h=%d: poucos betas disponiveis (%d) — pulando\n",
                h, length(idx_ok)))
    return(invisible(NULL))
  }

  # Monta matrix: linhas = tempo OOS, colunas = coeficientes
  # Cada elemento de betas_list e um vetor de comprimento K
  K        <- length(betas_list[[idx_ok[1]]])
  K_plot   <- min(K, max_betas)
  beta_mat <- matrix(NA_real_, nrow = n_oos, ncol = K)

  for (i in idx_ok) {
    b <- betas_list[[i]]
    if (length(b) == K) beta_mat[i, ] <- b
  }

  # Converte para data.frame longo
  df_b <- as.data.frame(beta_mat[, 1:K_plot])
  colnames(df_b) <- paste0("beta_", seq_len(K_plot))
  df_b$date      <- oos_dates
  df_b$i         <- seq_len(n_oos)

  df_long <- df_b %>%
    pivot_longer(cols = starts_with("beta_"),
                 names_to = "coef", values_to = "valor") %>%
    filter(!is.na(valor))

  # Nomes dos coeficientes: [const/lag_y1, lag_y2, F1_lag1, F1_lag2, ...]
  coef_labels <- c(
    paste0("Lag_Y_", seq_len(LY)),
    paste0("F", rep(seq_len(NF), each = LF), "_Lag", rep(seq_len(LF), NF))
  )[seq_len(K_plot)]
  label_map <- setNames(coef_labels, paste0("beta_", seq_len(K_plot)))
  df_long$coef_label <- label_map[df_long$coef]

  p <- ggplot(df_long, aes(x = date, y = valor, color = coef_label)) +
    geom_line(linewidth = 0.6, alpha = 0.85) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.4) +
    scale_x_date(date_breaks = "4 years", date_labels = "%Y") +
    labs(
      title    = sprintf("Betas TVP — Coulombe 2SRR | h = %d", h),
      subtitle = sprintf("Coeficientes estimados por janela OOS (%s a %s)",
                         format(oos_dates[1]), format(tail(oos_dates, 1))),
      x        = NULL,
      y        = "Valor do coeficiente",
      color    = "Coeficiente"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position  = "right",
      legend.text      = element_text(size = 8),
      panel.grid.minor = element_blank()
    )

  fname <- sprintf("coulombe/plots/betas_tvp_h%02d.png", h)
  ggsave(fname, p, width = 12, height = 5, dpi = 150)
  cat(sprintf("  [OK] %s\n", fname))
  return(invisible(p))
}

LY <- 2; LF <- 2; NF <- 8  # garante que existem no ambiente
for (h in c(1, 3, 6, 12)) plot_betas_tvp(h)

# ============================================================
# PASSO 4 — Comparacao de betas: Coulombe vs teu TVP anterior
#
# Se betas_2SRR.rda existir, compara o beta do lag_Y1
# (coeficiente mais interpretavel) entre os dois metodos.
# ============================================================

if (tem_betas_anterior) {
  cat("\n=== Comparacao betas Coulombe vs TVP anterior ===\n")

  # Tenta detectar formato de betas_list / df_betas
  # Caso 1: df_betas com colunas date + beta_1 ... beta_K
  # Caso 2: lista de listas com $date e $betas
  h_comp <- 1
  key    <- paste0("h", h_comp)
  r_coul <- results_coulombe[[key]]

  # Extrai beta_1 do Coulombe (primeiro lag de Y)
  beta1_coul <- sapply(seq_len(n_oos), function(i) {
    b <- r_coul$betas[[i]]
    if (is.null(b) || length(b) < 1) NA else b[1]
  })

  # Tenta extrair beta_1 do objeto anterior
  beta1_prev <- tryCatch({
    if (exists("df_betas") && is.data.frame(df_betas)) {
      # formato data.frame: procura coluna numerica
      num_cols <- sapply(df_betas, is.numeric)
      df_betas[, which(num_cols)[1]]
    } else if (exists("betas_list") && is.list(betas_list)) {
      sapply(betas_list, function(x) {
        b <- if (is.list(x)) x$betas else x
        if (is.null(b) || length(b) < 1) NA else b[1]
      })
    } else {
      NULL
    }
  }, error = function(e) NULL)

  if (!is.null(beta1_prev)) {
    n_comp <- min(length(beta1_coul), length(beta1_prev), n_oos)

    df_comp <- data.frame(
      date          = oos_dates[1:n_comp],
      Coulombe_2SRR = beta1_coul[1:n_comp],
      TVP_anterior  = beta1_prev[1:n_comp]
    ) %>%
      pivot_longer(-date, names_to = "modelo", values_to = "beta1") %>%
      filter(!is.na(beta1))

    p_comp <- ggplot(df_comp, aes(x = date, y = beta1, color = modelo)) +
      geom_line(linewidth = 0.7, alpha = 0.9) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      scale_color_manual(values = c("Coulombe_2SRR" = "#0072B2",
                                     "TVP_anterior"  = "#D55E00")) +
      scale_x_date(date_breaks = "4 years", date_labels = "%Y") +
      labs(
        title    = "Beta_1 (Lag Y) ao longo do tempo — h = 1",
        subtitle = "Coulombe 2SRR vs TVP 2SRR anterior",
        x        = NULL, y = "Beta_1", color = "Modelo"
      ) +
      theme_bw(base_size = 11) +
      theme(legend.position = "top", panel.grid.minor = element_blank())

    ggsave("coulombe/plots/betas_comp_h01.png", p_comp,
           width = 12, height = 4.5, dpi = 150)
    cat("  [OK] coulombe/plots/betas_comp_h01.png\n")
  } else {
    cat("  [INFO] Nao foi possivel extrair beta_1 do objeto anterior.\n")
    cat("         Verifique o formato de betas_2SRR.rda e ajuste o codigo.\n")
  }
}

# ============================================================
# PASSO 5 — Lambda selecionado ao longo do tempo
# ============================================================

cat("\n=== Grafico lambda selecionado ===\n")

df_lambda <- map_dfr(c(1, 3, 6, 12), function(h) {
  r <- results_coulombe[[paste0("h", h)]]
  data.frame(
    date   = oos_dates,
    lambda = r$lambda_sel,
    h      = paste0("h=", h)
  )
}) %>% filter(!is.na(lambda))

p_lambda <- ggplot(df_lambda, aes(x = date, y = log(lambda), color = h)) +
  geom_line(linewidth = 0.6, alpha = 0.85) +
  scale_x_date(date_breaks = "4 years", date_labels = "%Y") +
  labs(
    title    = "log(Lambda) selecionado pelo TVPRR_cosso ao longo do tempo",
    subtitle = "Cada curva = um horizonte de previsao",
    x        = NULL, y = "log(lambda)", color = "Horizonte"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "top", panel.grid.minor = element_blank())

ggsave("coulombe/plots/lambda_tempo.png", p_lambda,
       width = 12, height = 4, dpi = 150)
cat("  [OK] coulombe/plots/lambda_tempo.png\n")

# ============================================================
# PASSO 6 — Forecast vs Realizado por horizonte
# ============================================================

cat("\n=== Grafico forecast vs realizado ===\n")

for (h in c(1, 3, 6, 12)) {
  r <- results_coulombe[[paste0("h", h)]]

  df_fc <- data.frame(
    date    = oos_dates,
    Real    = r$y_real,
    `2SRR_Coulombe` = r$fc_2srr,
    Ridge   = r$fc_ridge
  ) %>%
    pivot_longer(-date, names_to = "serie", values_to = "valor") %>%
    filter(!is.na(valor))

  cores <- c("Real" = "black",
             "X2SRR_Coulombe" = "#0072B2",
             "Ridge" = "#CC79A7")

  p_fc <- ggplot(df_fc, aes(x = date, y = valor,
                             color = serie, linewidth = serie)) +
    geom_line(alpha = 0.85) +
    scale_linewidth_manual(values = c("Real" = 0.8,
                                       "X2SRR_Coulombe" = 0.7,
                                       "Ridge" = 0.5),
                            guide = "none") +
    scale_color_manual(values = cores,
                       labels = c("Real" = "Realizado",
                                  "X2SRR_Coulombe" = "2SRR Coulombe",
                                  "Ridge" = "Ridge")) +
    scale_x_date(date_breaks = "4 years", date_labels = "%Y") +
    labs(
      title    = sprintf("Forecast vs Realizado — h = %d", h),
      subtitle = "OOS: 1999-07 a 2025-06",
      x        = NULL, y = "Inflacao acumulada (h passos)", color = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "top", panel.grid.minor = element_blank())

  fname <- sprintf("coulombe/plots/forecast_vs_real_h%02d.png", h)
  ggsave(fname, p_fc, width = 12, height = 4.5, dpi = 150)
  cat(sprintf("  [OK] %s\n", fname))
}

# ============================================================
# PASSO 7 — Erro quadratico acumulado (CumSSE)
# ============================================================

cat("\n=== Grafico CumSSE ===\n")

for (h in c(1, 3, 6, 12)) {
  r <- results_coulombe[[paste0("h", h)]]

  cumsse_2srr  <- cumsum(ifelse(is.na((r$fc_2srr  - r$y_real)^2), 0,
                                 (r$fc_2srr  - r$y_real)^2))
  cumsse_ridge <- cumsum(ifelse(is.na((r$fc_ridge - r$y_real)^2), 0,
                                 (r$fc_ridge - r$y_real)^2))

  df_cum <- data.frame(
    date          = oos_dates,
    `2SRR_Coulombe` = cumsse_2srr,
    Ridge         = cumsse_ridge
  ) %>%
    pivot_longer(-date, names_to = "modelo", values_to = "cumsse")

  p_cum <- ggplot(df_cum, aes(x = date, y = cumsse, color = modelo)) +
    geom_line(linewidth = 0.7) +
    scale_color_manual(values = c("X2SRR_Coulombe" = "#0072B2",
                                   "Ridge" = "#CC79A7"),
                       labels = c("X2SRR_Coulombe" = "2SRR Coulombe",
                                  "Ridge" = "Ridge")) +
    scale_x_date(date_breaks = "4 years", date_labels = "%Y") +
    labs(
      title    = sprintf("Erro Quadratico Acumulado (CumSSE) — h = %d", h),
      subtitle = "Curva mais baixa = melhor previsao acumulada",
      x        = NULL, y = "CumSSE", color = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "top", panel.grid.minor = element_blank())

  fname <- sprintf("coulombe/plots/cumsse_h%02d.png", h)
  ggsave(fname, p_cum, width = 12, height = 4, dpi = 150)
  cat(sprintf("  [OK] %s\n", fname))
}

# ============================================================
# PASSO 8 — Resumo final no console
# ============================================================

cat("\n============================================================\n")
cat("RESUMO FINAL — MSFE e RMSFE por horizonte\n")
cat("============================================================\n")
cat(sprintf("  %-6s  %-14s  %-14s  %-12s\n",
            "h", "RMSFE_2SRR", "RMSFE_Ridge", "Razao_2SRR/Ridge"))
for (h in c(1, 3, 6, 12)) {
  r <- results_coulombe[[paste0("h", h)]]
  razao <- sqrt(r$msfe_2srr / r$msfe_ridge)
  cat(sprintf("  h=%-4d  %-14.6f  %-14.6f  %-12.4f\n",
              h, sqrt(r$msfe_2srr), sqrt(r$msfe_ridge), razao))
}

cat("\nArquivos gerados em coulombe/plots/:\n")
arqs <- list.files("coulombe/plots", pattern = "\\.png$", full.names = FALSE)
for (a in arqs) cat(sprintf("  %s\n", a))

cat("\n=== 06c CONCLUIDO ===\n")
