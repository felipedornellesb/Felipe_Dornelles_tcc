# 11_betas_descriptive_analysis.R
#
# Analise descritiva dos betas:
#   PARTE A: Primeira janela in-sample (t=tau)
#            - Trajetoria TVP dos betas dentro da janela
#            - Comparacao com Ridge constante
#   PARTE B: Evolucao dos betas finais ao longo do OOS
#            - Ultimo beta do 2SRR em cada janela vs Ridge
#   PARTE C: Analise de parcimonia
#            - Quantos betas sao "efetivamente zero"?
#            - 2SRR e mais parcimonioso que Ridge?
# ============================================================

rm(list = ls())
setwd("~/TCC/tcc/forecasting_inflation")

library(ggplot2)
library(reshape2)
library(gridExtra)

dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

cat("  ANALISE DESCRITIVA DOS BETAS\n")

# 0. CARREGA DADOS

load("forecasts/coulombe_betas_2SRR.rda")   # betas_2srr
load("forecasts/coulombe_betas_ridge.rda")   # betas_ridge
load("data/data.rda")

fred_raw <- as.data.frame(data)
date_col <- grep("^date$", colnames(fred_raw), ignore.case = TRUE)[1]
dates <- fred_raw[, date_col]

hor <- c(1, 3, 6, 12)

# Nomes dos regressores (2 lags y + 2 lags x 8 fatores = 18 regressores)
reg_names <- c("intercept",
               paste0("y_lag", 1:2),
               paste0("F", rep(1:8, 2), "_lag", rep(1:2, each = 8)))

cat("Objetos carregados.\n")
cat(sprintf("betas_2srr:  %d horizontes, %d janelas cada\n",
            length(betas_2srr), length(betas_2srr[[1]])))
cat(sprintf("betas_ridge: %d horizontes, %d janelas cada\n\n",
            length(betas_ridge), length(betas_ridge[[1]])))

# PARTE A: PRIMEIRA JANELA IN-SAMPLE
#
# "Pegar so a primeira janela. Rodar o TVP e ver se o beta
#  da muito parecido com o constante (Ridge). Se sim, nao
#  precisa do TVP."

cat("  PARTE A: Primeira Janela In-Sample (t = tau)\n")

hi <- 1  # h=1 para analise principal
h  <- hor[hi]

# Pega a primeira janela
first_2srr  <- betas_2srr[[hi]][[1]]
first_ridge <- betas_ridge[[hi]][[1]]

if (!is.null(first_2srr) && !is.null(first_ridge)) {
  
  cat(sprintf("  Data da primeira janela: %s (t=%d)\n",
              as.character(first_ridge$date), first_ridge$t))
  
  # A1: Extrair trajetoria TVP da primeira janela
  bm <- first_2srr$betas
  cat(sprintf("  Dimensao de betas_2srr: %s\n",
              paste(dim(bm), collapse = " x ")))
  
  if (is.array(bm) && length(dim(bm)) == 3) {
    # bm e [1 x K+1 x T] -> extrair serie temporal
    K <- dim(bm)[2]  # numero de betas (incluindo intercepto)
    T_is <- dim(bm)[3]  # tamanho da janela in-sample
    
    cat(sprintf("  K=%d betas | T=%d obs na janela\n\n", K, T_is))
    
    # Cria data.frame com trajetoria de cada beta
    beta_traj <- matrix(NA, nrow = T_is, ncol = K)
    for (k in 1:K) {
      beta_traj[, k] <- bm[1, k, ]
    }
    
    # Nomes
    col_names <- if (K <= length(reg_names)) reg_names[1:K] else paste0("beta", 0:(K-1))
    colnames(beta_traj) <- col_names
    beta_traj_df <- as.data.frame(beta_traj)
    beta_traj_df$t <- 1:T_is
    
    # Ridge constante (primeira janela)
    ridge_beta0 <- first_ridge$beta0
    ridge_betas <- first_ridge$betas
    ridge_all <- c(ridge_beta0, ridge_betas)
    
    cat("  A1. Comparacao: Ultimo beta TVP vs Ridge constante\n")
    cat(sprintf("  %-15s %10s %10s %10s %10s\n",
                "Regressor", "Ridge", "2SRR_t1", "2SRR_tT", "Diff_final"))
    cat(paste(rep("-", 60), collapse = ""), "\n")
    
    tvp_first <- beta_traj[1, ]      # primeiro instante
    tvp_last  <- beta_traj[T_is, ]   # ultimo instante
    
    parcimonia_df <- data.frame(
      regressor = col_names,
      ridge = ridge_all[1:K],
      tvp_first = tvp_first,
      tvp_last = tvp_last,
      diff_abs = abs(tvp_last - ridge_all[1:K]),
      stringsAsFactors = FALSE
    )
    
    for (i in 1:min(K, 19)) {
      cat(sprintf("  %-15s %10.4f %10.4f %10.4f %10.4f\n",
                  parcimonia_df$regressor[i],
                  parcimonia_df$ridge[i],
                  parcimonia_df$tvp_first[i],
                  parcimonia_df$tvp_last[i],
                  parcimonia_df$diff_abs[i]))
    }
    
    # Correlacao entre Ridge e ultimo TVP
    cor_val <- cor(parcimonia_df$ridge, parcimonia_df$tvp_last)
    cat(sprintf("\n  Correlacao entre Ridge e ultimo beta TVP: %.4f\n", cor_val))
    
    # Distancia euclidiana relativa
    dist_rel <- sqrt(sum((parcimonia_df$ridge - parcimonia_df$tvp_last)^2)) /
      sqrt(sum(parcimonia_df$ridge^2))
    cat(sprintf("  Distancia euclidiana relativa: %.4f\n", dist_rel))
    
    if (cor_val > 0.95) {
      cat("  >>> ALTA correlacao: betas TVP similares ao Ridge -> TVP agrega pouco\n")
    } else if (cor_val > 0.7) {
      cat("  >>> Correlacao MODERADA: TVP diverge parcialmente do Ridge\n")
    } else {
      cat("  >>> BAIXA correlacao: TVP diverge significativamente do Ridge\n")
    }
    
    # A2: Grafico da trajetoria TVP na primeira janela
    
    # Seleciona top 6 betas mais variaveis
    var_betas <- apply(beta_traj, 2, var)
    top_idx <- order(var_betas, decreasing = TRUE)[1:min(6, K)]
    top_names <- col_names[top_idx]
    
    # Melt para ggplot
    plot_df <- beta_traj_df[, c("t", top_names)]
    plot_long <- melt(plot_df, id.vars = "t",
                      variable.name = "beta", value.name = "valor")
    
    # Adiciona linhas horizontais do Ridge
    ridge_lines <- data.frame(
      beta = top_names,
      ridge_val = ridge_all[top_idx]
    )
    
    p_traj <- ggplot(plot_long, aes(x = t, y = valor)) +
      geom_line(color = "steelblue", linewidth = 0.5) +
      geom_hline(data = ridge_lines,
                 aes(yintercept = ridge_val),
                 linetype = "dashed", color = "red", linewidth = 0.6) +
      facet_wrap(~ beta, scales = "free_y", ncol = 2) +
      labs(title = sprintf("Trajetoria TVP (2SRR) na Primeira Janela (h=%d, T=%d)", h, T_is),
           subtitle = "Linha vermelha = Ridge constante. Se TVP coincide com Ridge, TVP nao agrega.",
           x = "Observacao dentro da janela",
           y = "Valor do beta") +
      theme_minimal() +
      theme(strip.text = element_text(face = "bold"))
    
    ggsave("results/figures/betas_primeira_janela_h01.pdf",
           p_traj, width = 12, height = 10)
    cat("\n  Grafico salvo: results/figures/betas_primeira_janela_h01.pdf\n")
    
    # A3: Barplot comparativo Ridge vs TVP (primeira janela)
    
    comp_df <- data.frame(
      regressor = factor(col_names, levels = col_names),
      Ridge = ridge_all[1:K],
      TVP_2SRR = tvp_last
    )
    comp_long <- melt(comp_df, id.vars = "regressor",
                      variable.name = "modelo", value.name = "beta")
    
    p_bar <- ggplot(comp_long, aes(x = regressor, y = beta, fill = modelo)) +
      geom_bar(stat = "identity", position = "dodge", width = 0.7) +
      labs(title = sprintf("Betas: Ridge vs 2SRR (ultima obs) — Primeira Janela h=%d", h),
           x = "", y = "Valor do beta", fill = "") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
            legend.position = "top")
    
    ggsave("results/figures/betas_barplot_primeira_janela_h01.pdf",
           p_bar, width = 14, height = 6)
    cat("  Grafico salvo: results/figures/betas_barplot_primeira_janela_h01.pdf\n")
    
  }
} else {
  cat("  [AVISO] Primeira janela nula. Verifique betas_2srr[[1]][[1]]\n")
}

# PARTE B: EVOLUCAO OOS — Ultimo beta TVP vs Ridge
#
# "Comparar esses graficos com o ridge fora da amostra.
#  Porque ele tem varios betas, comparando o ultimo beta
#  de cada [janela] e contra o ridge."

cat("  PARTE B: Evolucao OOS dos betas (h=1)\n")

# Extrai ultimo beta do 2SRR de cada janela
valid_2srr  <- Filter(Negate(is.null), betas_2srr[[hi]])
valid_ridge <- Filter(Negate(is.null), betas_ridge[[hi]])

cat(sprintf("  Janelas validas: 2SRR=%d | Ridge=%d\n",
            length(valid_2srr), length(valid_ridge)))

if (length(valid_2srr) > 5 && length(valid_ridge) > 5) {
  
  # Extrai serie temporal do ultimo beta TVP
  tvp_series <- do.call(rbind, lapply(valid_2srr, function(b) {
    bm <- b$betas
    if (is.array(bm) && length(dim(bm)) == 3) bvec <- bm[1, , dim(bm)[3]]
    else if (is.matrix(bm)) bvec <- bm[nrow(bm), ]
    else bvec <- as.numeric(bm)
    c(t = b$t, bvec)
  }))
  
  # Extrai serie temporal do Ridge
  ridge_series <- do.call(rbind, lapply(valid_ridge, function(b) {
    c(t = b$t, b$beta0, b$betas)
  }))
  
  tvp_df <- as.data.frame(tvp_series)
  ridge_df <- as.data.frame(ridge_series)
  
  n_betas <- min(ncol(tvp_df) - 1, ncol(ridge_df) - 1)
  col_nm <- if (n_betas <= length(reg_names)) reg_names[1:n_betas] else paste0("b", 0:(n_betas-1))
  
  colnames(tvp_df)   <- c("t", col_nm)
  colnames(ridge_df) <- c("t", col_nm)
  
  for (j in seq_len(ncol(tvp_df)))
    tvp_df[, j] <- as.numeric(as.character(tvp_df[, j]))
  for (j in seq_len(ncol(ridge_df)))
    ridge_df[, j] <- as.numeric(as.character(ridge_df[, j]))
  
  # Adiciona datas
  tvp_df$date   <- dates[tvp_df$t]
  ridge_df$date <- dates[ridge_df$t]
  
  # B1: Grafico para top 6 betas mais variaveis
  var_tvp <- sapply(col_nm, function(c) var(tvp_df[[c]], na.rm = TRUE))
  top6 <- names(sort(var_tvp, decreasing = TRUE))[1:min(6, length(var_tvp))]
  
  plots_list <- list()
  for (bn in top6) {
    df_plot <- data.frame(
      date = as.Date(c(tvp_df$date, ridge_df$date)),
      valor = c(tvp_df[[bn]], ridge_df[[bn]]),
      modelo = c(rep("2SRR (ultimo beta)", nrow(tvp_df)),
                 rep("Ridge (constante)", nrow(ridge_df)))
    )
    
    p <- ggplot(df_plot, aes(x = date, y = valor, color = modelo)) +
      geom_line(linewidth = 0.5) +
      labs(title = bn, x = "", y = "", color = "") +
      theme_minimal() +
      theme(legend.position = "bottom",
            plot.title = element_text(face = "bold", size = 10))
    
    plots_list[[bn]] <- p
  }
  
  p_combined <- do.call(gridExtra::grid.arrange,
                        c(plots_list, ncol = 2,
                          top = "Evolucao OOS: Ultimo beta 2SRR vs Ridge (h=1)"))
  
  ggsave("results/figures/betas_oos_evolution_h01.pdf",
         p_combined, width = 14, height = 12)
  cat("  Grafico salvo: results/figures/betas_oos_evolution_h01.pdf\n")
  
  # B2: Salva tabela comparativa para artigo
  write.csv(tvp_df, "results/betas_tvp_last_oos_h01.csv", row.names = FALSE)
  write.csv(ridge_df, "results/betas_ridge_oos_h01.csv", row.names = FALSE)
  cat("  CSVs salvos: betas_tvp_last_oos_h01.csv e betas_ridge_oos_h01.csv\n")
}

# PARTE C: ANALISE DE PARCIMONIA
#
# "Ver se o 2SRR, ou outro TVP, foi mais ou menos parcimonioso."
#
# Parcimonia = quantos betas sao efetivamente zero (ou muito
# pequenos). Um modelo mais parcimonioso "desliga" variaveis
# irrelevantes, concentrando a informacao em poucos preditores.

cat("  PARTE C: Analise de Parcimonia\n")

if (exists("tvp_df") && exists("ridge_df")) {
  
  # Para cada janela OOS, conta quantos betas sao "efetivamente zero"
  # Criterio: |beta_k| < threshold * max(|beta|)
  threshold <- 0.05  # beta < 5% do maior em magnitude -> "quase zero"
  
  # Media dos betas ao longo do OOS
  tvp_mean  <- colMeans(tvp_df[, col_nm], na.rm = TRUE)
  ridge_mean <- colMeans(ridge_df[, col_nm], na.rm = TRUE)
  
  tvp_absmax  <- max(abs(tvp_mean))
  ridge_absmax <- max(abs(ridge_mean))
  
  tvp_near_zero  <- sum(abs(tvp_mean) < threshold * tvp_absmax)
  ridge_near_zero <- sum(abs(ridge_mean) < threshold * ridge_absmax)
  
  cat("  C1. Parcimonia media ao longo do OOS:\n")
  cat(sprintf("    2SRR:  %d/%d betas 'quase zero' (< 5%% do max)\n",
              tvp_near_zero, length(col_nm)))
  cat(sprintf("    Ridge: %d/%d betas 'quase zero' (< 5%% do max)\n",
              ridge_near_zero, length(col_nm)))
  
  if (tvp_near_zero > ridge_near_zero) {
    cat("    >>> 2SRR e MAIS parcimonioso (desliga mais variaveis)\n")
  } else if (tvp_near_zero < ridge_near_zero) {
    cat("    >>> Ridge e mais parcimonioso\n")
  } else {
    cat("    >>> Parcimonia similar\n")
  }
  
  # C2: Ratio de parcimonia por janela
  cat("\n  C2. Parcimonia por janela ao longo do tempo:\n")
  
  parc_tvp  <- numeric(nrow(tvp_df))
  parc_ridge <- numeric(nrow(ridge_df))
  
  for (i in 1:nrow(tvp_df)) {
    betas_i <- as.numeric(tvp_df[i, col_nm])
    absmax_i <- max(abs(betas_i))
    if (absmax_i > 0) parc_tvp[i] <- sum(abs(betas_i) < threshold * absmax_i)
    else parc_tvp[i] <- length(col_nm)
  }
  
  for (i in 1:nrow(ridge_df)) {
    betas_i <- as.numeric(ridge_df[i, col_nm])
    absmax_i <- max(abs(betas_i))
    if (absmax_i > 0) parc_ridge[i] <- sum(abs(betas_i) < threshold * absmax_i)
    else parc_ridge[i] <- length(col_nm)
  }
  
  cat(sprintf("    2SRR:  media %.1f betas quase-zero por janela (de %d)\n",
              mean(parc_tvp), length(col_nm)))
  cat(sprintf("    Ridge: media %.1f betas quase-zero por janela (de %d)\n",
              mean(parc_ridge), length(col_nm)))
  
  # C3: Concentracao de informacao (indice de Herfindahl)
  cat("\n  C3. Concentracao de informacao (Herfindahl):\n")
  cat("    HHI = sum(share_k^2) onde share_k = |beta_k| / sum(|beta|)\n")
  cat("    HHI alto = informacao concentrada em poucos betas (mais parcimonioso)\n")
  cat("    HHI = 1/K = totalmente disperso\n\n")
  
  hhi_tvp <- numeric(nrow(tvp_df))
  hhi_ridge <- numeric(nrow(ridge_df))
  
  for (i in 1:nrow(tvp_df)) {
    betas_i <- abs(as.numeric(tvp_df[i, col_nm]))
    total <- sum(betas_i)
    if (total > 0) {
      shares <- betas_i / total
      hhi_tvp[i] <- sum(shares^2)
    }
  }
  
  for (i in 1:nrow(ridge_df)) {
    betas_i <- abs(as.numeric(ridge_df[i, col_nm]))
    total <- sum(betas_i)
    if (total > 0) {
      shares <- betas_i / total
      hhi_ridge[i] <- sum(shares^2)
    }
  }
  
  cat(sprintf("    2SRR  HHI medio: %.4f (1/K = %.4f)\n",
              mean(hhi_tvp), 1/length(col_nm)))
  cat(sprintf("    Ridge HHI medio: %.4f\n", mean(hhi_ridge)))
  
  if (mean(hhi_tvp) > mean(hhi_ridge)) {
    cat("    >>> 2SRR concentra mais informacao (mais parcimonioso)\n")
  } else {
    cat("    >>> Ridge concentra mais (ou dispersao similar)\n")
  }
  
  # C4: Grafico de parcimonia ao longo do tempo
  n_min <- min(nrow(tvp_df), nrow(ridge_df))
  parc_plot <- data.frame(
    date = as.Date(tvp_df$date[1:n_min]),
    HHI_2SRR = hhi_tvp[1:n_min],
    HHI_Ridge = hhi_ridge[1:n_min]
  )
  parc_long <- melt(parc_plot, id.vars = "date",
                    variable.name = "modelo", value.name = "HHI")
  
  p_hhi <- ggplot(parc_long, aes(x = date, y = HHI, color = modelo)) +
    geom_line(linewidth = 0.5) +
    geom_hline(yintercept = 1/length(col_nm), linetype = "dashed",
               color = "gray50", linewidth = 0.3) +
    annotate("text", x = min(parc_plot$date), y = 1/length(col_nm) + 0.005,
             label = "1/K (dispersao maxima)", hjust = 0, size = 3, color = "gray40") +
    labs(title = "Indice de Concentracao (HHI) ao longo do OOS — h=1",
         subtitle = "HHI maior = informacao concentrada em poucos betas (mais parcimonioso)",
         x = "", y = "HHI", color = "") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  ggsave("results/figures/parcimonia_hhi_h01.pdf",
         p_hhi, width = 12, height = 5)
  cat("\n  Grafico salvo: results/figures/parcimonia_hhi_h01.pdf\n")
  
  # C5: Ratio de magnitude |beta_2SRR| / |beta_Ridge|
  cat("\n  C4. Shrinkage relativo (media OOS):\n")
  cat(sprintf("    %-15s %10s %10s %10s\n",
              "Regressor", "|Ridge|", "|2SRR|", "Ratio"))
  cat(paste(rep("-", 50), collapse = ""), "\n")
  
  shrink_df <- data.frame(
    regressor = col_nm,
    abs_ridge = abs(ridge_mean),
    abs_tvp   = abs(tvp_mean),
    ratio = abs(tvp_mean) / pmax(abs(ridge_mean), 1e-10),
    stringsAsFactors = FALSE
  )
  shrink_df <- shrink_df[order(-shrink_df$abs_ridge), ]
  
  for (i in 1:nrow(shrink_df)) {
    cat(sprintf("    %-15s %10.4f %10.4f %10.3f\n",
                shrink_df$regressor[i], shrink_df$abs_ridge[i],
                shrink_df$abs_tvp[i], shrink_df$ratio[i]))
  }
  
  avg_ratio <- mean(shrink_df$ratio, na.rm = TRUE)
  cat(sprintf("\n    Ratio medio: %.3f\n", avg_ratio))
  if (avg_ratio < 0.9) {
    cat("    >>> 2SRR aplica MAIS shrinkage que Ridge (betas menores)\n")
  } else if (avg_ratio > 1.1) {
    cat("    >>> 2SRR aplica MENOS shrinkage (betas maiores — menos parcimonioso)\n")
  } else {
    cat("    >>> Shrinkage similar entre os modelos\n")
  }
  
  write.csv(shrink_df, "results/parcimonia_shrinkage_h01.csv", row.names = FALSE)
}

# PARTE D: RESUMO PARA O ARTIGO

cat("  RESUMO: Outputs para o artigo\n")

cat("  Figuras geradas:\n")
cat("    1. betas_primeira_janela_h01.pdf — Trajetoria TVP vs Ridge (1a janela)\n")
cat("    2. betas_barplot_primeira_janela_h01.pdf — Barplot comparativo\n")
cat("    3. betas_oos_evolution_h01.pdf — Evolucao OOS: ultimo beta TVP vs Ridge\n")
cat("    4. parcimonia_hhi_h01.pdf — Concentracao de informacao ao longo do tempo\n")
cat("\n")
cat("  Tabelas geradas:\n")
cat("    1. betas_tvp_last_oos_h01.csv — Serie temporal dos ultimos betas TVP\n")
cat("    2. betas_ridge_oos_h01.csv — Serie temporal dos betas Ridge\n")
cat("    3. parcimonia_shrinkage_h01.csv — Shrinkage relativo por regressor\n")

cat("  11_betas_descriptive_analysis.R --- COMPLETO\n")
