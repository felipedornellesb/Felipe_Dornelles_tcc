# ============================================================
# 09_tvp_analysis.R 
#
# | Caso | Especificação | Viabilidade | O que demonstra |
# |------|--------------|-------------|-----------------|
# | 1 | TVP-AR (só lags y) | ✅ Fácil, rápido | Valor puro do TVP sem informação extra |
# | 2 | TVP-Factor (só PCA) | ✅ Fácil | Valor dos fatores sem persistência AR |
# | 3 | TVP-FAVAR (fatores + y) | ✅ **06_coulombe_2SRR_pipeline.R** | Modelo completo (benchmark) |
# | 3+ | TVP + 116 vars brutas | ❌ Precisa shrinkTVP | Não viável com TVPRR_cosso |
# Compara Caso 1, 2 e 3
# ============================================================

rm(list = ls())
setwd("~/TCC/tcc/forecasting_inflation")

library(forecast)
library(lmtest)
library(sandwich)
library(ggplot2)
library(reshape2)

hor <- c(1, 3, 6, 12)
case_names <- c("TVP_AR", "TVP_Factor", "TVP_FAVAR")

# Arquivos de cada caso
case_files <- list(
  TVP_AR     = "forecasts/tvp_TVP_AR_h%02d.csv",
  TVP_Factor = "forecasts/tvp_TVP_Factor_h%02d.csv",
  TVP_FAVAR  = "forecasts/coulombe_fc_h%02d.csv"  # Caso 3 ja rodado
)

cat("  COMPARACAO TVP: AR vs Factor vs FAVAR\n")

dir.create("results", showWarnings = FALSE)
dir.create("results/figures", showWarnings = FALSE)

# ============================================================
# 1. TABELA COMPARATIVA DE RMSE
# ============================================================
cat("=== 1. RMSE por caso e horizonte ===\n\n")

all_results <- list()

for (h in hor) {
  cat(sprintf("  h = %d:\n", h))
  cat(sprintf("  %-15s %10s %10s %10s %10s\n",
              "Caso", "RMSE_Ridge", "RMSE_2SRR", "Ratio", "Melhor"))
  cat(paste(rep("-", 60), collapse = ""), "\n")

  for (cn in case_names) {
    fname <- sprintf(case_files[[cn]], h)
    if (!file.exists(fname)) {
      cat(sprintf("  %-15s [arquivo nao encontrado]\n", cn))
      next
    }

    df <- read.csv(fname, stringsAsFactors = FALSE)
    df <- df[complete.cases(df$realized, df$fc_ridge, df$fc_2srr), ]
    if (nrow(df) < 20) next

    rmse_r <- sqrt(mean((df$fc_ridge - df$realized)^2))
    rmse_2 <- sqrt(mean((df$fc_2srr  - df$realized)^2))
    ratio  <- rmse_2 / rmse_r
    melhor <- ifelse(ratio < 0.95, "2SRR", ifelse(ratio > 1.05, "Ridge", "Empate"))

    all_results[[paste0(cn, "_h", h)]] <- data.frame(
      h = h, case = cn,
      RMSE_Ridge = rmse_r, RMSE_2SRR = rmse_2,
      Ratio = ratio, melhor = melhor)

    cat(sprintf("  %-15s %10.4f %10.4f %10.4f %10s\n",
                cn, rmse_r, rmse_2, ratio, melhor))
  }
  cat("\n")
}

# Tabela consolidada
if (length(all_results) > 0) {
  tab <- do.call(rbind, all_results)
  write.csv(tab, "results/tvp_comparison_table.csv", row.names = FALSE)
  cat("Tabela salva: results/tvp_comparison_table.csv\n\n")
}

# ============================================================
# 2. DM TEST ENTRE CASOS (2SRR de cada caso)
# ============================================================
cat("=== 2. Diebold-Mariano entre especificacoes 2SRR ===\n\n")

for (h in hor) {
  cat(sprintf("  h = %d:\n", h))

  # Carrega erros de cada caso
  errors <- list()
  realized <- NULL

  for (cn in case_names) {
    fname <- sprintf(case_files[[cn]], h)
    if (!file.exists(fname)) next
    df <- read.csv(fname, stringsAsFactors = FALSE)
    df <- df[complete.cases(df$realized, df$fc_2srr), ]
    if (nrow(df) < 20) next
    errors[[cn]] <- df$fc_2srr - df$realized
    if (is.null(realized)) realized <- df$realized
  }

  # Pairwise DM
  cn_list <- names(errors)
  for (i in seq_along(cn_list)) {
    for (j in seq_along(cn_list)) {
      if (i >= j) next
      cn1 <- cn_list[i]; cn2 <- cn_list[j]

      # Alinhar tamanhos
      n_min <- min(length(errors[[cn1]]), length(errors[[cn2]]))
      e1 <- errors[[cn1]][1:n_min]
      e2 <- errors[[cn2]][1:n_min]

      dm <- tryCatch(
        dm.test(e1, e2, h = h, alternative = "two.sided", power = 2),
        error = function(e) NULL)

      if (!is.null(dm)) {
        sig <- ifelse(dm$p.value < 0.01, "***",
               ifelse(dm$p.value < 0.05, "**",
               ifelse(dm$p.value < 0.10, "*", "")))
        melhor <- ifelse(dm$statistic > 0, cn2, cn1)
        cat(sprintf("    %s vs %s: DM=%.3f p=%.4f %s -> %s melhor\n",
                    cn1, cn2, dm$statistic, dm$p.value, sig, melhor))
      }
    }
  }
  cat("\n")
}

# ============================================================
# 3. GRAFICO COMPARATIVO DE RMSE
# ============================================================
cat("=== 3. Graficos ===\n")

if (length(all_results) > 0) {
  tab <- do.call(rbind, all_results)

  # RMSE do 2SRR por caso e horizonte
  p1 <- ggplot(tab, aes(x = factor(h), y = RMSE_2SRR, fill = case)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = "RMSE do 2SRR por especificacao e horizonte",
         x = "Horizonte (meses)", y = "RMSE", fill = "Especificacao") +
    theme_minimal()
  ggsave("results/figures/tvp_comparison_rmse.pdf", p1, width = 10, height = 6)

  # Ratio por caso e horizonte
  p2 <- ggplot(tab, aes(x = factor(h), y = Ratio, fill = case)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    labs(title = "Ratio RMSE(2SRR/Ridge) por especificacao",
         subtitle = "Abaixo de 1 = 2SRR melhor",
         x = "Horizonte (meses)", y = "Ratio", fill = "Especificacao") +
    theme_minimal()
  ggsave("results/figures/tvp_comparison_ratio.pdf", p2, width = 10, height = 6)

  cat("  Graficos salvos em results/figures/\n")
}

cat("  09_tvp_analysis.R --- COMPLETO\n")
