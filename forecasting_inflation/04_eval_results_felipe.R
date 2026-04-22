# ============================================================
# 04_eval_results_felipe.R
#
# Avaliação dos resultados do 2SRR comparado com os demais
# modelos do repositório ForecastingInflation (Medeiros).
#
# Produz (todos salvos em results/):
#   1. rmsfe_comparativo.csv     — tabela RMSFE relativo ao RW
#   2. fig_rmsfe_comparativo.png — RMSFE por horizonte, todos modelos
#   3. fig_betas_2SRR_tempo.png  — betas time-varying do 2SRR (h=1)
#   4. fig_2srr_vs_ridge.png     — comparação direta 2SRR vs Ridge
#
# Pré-requisito:
#   - forecasts/rw.rda e forecasts/yout.rda (02_random_walk_oos_y.R)
#   - forecasts/2SRR.rda e forecasts/betas_2SRR.rda (03_call_model_felipe.R)
#   - forecasts dos outros modelos (03_call_model.R)
# ============================================================

#library(here)
#setwd(here("forecasting_inflation"))

setwd("~/tcc/Felipe_Dornelles_tcc/forecasting_inflation")

library(tidyverse)

source("functions/functions.R")
source("functions/felipe_functions.R")


# ============================================================
# 1. Carrega benchmark e realizados
# ============================================================

load("forecasts/yout.rda")
load("forecasts/rw.rda")


# ============================================================
# 2. Carrega todos os modelos da pasta forecasts/
# ============================================================

model_files <- setdiff(
  list.files("forecasts/"),
  c("rw.rda", "yout.rda", "betas_2SRR.rda")
)

models_list <- list()
for (i in seq_along(model_files)) {
  load(paste0("forecasts/", model_files[i]))
  models_list[[i]] <- forecasts
}
names(models_list) <- sub("\\.rda$", "", model_files)

cat("Modelos carregados:", paste(names(models_list), collapse = ", "), "\n\n")


# ============================================================
# 3. RMSFE — horizontes individuais h=1 a h=12
# ============================================================

rwe <- sqrt(colMeans((rw[, 1:12] - yout[, 1])^2, na.rm = TRUE))

errors <- lapply(models_list, function(x) {
  sqrt(colMeans((x[, 1:12] - yout[, 1])^2, na.rm = TRUE))
}) %>% Reduce(f = cbind)
colnames(errors) <- names(models_list)

res_h <- errors / rwe
rownames(res_h) <- paste0("h=", 1:12)


# ============================================================
# 4. RMSFE — horizontes acumulados (acc3, acc6, acc12)
# ============================================================

rweacc <- sqrt(colMeans((rw[, 13:15] - yout[, 2:4])^2, na.rm = TRUE))

errorsacc <- lapply(models_list, function(x) {
  sqrt(colMeans((x[, 13:15] - yout[, 2:4])^2, na.rm = TRUE))
}) %>% Reduce(f = cbind)
colnames(errorsacc) <- names(models_list)

res_acc <- errorsacc / rweacc
rownames(res_acc) <- c("acc3", "acc6", "acc12")

# ============================================================
# 5. Tabela final
# ============================================================

res_full <- rbind(res_h, res_acc)

cat("=== RMSFE relativo ao Random Walk ===\n")
print(round(res_full, 4))

if ("2SRR" %in% colnames(res_full)) {
  cat("\n=== 2SRR vs Random Walk ===\n")
  print(round(res_full[, "2SRR", drop = FALSE], 4))

  n_beat <- sum(res_full[1:12, "2SRR"] < 1.0, na.rm = TRUE)
  cat(sprintf("\n2SRR supera o Random Walk em %d de 12 horizontes.\n", n_beat))
}

write.csv(
  as.data.frame(round(res_full, 4)),
  file = "results/rmsfe_comparativo.csv"
)
cat("\nTabela salva em results/rmsfe_comparativo.csv\n")


# ============================================================
# 6. Gráfico comparativo de RMSFE por horizonte (h=1 a 12)
# ============================================================

df_plot_rmsfe <- as.data.frame(res_h) %>%
  mutate(horizonte = 1:12) %>%
  pivot_longer(-horizonte, names_to = "modelo", values_to = "rmsfe_rel") %>%
  mutate(
    destaque  = ifelse(modelo == "2SRR", "2SRR", "outros"),
    alpha_val = ifelse(modelo == "2SRR", 1.0, 0.4)
  )

p_rmsfe <- ggplot(df_plot_rmsfe,
       aes(x = horizonte, y = rmsfe_rel,
           colour    = modelo,
           group     = modelo,
           alpha     = alpha_val,
           linewidth = destaque)) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed",
             colour = "black", linewidth = 0.5) +
  scale_x_continuous(breaks = 1:12) +
  scale_alpha_identity() +
  scale_linewidth_manual(
    values = c("2SRR" = 1.2, "outros" = 0.5),
    guide  = "none"
  ) +
  labs(
    title    = "RMSFE relativo ao Random Walk — todos os modelos",
    subtitle = "Linha tracejada = benchmark (RW = 1). Abaixo de 1 = melhor que RW.",
    x        = "Horizonte (meses)",
    y        = "RMSFE / RMSFE(RW)",
    colour   = "Modelo"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold"),
    legend.position  = "right"
  )

print(p_rmsfe)
ggsave("results/fig_rmsfe_comparativo.png",
       p_rmsfe, width = 10, height = 5, dpi = 150)
cat("Gráfico RMSFE salvo em results/fig_rmsfe_comparativo.png\n")


# ============================================================
# 7. Gráfico de betas time-varying do 2SRR (horizonte h=1)
# ============================================================

load("forecasts/betas_2SRR.rda")

load("data/data.rda")
var_names <- setdiff(colnames(data), "date")

p_betas <- plot_betas_over_time(
  df_betas   = df_betas,
  var_names  = var_names,
  top_n      = 10,
  variable   = "CPIAUCSL",
  save_path  = "results/fig_betas_2SRR_tempo.png"
)

print(p_betas)
cat("Gráfico de betas salvo em results/fig_betas_2SRR_tempo.png\n")


# ============================================================
# 8. Comparação direta: 2SRR vs Ridge (se Ridge disponível)
# ============================================================

if (all(c("2SRR", "Ridge") %in% names(models_list))) {

  df_compare <- data.frame(
    horizonte   = 1:12,
    rmsfe_2SRR  = res_h[1:12, "2SRR"],
    rmsfe_Ridge = res_h[1:12, "Ridge"]
  ) %>%
    mutate(razao_2SRR_Ridge = rmsfe_2SRR / rmsfe_Ridge)

  cat("\n=== 2SRR vs Ridge (RMSFE relativo ao RW) ===\n")
  print(round(df_compare, 4))

  p_compare <- ggplot(df_compare, aes(x = horizonte)) +
    geom_line(aes(y = rmsfe_2SRR,  colour = "2SRR"),
              linewidth = 1.1) +
    geom_line(aes(y = rmsfe_Ridge, colour = "Ridge"),
              linewidth = 1.1, linetype = "dashed") +
    geom_hline(yintercept = 1, linetype = "dotted",
               colour = "grey40", linewidth = 0.4) +
    scale_x_continuous(breaks = 1:12) +
    scale_colour_manual(
      values = c("2SRR" = "#01696f", "Ridge" = "#964219")
    ) +
    labs(
      title    = "2SRR vs Ridge — RMSFE relativo ao Random Walk",
      subtitle = "Abaixo de 1 = melhor que RW. Abaixo da linha Ridge = 2SRR ganha.",
      x        = "Horizonte (meses)",
      y        = "RMSFE / RMSFE(RW)",
      colour   = "Modelo"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title       = element_text(face = "bold")
    )

  print(p_compare)
  ggsave("results/fig_2srr_vs_ridge.png",
         p_compare, width = 9, height = 5, dpi = 150)
  cat("Gráfico 2SRR vs Ridge salvo em results/fig_2srr_vs_ridge.png\n")
}

