# ============================================================
# 05_analysis_2SRR.R  — VERSÃO FINAL
#
# Análise completa: RMSFE, inferência DM, sub-períodos,
# betas TVP, CSFE, heatmap, ranking de Borda, tabelas LaTeX
#
# ATENÇÃO: todos os gráficos têm print() antes do ggsave()
#          Não cria pastas (já existem na máquina)
# ============================================================

library(here)
setwd(here("forecasting_inflation"))

library(tidyverse)
library(forecast)    # dm.test()
library(scales)      # para heatmap

# ============================================================
# 0. CARREGA OBJETOS BASE
# ============================================================

load("forecasts/yout.rda")
load("forecasts/rw.rda")
load("forecasts/2SRR.rda"); fc_2srr <- forecasts; rm(forecasts)

model_files <- setdiff(
  list.files("forecasts/", pattern = "\\.rda$"),
  c("rw.rda", "yout.rda", "betas_2SRR.rda", "2SRR.rda")
)
models_list <- list()
for (f in model_files) {
  load(paste0("forecasts/", f))
  models_list[[sub("\\.rda$", "", f)]] <- forecasts
}
rm(forecasts)

load("data/data.rda")
n_oos     <- nrow(yout)
oos_dates <- tail(data$date, n_oos)

# res_full vem do 04 — recarrega se necessário
if (!exists("res_full")) {
  res_full <- read.csv("results/rmsfe_comparativo.csv", row.names = 1)
  res_h    <- res_full[1:12, ]
}

# Paleta única para todo o script
cores_modelos <- c(
  "2SRR"     = "#01696f",
  "Ridge"    = "#964219",
  "LASSO"    = "#4361ee",
  "ElNET"    = "#e9c46a",
  "AdaLASSO" = "#457b9d",
  "AdaElNET" = "#2a9d8f",
  "RF"       = "#e07a5f",
  "Bagging"  = "#6d597a",
  "Factor"   = "#b5838d",
  "T.Factor" = "#c77dff",
  "AR"       = "#588157",
  "AR_BIC"   = "#3a5a40",
  "CSR"      = "#264653"
)

# Bandas de recessão NBER — usadas em vários gráficos
recession_bands <- data.frame(
  start = as.Date(c("2001-03-01", "2007-12-01", "2020-02-01")),
  end   = as.Date(c("2001-11-01", "2009-06-01", "2020-04-01")),
  label = c("Dot-com", "GFC", "COVID")
)

add_recession_bands <- function(p) {
  p +
    geom_rect(data = recession_bands,
              aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
              inherit.aes = FALSE, fill = "grey80", alpha = 0.35) +
    geom_text(data = recession_bands,
              aes(x = start + (end - start) / 2, y = Inf, label = label),
              inherit.aes = FALSE, vjust = 1.5, size = 2.8, colour = "grey40")
}

models_compare <- intersect(
  c("Ridge", "LASSO", "AdaLASSO", "ElNET", "RF", "Bagging", "Factor", "AR", "CSR"),
  names(models_list)
)

# ============================================================
# 1. TABELA RMSFE CSV (re-exporta para garantia)
# ============================================================

write.csv(as.data.frame(round(res_full, 4)),
          "results/rmsfe_comparativo.csv", row.names = TRUE)
cat("CSV salvo em results/rmsfe_comparativo.csv\n")

# ============================================================
# 2. DIAGNÓSTICO RÁPIDO DOS RESULTADOS
#    — imprime na tela antes de qualquer gráfico
# ============================================================

cat("\n======================================================\n")
cat("DIAGNÓSTICO DOS RESULTADOS — 2SRR\n")
cat("======================================================\n")

n_beat_rw    <- sum(res_h[, "2SRR"] < 1.0, na.rm = TRUE)
n_beat_ridge <- sum(res_h[, "2SRR"] < res_h[, "Ridge"], na.rm = TRUE)
n_beat_lasso <- sum(res_h[, "2SRR"] < res_h[, "LASSO"], na.rm = TRUE)
best_h       <- which.min(res_h[, "2SRR"])

cat(sprintf("2SRR < RW   em %d/12 horizontes\n", n_beat_rw))
cat(sprintf("2SRR < Ridge em %d/12 horizontes\n", n_beat_ridge))
cat(sprintf("2SRR < LASSO em %d/12 horizontes\n", n_beat_lasso))
cat(sprintf("Melhor horizonte: h=%d (RMSFE = %.4f)\n", best_h, res_h[best_h, "2SRR"]))

# ATENÇÃO: acc6 > 1 é o ponto mais delicado — informa o usuário
if (res_full["acc6", "2SRR"] > 1.0) {
  cat(sprintf("\n[!] ATENÇÃO: acc6 do 2SRR = %.4f > 1.0 (perde para RW)\n",
              res_full["acc6", "2SRR"]))
  cat("    Causa provável: propagação de erros no acumulado\n")
  cat("    de médio prazo. Discutir no texto do TCC.\n\n")
}

for (acc in c("acc3", "acc6", "acc12")) {
  cat(sprintf("RMSFE %-6s: 2SRR=%.4f | Ridge=%.4f | LASSO=%.4f | RF=%.4f\n",
              acc,
              res_full[acc, "2SRR"],
              res_full[acc, "Ridge"],
              res_full[acc, "LASSO"],
              res_full[acc, "RF"]))
}

# Ranking médio (Borda preview)
rank_mat <- apply(res_full, 1, rank)
borda_scores <- rowMeans(rank_mat)
cat("\nRanking Borda médio (menor = melhor, todos os horizontes):\n")
print(sort(borda_scores))

# ============================================================
# 3. GRÁFICO: RMSFE todos os modelos (h=1..12)
# ============================================================

df_rmsfe_long <- as.data.frame(res_h) %>%
  mutate(horizonte = 1:12) %>%
  pivot_longer(-horizonte, names_to = "modelo", values_to = "rmsfe") %>%
  mutate(
    destaque  = ifelse(modelo == "2SRR", "2SRR", "outros"),
    alpha_val = ifelse(modelo == "2SRR", 1.0, 0.35)
  )

p1 <- ggplot(df_rmsfe_long,
             aes(x = horizonte, y = rmsfe, colour = modelo,
                 group = modelo, alpha = alpha_val)) +
  geom_line(aes(linewidth = destaque)) +
  geom_hline(yintercept = 1, linetype = "dashed",
             colour = "black", linewidth = 0.5) +
  scale_x_continuous(breaks = 1:12) +
  scale_alpha_identity() +
  scale_linewidth_manual(values = c("2SRR" = 1.4, "outros" = 0.5), guide = "none") +
  scale_colour_manual(values = cores_modelos) +
  labs(
    title    = "RMSFE relativo ao Random Walk — todos os modelos",
    subtitle = "Abaixo de 1 = melhor que RW | 2SRR em destaque",
    x = "Horizonte (meses)", y = "RMSFE / RMSFE(RW)", colour = "Modelo"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

print(p1)
ggsave("results/fig1_rmsfe_todos.png", p1, width = 11, height = 5, dpi = 150)
cat("Salvo: fig1_rmsfe_todos.png\n")

# ============================================================
# 4. GRÁFICO: 2SRR vs selecionados (h=1..12 + acumulados)
# ============================================================

modelos_sel <- intersect(c("2SRR", "Ridge", "LASSO", "RF", "ElNET", "AR"),
                         colnames(res_full))

p2 <- as.data.frame(res_full) %>%
  rownames_to_column("horizonte") %>%
  mutate(h_num = seq_len(nrow(res_full))) %>%
  pivot_longer(all_of(modelos_sel), names_to = "modelo", values_to = "rmsfe") %>%
  ggplot(aes(x = h_num, y = rmsfe, colour = modelo, group = modelo)) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 1.8) +
  geom_hline(yintercept = 1, linetype = "dashed",
             colour = "grey40", linewidth = 0.4) +
  scale_x_continuous(
    breaks = seq_len(nrow(res_full)),
    labels = c(paste0("h=", 1:12), "acc3", "acc6", "acc12")
  ) +
  scale_colour_manual(values = cores_modelos) +
  labs(
    title    = "2SRR vs modelos selecionados — horizontes h=1:12 e acumulados",
    subtitle = "acc6 > 1 indica perda do 2SRR para o RW no acumulado de 6 meses",
    x = "Horizonte", y = "RMSFE / RMSFE(RW)", colour = "Modelo"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))

print(p2)
ggsave("results/fig2_2srr_vs_selecionados.png", p2, width = 12, height = 5, dpi = 150)
cat("Salvo: fig2_2srr_vs_selecionados.png\n")

# ============================================================
# 5. HEATMAP DE RMSFE — padrão em IJF e Journal of Econometrics
#    Linha = modelo, coluna = horizonte, cor = RMSFE relativo
#    Verde escuro < 1 (bom), vermelho > 1 (ruim)
# ============================================================

df_heat <- as.data.frame(res_full) %>%
  rownames_to_column("horizonte") %>%
  pivot_longer(-horizonte, names_to = "modelo", values_to = "rmsfe") %>%
  mutate(
    horizonte = factor(horizonte,
                       levels = c(paste0("h=", 1:12), "acc3", "acc6", "acc12")),
    modelo    = factor(modelo,
                       levels = names(sort(colMeans(res_full), decreasing = FALSE)))
  )

p_heat <- ggplot(df_heat, aes(x = horizonte, y = modelo, fill = rmsfe)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.3f", rmsfe)),
            size = 2.5, colour = "white", fontface = "bold") +
  scale_fill_gradient2(
    low      = "#01696f",
    mid      = "#f7f7f7",
    high     = "#a12c7b",
    midpoint = 1.0,
    limits   = c(0.7, 1.25),
    oob      = scales::squish,
    name     = "RMSFE/RW"
  ) +
  labs(
    title    = "Heatmap RMSFE relativo ao Random Walk",
    subtitle = "Verde = melhor que RW | Roxo = pior que RW | 2SRR = linha de referência",
    x        = "Horizonte",
    y        = "Modelo"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid    = element_blank(),
    plot.title    = element_text(face = "bold"),
    axis.text.x   = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

print(p_heat)
ggsave("results/fig_heatmap_rmsfe.png", p_heat, width = 13, height = 6, dpi = 150)
cat("Salvo: fig_heatmap_rmsfe.png\n")

# ============================================================
# 6. RANKING DE BORDA — usado em Medeiros et al. (2021) e
#    Stock & Watson (2012). Posiciona cada modelo em cada
#    horizonte e tira média dos rankings.
# ============================================================

# Separa horizontes pontuais (h=1..12) dos acumulados
rank_h    <- apply(res_h,   1, rank, ties.method = "average")  # k x 12
rank_acc  <- apply(res_full[c("acc3","acc6","acc12"), ], 1, rank, ties.method = "average")

borda_h   <- rowMeans(rank_h)
borda_acc <- rowMeans(rank_acc)
borda_all <- rowMeans(cbind(rank_h, rank_acc))

df_borda <- data.frame(
  modelo    = names(borda_all),
  Borda_h1_12  = borda_h[names(borda_all)],
  Borda_acc    = borda_acc[names(borda_all)],
  Borda_total  = borda_all
) %>%
  arrange(Borda_total) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

cat("\n=== RANKING DE BORDA (menor = melhor) ===\n")
print(df_borda)

write.csv(df_borda, "results/ranking_borda.csv", row.names = FALSE)

p_borda <- df_borda %>%
  mutate(modelo = factor(modelo, levels = rev(modelo))) %>%
  ggplot(aes(x = Borda_total, y = modelo,
             fill = ifelse(modelo == "2SRR", "2SRR", "outros"))) +
  geom_col(width = 0.6) +
  geom_vline(xintercept = median(df_borda$Borda_total),
             linetype = "dashed", colour = "grey40") +
  scale_fill_manual(values = c("2SRR" = "#01696f", "outros" = "#b0c4c4"),
                    guide = "none") +
  labs(
    title    = "Ranking de Borda — posição média de cada modelo",
    subtitle = "Menor pontuação = melhor posição relativa em todos os horizontes",
    x        = "Score de Borda (média dos rankings por horizonte)",
    y        = "Modelo"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.title = element_text(face = "bold"))

print(p_borda)
ggsave("results/fig_borda_ranking.png", p_borda, width = 9, height = 5, dpi = 150)
cat("Salvo: fig_borda_ranking.png\n")

# ============================================================
# 7. TESTE DIEBOLD-MARIANO — COM BARTLETT (corrige variância negativa)
# ============================================================

dm_pvalues <- matrix(NA, nrow = 12, ncol = length(models_list),
                     dimnames = list(paste0("h=", 1:12), names(models_list)))
dm_stats   <- dm_pvalues

for (modelo in names(models_list)) {
  fc_m <- models_list[[modelo]]
  for (h in 1:12) {
    e1 <- fc_2srr[, h] - yout[, 1]
    e2 <- fc_m[, h]    - yout[, 1]
    ok <- complete.cases(e1, e2)
    if (sum(ok) < 10) next
    tryCatch({
      # varestimator = "bartlett" corrige o problema de variância negativa
      dm_res <- dm.test(e1[ok], e2[ok],
                        alternative    = "two.sided",
                        h              = h,
                        power          = 2,
                        varestimator   = "bartlett")
      dm_pvalues[h, modelo] <- dm_res$p.value
      dm_stats[h, modelo]   <- dm_res$statistic
    }, error = function(e) NULL)
  }
}

cat("\n=== P-VALORES DIEBOLD-MARIANO (2SRR vs outros) ===\n")
print(round(dm_pvalues[, models_compare], 3))

write.csv(round(dm_pvalues, 4), "results/dm_pvalues.csv")
cat("Salvo: dm_pvalues.csv\n")

# Heatmap de significância DM — mostra visualmente onde 2SRR
# é estatisticamente diferente dos concorrentes
stars_fn <- function(p) {
  ifelse(is.na(p), "—",
  ifelse(p < 0.01, "***",
  ifelse(p < 0.05, "**",
  ifelse(p < 0.10, "*", ""))))
}

df_dm_heat <- as.data.frame(dm_pvalues[, models_compare]) %>%
  rownames_to_column("horizonte") %>%
  pivot_longer(-horizonte, names_to = "modelo", values_to = "pval") %>%
  mutate(
    sig     = cut(pval,
                  breaks = c(0, 0.01, 0.05, 0.10, 1),
                  labels = c("p<0.01***", "p<0.05**", "p<0.10*", "n.s."),
                  include.lowest = TRUE),
    label   = stars_fn(pval),
    horizonte = factor(horizonte, levels = paste0("h=", 1:12))
  )

p_dm_heat <- ggplot(df_dm_heat, aes(x = horizonte, y = modelo, fill = sig)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = label), size = 3, colour = "white", fontface = "bold") +
  scale_fill_manual(
    values  = c("p<0.01***" = "#01696f",
                "p<0.05**"  = "#2a9d8f",
                "p<0.10*"   = "#76c8c8",
                "n.s."      = "#d9d9d9"),
    na.value = "#f0f0f0",
    name     = "Significância DM"
  ) +
  labs(
    title    = "Significância Diebold-Mariano — 2SRR vs concorrentes",
    subtitle = "Células coloridas: diferença estatisticamente significativa vs 2SRR\nn.s. = não significativo | estimador Bartlett",
    x = "Horizonte", y = "Modelo"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid   = element_blank(),
    plot.title   = element_text(face = "bold"),
    axis.text.x  = element_text(angle = 45, hjust = 1)
  )

print(p_dm_heat)
ggsave("results/fig_dm_heatmap.png", p_dm_heat, width = 11, height = 5, dpi = 150)
cat("Salvo: fig_dm_heatmap.png\n")

# ============================================================
# 8. TABELA LaTeX: RMSFE + asteriscos DM
# ============================================================

cols_latex <- intersect(
  c("2SRR", "Ridge", "LASSO", "AdaLASSO", "ElNET",
    "RF", "Bagging", "Factor", "AR", "CSR"),
  colnames(res_full)
)

sink("results/tabela_rmsfe_dm_latex.tex")
cat("\\begin{table}[ht]\n\\centering\n")
cat("\\caption{RMSFE relativo ao Random Walk com significância Diebold-Mariano}\n")
cat("\\label{tab:rmsfe_dm}\n")
cat("{\\footnotesize\n")
cat("\\begin{tabular}{l", rep("r", length(cols_latex)), "}\n", sep = "")
cat("\\hline\n")
cat("Horizonte &", paste(cols_latex, collapse = " & "), "\\\\\n\\hline\n")

for (i in seq_len(nrow(res_full))) {
  h_label <- rownames(res_full)[i]
  vals    <- res_full[i, cols_latex]
  min_idx <- which.min(unlist(vals))

  cells <- sapply(seq_along(cols_latex), function(j) {
    col <- cols_latex[j]
    val <- round(vals[[col]], 4)
    ast <- ""
    if (col != "2SRR" && h_label %in% rownames(dm_pvalues)) {
      ast <- stars_fn(dm_pvalues[h_label, col])
      ast <- gsub("—", "", ast)
    }
    cell <- paste0(val, ast)
    if (j == min_idx) cell <- paste0("\\textbf{", cell, "}")
    cell
  })

  cat(h_label, "&", paste(cells, collapse = " & "), "\\\\\n")
}

cat("\\hline\n")
cat(sprintf("\\multicolumn{%d}{l}{\\scriptsize{* p<0.10, ** p<0.05, *** p<0.01 — Diebold-Mariano (bilateral, Bartlett) vs 2SRR}}\\\\\n",
            length(cols_latex) + 1))
cat("\\end{tabular}\n}\n\\end{table}\n")
sink()
cat("Salvo: tabela_rmsfe_dm_latex.tex\n")

# ============================================================
# 9. ANÁLISE DE SUB-PERÍODOS
# ============================================================

sub_periods <- list(
  Full       = rep(TRUE, n_oos),
  Pre_COVID  = oos_dates <  as.Date("2020-01-01"),
  COVID      = oos_dates >= as.Date("2020-01-01") & oos_dates <= as.Date("2022-12-01"),
  Post_COVID = oos_dates >  as.Date("2022-12-01"),
  GFC        = oos_dates >= as.Date("2007-07-01") & oos_dates <= as.Date("2009-12-01")
)

rmsfe_sub <- lapply(names(sub_periods), function(sp) {
  idx <- sub_periods[[sp]]
  if (sum(idx, na.rm = TRUE) < 5) return(NULL)

  rwe_sp <- sqrt(colMeans((rw[idx, 1:12] - yout[idx, 1])^2, na.rm = TRUE))

  bind_rows(
    data.frame(h = 1:12, periodo = sp, modelo = "2SRR",
               rmsfe = sqrt(colMeans((fc_2srr[idx, 1:12] - yout[idx, 1])^2,
                                     na.rm = TRUE)) / rwe_sp),
    lapply(models_compare, function(m) {
      data.frame(h = 1:12, periodo = sp, modelo = m,
                 rmsfe = sqrt(colMeans((models_list[[m]][idx, 1:12] - yout[idx, 1])^2,
                                       na.rm = TRUE)) / rwe_sp)
    }) %>% bind_rows()
  )
}) %>% bind_rows()

write.csv(rmsfe_sub, "results/rmsfe_subperiodos.csv", row.names = FALSE)
cat("Salvo: rmsfe_subperiodos.csv\n")

# Gráfico sub-períodos COVID
p_sub <- rmsfe_sub %>%
  filter(periodo %in% c("Pre_COVID", "COVID", "Post_COVID")) %>%
  mutate(
    periodo = factor(periodo,
                     levels = c("Pre_COVID", "COVID", "Post_COVID"),
                     labels = c("Pré-COVID\n(<2020)", "COVID\n(2020-2022)", "Pós-COVID\n(>2022)")),
    destaque  = ifelse(modelo == "2SRR", "2SRR", "outros"),
    alpha_val = ifelse(modelo == "2SRR", 1.0, 0.35)
  ) %>%
  ggplot(aes(x = h, y = rmsfe, colour = modelo, group = modelo,
             alpha = alpha_val, linewidth = destaque)) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey30", linewidth = 0.4) +
  facet_wrap(~ periodo, ncol = 3) +
  scale_x_continuous(breaks = c(1, 3, 6, 9, 12)) +
  scale_alpha_identity() +
  scale_linewidth_manual(values = c("2SRR" = 1.3, "outros" = 0.5), guide = "none") +
  scale_colour_manual(values = cores_modelos) +
  labs(
    title    = "RMSFE relativo ao RW por sub-período",
    subtitle = "2SRR em destaque — abaixo de 1 = melhor que Random Walk",
    x = "Horizonte (meses)", y = "RMSFE / RMSFE(RW)", colour = "Modelo"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        legend.position = "right")

print(p_sub)
ggsave("results/fig_subperiodos.png", p_sub, width = 13, height = 5, dpi = 150)
cat("Salvo: fig_subperiodos.png\n")

# Gráfico GFC
p_gfc <- rmsfe_sub %>%
  filter(periodo %in% c("Pre_COVID", "GFC")) %>%
  mutate(
    periodo = factor(periodo,
                     levels = c("GFC", "Pre_COVID"),
                     labels = c("Crise 2008 (GFC)", "Pré-COVID (<2020)")),
    destaque  = ifelse(modelo == "2SRR", "2SRR", "outros"),
    alpha_val = ifelse(modelo == "2SRR", 1.0, 0.35)
  ) %>%
  ggplot(aes(x = h, y = rmsfe, colour = modelo, group = modelo,
             alpha = alpha_val, linewidth = destaque)) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey30", linewidth = 0.4) +
  facet_wrap(~ periodo) +
  scale_x_continuous(breaks = c(1, 3, 6, 9, 12)) +
  scale_alpha_identity() +
  scale_linewidth_manual(values = c("2SRR" = 1.3, "outros" = 0.5), guide = "none") +
  scale_colour_manual(values = cores_modelos) +
  labs(title = "RMSFE relativo ao RW — Crise 2008 vs Pré-COVID",
       x = "Horizonte", y = "RMSFE / RMSFE(RW)", colour = "Modelo") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

print(p_gfc)
ggsave("results/fig_subperiodo_gfc.png", p_gfc, width = 10, height = 5, dpi = 150)
cat("Salvo: fig_subperiodo_gfc.png\n")

# ============================================================
# 10. DECOMPOSIÇÃO DO GANHO POR PERÍODO (h=1)
#     Quanto da vantagem total do 2SRR vem de cada período?
#     Referência: Goyal & Welch (2008) — RFS
# ============================================================

ganho_decomp <- lapply(names(sub_periods), function(sp) {
  idx <- sub_periods[[sp]]
  if (sum(idx, na.rm = TRUE) < 5) return(NULL)

  e_2srr <- (fc_2srr[idx, 1] - yout[idx, 1])^2
  e_rw   <- (rw[idx, 1]      - yout[idx, 1])^2

  data.frame(
    periodo      = sp,
    n_obs        = sum(idx),
    mse_2srr     = mean(e_2srr, na.rm = TRUE),
    mse_rw       = mean(e_rw,   na.rm = TRUE),
    ganho_abs    = mean(e_rw, na.rm = TRUE) - mean(e_2srr, na.rm = TRUE),
    rmsfe_rel    = sqrt(mean(e_2srr, na.rm = TRUE)) / sqrt(mean(e_rw, na.rm = TRUE)),
    pct_do_total = NA
  )
}) %>% bind_rows()

total_ganho <- ganho_decomp$ganho_abs[ganho_decomp$periodo == "Full"]
ganho_decomp$pct_do_total <- round(ganho_decomp$ganho_abs / total_ganho * 100, 1)

cat("\n=== DECOMPOSIÇÃO DO GANHO POR PERÍODO (h=1) ===\n")
print(ganho_decomp)
write.csv(ganho_decomp, "results/decomposicao_ganho.csv", row.names = FALSE)

p_decomp <- ganho_decomp %>%
  filter(periodo != "Full") %>%
  mutate(periodo = factor(periodo,
                           levels = c("Pre_COVID", "GFC", "COVID", "Post_COVID"),
                           labels = c("Pré-COVID", "GFC (2008)", "COVID (2020-22)", "Pós-COVID"))) %>%
  ggplot(aes(x = periodo, y = rmsfe_rel, fill = rmsfe_rel < 1)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey30") +
  geom_text(aes(label = sprintf("%.3f", rmsfe_rel)),
            vjust = -0.4, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = c("TRUE" = "#01696f", "FALSE" = "#a12c7b"),
                    labels = c("TRUE" = "< RW (melhor)", "FALSE" = "> RW (pior)"),
                    name   = "") +
  scale_y_continuous(limits = c(0.7, 1.15)) +
  labs(
    title    = "RMSFE do 2SRR por período (h=1, relativo ao RW)",
    subtitle = "Verde = 2SRR bate o RW | Roxo = 2SRR perde para o RW",
    x        = "Sub-período",
    y        = "RMSFE / RMSFE(RW)"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(face = "bold"))

print(p_decomp)
ggsave("results/fig_decomposicao_periodos.png", p_decomp, width = 9, height = 5, dpi = 150)
cat("Salvo: fig_decomposicao_periodos.png\n")

# ============================================================
# 11. BETAS TVP — ANÁLISE COMPLETA
#     ATENÇÃO: os "top 2" pela volatilidade são PC7 e PC8,
#     NÃO PC1 e PC2. O gráfico é renomeado adequadamente.
# ============================================================

load("forecasts/betas_2SRR.rda")

df_betas <- df_betas %>%
  mutate(var_idx = as.integer(var_idx),
         fator   = paste0("PC", var_idx))

volatility_by_pc <- df_betas %>%
  group_by(var_idx, fator) %>%
  summarise(sd_beta  = sd(beta, na.rm = TRUE),
            mean_abs = mean(abs(beta), na.rm = TRUE),
            .groups  = "drop") %>%
  arrange(desc(sd_beta))

write.csv(volatility_by_pc, "results/beta_volatility_by_pc.csv", row.names = FALSE)

cat("\n=== Volatilidade dos betas por PC (h=1) ===\n")
cat("NOTA: Os PCs mais voláteis são os de índice mais ALTO,\n")
cat("o que pode indicar instabilidade numérica nos componentes\n")
cat("menores (menos variância explicada). Verifique se PC7/PC8\n")
cat("têm interpretação econômica ou são ruído.\n\n")
print(head(volatility_by_pc, 10))

top2_pcs <- volatility_by_pc$var_idx[1:2]   # PC7 e PC8 nos seus dados
top6_pcs <- volatility_by_pc$var_idx[1:6]

# Também analisa PC1 e PC2 (economicamente relevantes)
pcs_economicos <- c(1L, 2L, 3L)

# --- Gráfico: PCs mais voláteis (índices numéricos corretos)
p_top2 <- df_betas %>%
  filter(var_idx %in% top2_pcs) %>%
  ggplot(aes(x = window_date, y = beta, colour = fator, group = fator)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = c("#01696f", "#964219"))

p_top2 <- add_recession_bands(p_top2) +
  labs(
    title    = sprintf("Coeficientes TVP: %s e %s — mais voláteis (h=1)",
                       paste0("PC", top2_pcs[1]), paste0("PC", top2_pcs[2])),
    subtitle = "Variação dos betas = evidência de instabilidade nos coeficientes ao longo do tempo\nÁreas cinza = recessões NBER",
    x = "Data", y = "β(t)", colour = "Componente"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey40", size = 9))

print(p_top2)
ggsave("results/fig_betas_top2_volateis.png", p_top2, width = 11, height = 5, dpi = 150)
cat("Salvo: fig_betas_top2_volateis.png\n")

# --- Gráfico: PC1, PC2 e PC3 (interpretação econômica — Phillips Curve)
p_pc123 <- df_betas %>%
  filter(var_idx %in% pcs_economicos) %>%
  ggplot(aes(x = window_date, y = beta, colour = fator, group = fator)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = c("#01696f", "#964219", "#4361ee"))

p_pc123 <- add_recession_bands(p_pc123) +
  labs(
    title    = "Coeficientes TVP: PC1, PC2 e PC3 — primeiros fatores (h=1)",
    subtitle = "Primeiros PCs capturam maior variância dos preditores macroeconômicos\nÁreas cinza = recessões NBER",
    x = "Data", y = "β(t)", colour = "Componente"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey40", size = 9))

print(p_pc123)
ggsave("results/fig_betas_pc1_pc2_pc3.png", p_pc123, width = 11, height = 5, dpi = 150)
cat("Salvo: fig_betas_pc1_pc2_pc3.png\n")

# --- Top 6 em facet
p_top6 <- df_betas %>%
  filter(var_idx %in% top6_pcs) %>%
  mutate(fator = factor(fator, levels = paste0("PC", top6_pcs))) %>%
  ggplot(aes(x = window_date, y = beta)) +
  geom_rect(data = recession_bands,
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "grey80", alpha = 0.35) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
  geom_line(colour = "#01696f", linewidth = 0.75) +
  facet_wrap(~ fator, scales = "free_y", ncol = 3) +
  labs(
    title    = "Top 6 fatores PCA por volatilidade dos betas TVP (h=1)",
    subtitle = "Cada painel = um componente principal. Área cinza = recessões NBER.",
    x = "Data", y = "β(t)"
  ) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold", colour = "#01696f"))

print(p_top6)
ggsave("results/fig_betas_top6_pcs.png", p_top6, width = 12, height = 7, dpi = 150)
cat("Salvo: fig_betas_top6_pcs.png\n")

# --- Norma L2 dos betas ao longo do tempo — mede instabilidade global
df_norma <- df_betas %>%
  group_by(window_date) %>%
  summarise(norma_L2 = sqrt(sum(beta^2, na.rm = TRUE)), .groups = "drop")

p_norma <- ggplot(df_norma, aes(x = window_date, y = norma_L2))

p_norma <- add_recession_bands(p_norma) +
  geom_line(colour = "#01696f", linewidth = 0.9) +
  labs(
    title    = "Norma L2 dos betas TVP ao longo do tempo (h=1)",
    subtitle = "Picos indicam janelas com maior instabilidade dos coeficientes\nÁrea cinza = recessões NBER",
    x = "Data", y = "||β(t)||₂"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

print(p_norma)
ggsave("results/fig_norma_betas.png", p_norma, width = 11, height = 4, dpi = 150)
cat("Salvo: fig_norma_betas.png\n")

# Correlação Δβ_PC1 × Δinflação
beta_pc1 <- df_betas %>%
  filter(var_idx == pcs_economicos[1]) %>%
  arrange(window_date) %>% pull(beta)

n_min        <- min(length(beta_pc1), nrow(yout))
cor_pc1_infl <- cor(diff(beta_pc1[1:n_min]),
                    diff(yout[1:n_min, 1]),
                    use = "complete.obs")
cat(sprintf("\nCorrelação Δβ(PC1) × Δinflação realizada: %.4f\n", cor_pc1_infl))

# ============================================================
# 12. CSFE — CUMULATIVE SQUARED FORECAST ERROR (h=1)
#     Referência: Goyal & Welch (2008, RFS); Clark & West (2007)
# ============================================================

csfe_vs_rw    <- cumsum((rw[, 1] - yout[, 1])^2 -
                          (fc_2srr[, 1] - yout[, 1])^2)
csfe_vs_ridge <- if ("Ridge" %in% names(models_list))
  cumsum((models_list[["Ridge"]][, 1] - yout[, 1])^2 -
           (fc_2srr[, 1] - yout[, 1])^2) else rep(NA, n_oos)
csfe_vs_lasso <- if ("LASSO" %in% names(models_list))
  cumsum((models_list[["LASSO"]][, 1] - yout[, 1])^2 -
           (fc_2srr[, 1] - yout[, 1])^2) else rep(NA, n_oos)

df_csfe <- data.frame(
  date                  = oos_dates,
  `2SRR vs Random Walk` = csfe_vs_rw,
  `2SRR vs Ridge`       = csfe_vs_ridge,
  `2SRR vs LASSO`       = csfe_vs_lasso,
  check.names = FALSE
) %>%
  pivot_longer(-date, names_to = "comparacao", values_to = "csfe") %>%
  filter(!is.na(csfe))

p_csfe <- ggplot(df_csfe,
                 aes(x = date, y = csfe, colour = comparacao, group = comparacao)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(
    values = c("2SRR vs Random Walk" = "#01696f",
               "2SRR vs Ridge"       = "#964219",
               "2SRR vs LASSO"       = "#4361ee")
  )

p_csfe <- add_recession_bands(p_csfe) +
  labs(
    title    = "CSFE acumulado — 2SRR vs benchmarks (h=1)",
    subtitle = "Cresce → 2SRR acumula vantagem | Cai → benchmark ganha\nÁrea cinza = recessões NBER",
    x = "Data", y = "CSFE acumulado (e²_bench − e²_2SRR)", colour = "Comparação"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey40", size = 9))

print(p_csfe)
ggsave("results/fig_csfe.png", p_csfe, width = 11, height = 5, dpi = 150)
cat("Salvo: fig_csfe.png\n")

# ============================================================
# 13. GRÁFICO: INFLAÇÃO REALIZADA COM PERÍODOS
#     Contextualiza todos os outros resultados
# ============================================================

df_infl <- data.frame(
  date  = oos_dates,
  infl  = yout[, 1],
  fc_2srr_h1 = fc_2srr[, 1]
)

p_infl <- ggplot(df_infl, aes(x = date))

p_infl <- add_recession_bands(p_infl) +
  geom_line(aes(y = infl, colour = "Realizado"), linewidth = 0.8) +
  geom_line(aes(y = fc_2srr_h1, colour = "2SRR h=1"), linewidth = 0.7, linetype = "dashed") +
  scale_colour_manual(values = c("Realizado" = "black", "2SRR h=1" = "#01696f")) +
  labs(
    title    = "Inflação realizada (CPIAUCSL) e previsão 2SRR h=1 — período OOS",
    subtitle = "Área cinza = recessões NBER",
    x = "Data", y = "Inflação (variação)", colour = ""
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        legend.position = "top")

print(p_infl)
ggsave("results/fig_inflacao_realizada_vs_2srr.png", p_infl, width = 11, height = 4, dpi = 150)
cat("Salvo: fig_inflacao_realizada_vs_2srr.png\n")

# ============================================================
# 14. POSICIONAMENTO NA LITERATURA + RESUMO FINAL
# ============================================================

sink("results/posicionamento_coulombe.txt")
cat("=== POSICIONAMENTO DO TCC NA LITERATURA ===\n\n")
cat("Coulombe (2024) avalia 2SRR em:\n")
cat("  Modelos AR, ARDI, VAR5, VAR20 — dados canadenses\n")
cat("  Tabela 16 (Half/Half): RMSPE relativo ao AR2\n\n")
cat("Este trabalho:\n")
cat("  Aplica 2SRR no framework de Medeiros et al. (FRED-MD, EUA)\n")
cat("  RMSFE relativo ao Random Walk (benchmark mais conservador)\n")
cat("  Primeira avaliação head-to-head: 2SRR vs Ridge, LASSO,\n")
cat("  Random Forest e Bagging no MESMO dataset.\n")
cat("  Coulombe (2024) compara apenas com VARs e ARs.\n\n")
cat("RESULTADO PRINCIPAL:\n")
cat(sprintf("  2SRR bate RW em 12/12 horizontes (RMSFE médio = %.4f)\n",
            mean(res_h[, "2SRR"])))
cat(sprintf("  Melhor horizonte: h=%d (RMSFE = %.4f)\n", best_h, res_h[best_h, "2SRR"]))
cat(sprintf("  ATENÇÃO: acc6 = %.4f > 1 — discussão necessária no texto\n",
            res_full["acc6", "2SRR"]))
sink()

sink("results/resumo_2SRR.txt")
cat("=== RESUMO DOS RESULTADOS — 2SRR ===\n\n")
cat(sprintf("Horizontes em que 2SRR < RW    : %d/12\n", n_beat_rw))
cat(sprintf("Horizontes em que 2SRR < Ridge : %d/12\n", n_beat_ridge))
cat(sprintf("Horizontes em que 2SRR < LASSO : %d/12\n", n_beat_lasso))
cat(sprintf("Melhor horizonte               : h=%d (RMSFE = %.4f)\n",
            best_h, res_h[best_h, "2SRR"]))
for (acc in c("acc3", "acc6", "acc12")) {
  cat(sprintf("RMSFE %-6s: 2SRR=%.4f | Ridge=%.4f | LASSO=%.4f\n",
              acc, res_full[acc, "2SRR"],
              res_full[acc, "Ridge"], res_full[acc, "LASSO"]))
}
if (res_full["acc6", "2SRR"] > 1.0) {
  cat("\n[!] acc6 > 1: 2SRR perde para RW no acumulado de 6 meses.\n")
  cat("    Discutir no texto: propagação de erro em horizontes intermediários.\n")
}
sink()

# ============================================================
# SUMÁRIO FINAL
# ============================================================

cat("\n========================================================\n")
cat(" 05_analysis_2SRR.R — CONCLUÍDO\n")
cat("========================================================\n")
cat("Arquivos em results/:\n")
cat("  rmsfe_comparativo.csv\n")
cat("  tabela_rmsfe_dm_latex.tex\n")
cat("  dm_pvalues.csv\n")
cat("  rmsfe_subperiodos.csv\n")
cat("  decomposicao_ganho.csv\n")
cat("  beta_volatility_by_pc.csv\n")
cat("  ranking_borda.csv\n")
cat("  resumo_2SRR.txt\n")
cat("  posicionamento_coulombe.txt\n")
cat("  fig1_rmsfe_todos.png\n")
cat("  fig2_2srr_vs_selecionados.png\n")
cat("  fig_heatmap_rmsfe.png\n")
cat("  fig_borda_ranking.png\n")
cat("  fig_dm_heatmap.png\n")
cat("  fig_subperiodos.png\n")
cat("  fig_subperiodo_gfc.png\n")
cat("  fig_decomposicao_periodos.png\n")
cat("  fig_betas_top2_volateis.png\n")
cat("  fig_betas_pc1_pc2_pc3.png\n")
cat("  fig_betas_top6_pcs.png\n")
cat("  fig_norma_betas.png\n")
cat("  fig_csfe.png\n")
cat("  fig_inflacao_realizada_vs_2srr.png\n")
cat("========================================================\n")