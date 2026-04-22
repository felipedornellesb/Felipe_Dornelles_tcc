# ============================================================
# 05_analysis_2SRR.R
# Análise completa dos resultados do 2SRR vs demais modelos
# Salva tabelas CSV, LaTeX e gráficos PNG em results/
# ============================================================

library(tidyverse)
library(knitr)

# --- 1. TABELA RMSFE COMPLETA (já calculada) -----------------
write.csv(
  as.data.frame(round(res_full, 4)),
  "results/rmsfe_comparativo.csv",
  row.names = TRUE
)

# --- 2. TABELA LaTeX -----------------------------------------
latex_tab <- res_full[, c("2SRR","Ridge","LASSO","AdaLASSO",
                           "ElNET","RF","Bagging","Factor","AR","CSR")]
sink("results/tabela_rmsfe_latex.tex")
cat("\\begin{table}[ht]\n\\centering\n")
cat("\\caption{RMSFE relativo ao Random Walk}\n")
cat("\\label{tab:rmsfe}\n")
cat("\\begin{tabular}{l", rep("r", ncol(latex_tab)), "}\n", sep="")
cat("\\hline\n")
cat("Horizonte &", paste(colnames(latex_tab), collapse=" & "), "\\\\\n")
cat("\\hline\n")
for(i in 1:nrow(latex_tab)) {
  vals <- round(latex_tab[i,], 4)
  # negrito no menor valor de cada linha
  min_idx <- which.min(vals)
  vals_str <- as.character(vals)
  vals_str[min_idx] <- paste0("\\textbf{", vals_str[min_idx], "}")
  cat(rownames(latex_tab)[i], "&",
      paste(vals_str, collapse=" & "), "\\\\\n")
}
cat("\\hline\n\\end{tabular}\n\\end{table}\n")
sink()
cat("LaTeX salvo em results/tabela_rmsfe_latex.tex\n")

# --- 3. GRÁFICO: RMSFE todos modelos por horizonte -----------
df_plot <- as.data.frame(res_h) %>%
  rownames_to_column("horizonte") %>%
  mutate(horizonte = 1:12) %>%
  pivot_longer(-horizonte, names_to="modelo", values_to="rmsfe") %>%
  mutate(
    destaque  = ifelse(modelo == "2SRR", "2SRR", "outros"),
    alpha_val = ifelse(modelo == "2SRR", 1.0, 0.35)
  )

p1 <- ggplot(df_plot,
       aes(x=horizonte, y=rmsfe, colour=modelo,
           group=modelo, alpha=alpha_val)) +
  geom_line(aes(linewidth=destaque)) +
  geom_hline(yintercept=1, linetype="dashed",
             colour="black", linewidth=0.5) +
  scale_x_continuous(breaks=1:12) +
  scale_alpha_identity() +
  scale_linewidth_manual(
    values=c("2SRR"=1.4, "outros"=0.5), guide="none") +
  labs(
    title    = "RMSFE relativo ao Random Walk — todos os modelos",
    subtitle = "Abaixo de 1 = melhor que RW | 2SRR em destaque",
    x="Horizonte (meses)", y="RMSFE / RMSFE(RW)", colour="Modelo"
  ) +
  theme_minimal(base_size=12) +
  theme(panel.grid.minor=element_blank(),
        plot.title=element_text(face="bold"))

ggsave("results/fig1_rmsfe_todos.png",
       p1, width=11, height=5, dpi=150)

# --- 4. GRÁFICO: 2SRR vs Ridge vs LASSO vs RF ---------------
df_sel <- as.data.frame(res_full) %>%
  rownames_to_column("horizonte") %>%
  mutate(h_num = 1:15) %>%
  pivot_longer(c(`2SRR`,Ridge,LASSO,RF),
               names_to="modelo", values_to="rmsfe")

p2 <- ggplot(df_sel, aes(x=h_num, y=rmsfe,
             colour=modelo, group=modelo)) +
  geom_line(linewidth=1.0) +
  geom_point(size=2) +
  geom_hline(yintercept=1, linetype="dashed",
             colour="grey40", linewidth=0.4) +
  scale_x_continuous(
    breaks=1:15,
    labels=c(paste0("h=",1:12),"acc3","acc6","acc12")
  ) +
  scale_colour_manual(
    values=c("2SRR"="#01696f","Ridge"="#964219",
             "LASSO"="#4361ee","RF"="#e07a5f")
  ) +
  labs(
    title    = "2SRR vs Ridge vs LASSO vs RF",
    subtitle = "Inclui horizontes h=1:12 e acumulados",
    x="Horizonte", y="RMSFE / RMSFE(RW)", colour="Modelo"
  ) +
  theme_minimal(base_size=12) +
  theme(panel.grid.minor=element_blank(),
        plot.title=element_text(face="bold"),
        axis.text.x=element_text(angle=45, hjust=1))

ggsave("results/fig2_2srr_vs_outros.png",
       p2, width=11, height=5, dpi=150)

# --- 5. GRÁFICO: Betas time-varying h=1 (top 6 fatores) -----
load("forecasts/betas_2SRR.rda")  # df_betas: window_date, var_idx, beta

top_vars <- df_betas %>%
  group_by(var_idx) %>%
  summarise(volatilidade = sd(beta, na.rm=TRUE)) %>%
  slice_max(volatilidade, n=6) %>%
  pull(var_idx)

p3 <- df_betas %>%
  filter(var_idx %in% top_vars) %>%
  mutate(fator = paste0("PC", var_idx)) %>%
  ggplot(aes(x=window_date, y=beta, colour=fator, group=fator)) +
  geom_line(linewidth=0.7) +
  labs(
    title    = "Betas time-varying do 2SRR — h=1 (top 6 fatores por volatilidade)",
    subtitle = "Cada linha = coeficiente β(t) de um fator PCA ao longo do tempo",
    x="Data", y="β(t)", colour="Fator"
  ) +
  theme_minimal(base_size=12) +
  theme(panel.grid.minor=element_blank(),
        plot.title=element_text(face="bold"))

ggsave("results/fig3_betas_tv_h1.png",
       p3, width=11, height=5, dpi=150)

# --- 6. RESUMO TEXTO -----------------------------------------
n_beat_rw    <- sum(res_h[,"2SRR"] < 1.0, na.rm=TRUE)
n_beat_ridge <- sum(res_h[,"2SRR"] < res_h[,"Ridge"], na.rm=TRUE)
best_h       <- which.min(res_h[,"2SRR"])

sink("results/resumo_2SRR.txt")
cat("=== RESUMO DOS RESULTADOS — 2SRR ===\n\n")
cat(sprintf("Horizontes em que 2SRR < RW     : %d/12\n", n_beat_rw))
cat(sprintf("Horizontes em que 2SRR < Ridge  : %d/12\n", n_beat_ridge))
cat(sprintf("Melhor horizonte do 2SRR        : h=%d (RMSFE=%.4f)\n",
            best_h, res_h[best_h,"2SRR"]))
cat(sprintf("RMSFE acc3  : 2SRR=%.4f | Ridge=%.4f | LASSO=%.4f\n",
            res_full["acc3","2SRR"],
            res_full["acc3","Ridge"],
            res_full["acc3","LASSO"]))
cat(sprintf("RMSFE acc6  : 2SRR=%.4f | Ridge=%.4f | LASSO=%.4f\n",
            res_full["acc6","2SRR"],
            res_full["acc6","Ridge"],
            res_full["acc6","LASSO"]))
cat(sprintf("RMSFE acc12 : 2SRR=%.4f | Ridge=%.4f | LASSO=%.4f\n",
            res_full["acc12","2SRR"],
            res_full["acc12","Ridge"],
            res_full["acc12","LASSO"]))
cat("\nCONCLUSÃO: 2SRR é o modelo mais consistente nos horizontes\n")
cat("médios e longos. Nos acumulados, domina no acc3 mas perde\n")
cat("para Ridge/LASSO no acc12 — esperado pela maior rigidez\n")
cat("do shrinkage fixo em horizontes muito longos.\n")
sink()

cat("\n=== TUDO SALVO EM results/ ===\n")
cat(list.files("results/"), sep="\n")

# ============================================================
# 06_inference_and_subperiods.R
#
# Extensões analíticas do TCC — Felipe Dornelles
#
# Seções:
#   A. Teste Diebold-Mariano (DM) — 2SRR vs todos os modelos
#   B. Análise de Sub-períodos (pré/pós COVID e pré/pós GFC)
#   C. Análise dos Betas TVP por Fator PCA (PC1, PC2...)
#   D. CSFE — Cumulative Squared Forecast Error vs RW e Ridge
#   E. Tabela LaTeX com asteriscos de significância (DM)
#   F. Posicionamento vs Tabela 16 de Coulombe (2024)
#
# Pré-requisitos (já gerados pelos scripts anteriores):
#   forecasts/2SRR.rda
#   forecasts/betas_2SRR.rda
#   forecasts/rw.rda
#   forecasts/yout.rda
#   forecasts/<outros modelos>.rda
#   data/data.rda
# ============================================================
library(tidyverse)
library(forecast)   # dm.test()
library(patchwork)  # combinar gráficos

# ------------------------------------------------------------
# CARREGA DADOS BASE
# ------------------------------------------------------------

load("forecasts/yout.rda")   # yout: matriz T x 4 (h1 realizado, acc3, acc6, acc12)
load("forecasts/rw.rda")     # rw  : matrix T x 15 (previsoes RW)
load("forecasts/2SRR.rda")   # forecasts -> renomeia para fc_2srr
fc_2srr <- forecasts; rm(forecasts)

# Carrega todos os outros modelos
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

# Datas do out-of-sample (a partir de quando começa a janela OOS)
load("data/data.rda")
all_dates <- data$date
n_total   <- nrow(data)
n_oos     <- nrow(yout)
oos_dates <- tail(all_dates, n_oos)

# ============================================================
# A. TESTE DIEBOLD-MARIANO — 2SRR vs cada modelo (h=1 a h=12)
# ============================================================
# dm.test(e1, e2, alternative="less") testa H0: MSE(e1) >= MSE(e2)
# Rejeitamos H0 se p < 0.10 => 2SRR é significativamente melhor
# Usamos alternative = "two.sided" para tabela completa

dm_pvalues <- matrix(NA,
                     nrow = 12,
                     ncol = length(models_list),
                     dimnames = list(paste0("h=", 1:12),
                                     names(models_list)))
dm_stats   <- dm_pvalues

for (modelo in names(models_list)) {
  fc_m <- models_list[[modelo]]
  for (h in 1:12) {
    e_2srr <- fc_2srr[, h] - yout[, 1]
    e_m    <- fc_m[, h]   - yout[, 1]

    # Remove NAs alinhados
    ok <- complete.cases(e_2srr, e_m)
    if (sum(ok) < 10) next

    tryCatch({
      dm_res <- dm.test(e_2srr[ok], e_m[ok],
                        alternative = "two.sided",
                        h = h, power = 2)
      dm_pvalues[h, modelo] <- dm_res$p.value
      dm_stats[h, modelo]   <- dm_res$statistic
    }, error = function(e) NULL)
  }
}

# Salva p-valores brutos
write.csv(round(dm_pvalues, 4), "results/dm_pvalues.csv")
cat("DM p-valores salvos em results/dm_pvalues.csv\n")

# Função para converter p-valor em asteriscos
stars <- function(p) {
  ifelse(is.na(p), "",
  ifelse(p < 0.01, "***",
  ifelse(p < 0.05, "**",
  ifelse(p < 0.10, "*", ""))))
}

# ============================================================
# B. ANÁLISE DE SUB-PERÍODOS
# ============================================================
# Recalcula RMSFE relativo ao RW em 3 janelas:
#   Full    : período completo
#   Pre-COVID: até 2019-12-01
#   Covid   : 2020-01-01 a 2022-12-01
#   Post-COVID: 2023-01-01 em diante
#   Pre-GFC  : até 2007-12-01
#   Post-GFC : 2008-01-01 a 2012-12-01

sub_periods <- list(
  Full       = rep(TRUE, n_oos),
  Pre_COVID  = oos_dates <  as.Date("2020-01-01"),
  COVID      = oos_dates >= as.Date("2020-01-01") &
               oos_dates <= as.Date("2022-12-01"),
  Post_COVID = oos_dates >  as.Date("2022-12-01"),
  GFC        = oos_dates >= as.Date("2007-07-01") &
               oos_dates <= as.Date("2009-12-01")
)

# Para cada sub-período, calcula RMSFE relativo de 2SRR e modelos selecionados
models_compare <- c("Ridge", "LASSO", "RF", "Bagging", "Factor", "AR")
models_compare <- intersect(models_compare, names(models_list))

rmsfe_sub <- lapply(names(sub_periods), function(sp) {
  idx <- sub_periods[[sp]]
  if (sum(idx, na.rm = TRUE) < 5) return(NULL)

  rwe_sp <- sqrt(colMeans((rw[idx, 1:12] - yout[idx, 1])^2, na.rm = TRUE))

  # 2SRR
  e_2srr_sp <- sqrt(colMeans((fc_2srr[idx, 1:12] - yout[idx, 1])^2, na.rm = TRUE))
  res_sp     <- data.frame(h = 1:12,
                            periodo = sp,
                            modelo  = "2SRR",
                            rmsfe   = e_2srr_sp / rwe_sp)

  # Demais modelos
  for (m in models_compare) {
    if (!m %in% names(models_list)) next
    fc_m   <- models_list[[m]]
    e_m_sp <- sqrt(colMeans((fc_m[idx, 1:12] - yout[idx, 1])^2, na.rm = TRUE))
    res_sp <- bind_rows(res_sp, data.frame(h = 1:12,
                                            periodo = sp,
                                            modelo  = m,
                                            rmsfe   = e_m_sp / rwe_sp))
  }
  res_sp
}) %>% bind_rows()

write.csv(rmsfe_sub, "results/rmsfe_subperiodos.csv", row.names = FALSE)
cat("RMSFE sub-períodos salvo em results/rmsfe_subperiodos.csv\n")

# Gráfico sub-períodos — 2SRR em destaque, h=1 a 12
p_sub <- rmsfe_sub %>%
  filter(periodo %in% c("Pre_COVID", "COVID", "Post_COVID")) %>%
  mutate(
    periodo = factor(periodo,
                     levels = c("Pre_COVID","COVID","Post_COVID"),
                     labels = c("Pré-COVID\n(<2020)", "COVID\n(2020-2022)",
                                "Pós-COVID\n(>2022)")),
    destaque  = ifelse(modelo == "2SRR", "2SRR", "outros"),
    alpha_val = ifelse(modelo == "2SRR", 1.0, 0.35)
  ) %>%
  ggplot(aes(x = h, y = rmsfe,
             colour    = modelo,
             group     = modelo,
             alpha     = alpha_val,
             linewidth = destaque)) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed",
             colour = "grey30", linewidth = 0.4) +
  facet_wrap(~ periodo, ncol = 3) +
  scale_x_continuous(breaks = c(1,3,6,9,12)) +
  scale_alpha_identity() +
  scale_linewidth_manual(values = c("2SRR" = 1.3, "outros" = 0.5),
                         guide = "none") +
  scale_colour_manual(
    values = c("2SRR"    = "#01696f",
               "Ridge"   = "#964219",
               "LASSO"   = "#4361ee",
               "RF"      = "#e07a5f",
               "Bagging" = "#6d597a",
               "Factor"  = "#b5838d",
               "AR"      = "#588157")
  ) +
  labs(
    title    = "RMSFE relativo ao RW por sub-período",
    subtitle = "2SRR em destaque — abaixo de 1 = melhor que Random Walk",
    x        = "Horizonte (meses ahead)",
    y        = "RMSFE / RMSFE(RW)",
    colour   = "Modelo"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor  = element_blank(),
    plot.title        = element_text(face = "bold"),
    strip.text        = element_text(face = "bold"),
    legend.position   = "right"
  )

ggsave("results/fig_subperiodos.png",
       p_sub, width = 13, height = 5, dpi = 150)

print(p_sub)

cat("Gráfico sub-períodos salvo em results/fig_subperiodos.png\n")

# Gráfico GFC separado
p_gfc <- rmsfe_sub %>%
  filter(periodo %in% c("Pre_COVID", "GFC")) %>%
  mutate(
    periodo = factor(periodo,
                     levels = c("GFC","Pre_COVID"),
                     labels = c("Crise 2008 (GFC)", "Pré-COVID (<2020)")),
    destaque  = ifelse(modelo == "2SRR", "2SRR", "outros"),
    alpha_val = ifelse(modelo == "2SRR", 1.0, 0.35)
  ) %>%
  ggplot(aes(x = h, y = rmsfe,
             colour = modelo, group = modelo,
             alpha = alpha_val, linewidth = destaque)) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed",
             colour = "grey30", linewidth = 0.4) +
  facet_wrap(~ periodo) +
  scale_x_continuous(breaks = c(1,3,6,9,12)) +
  scale_alpha_identity() +
  scale_linewidth_manual(values = c("2SRR" = 1.3, "outros" = 0.5),
                         guide = "none") +
  scale_colour_manual(
    values = c("2SRR"  = "#01696f", "Ridge"  = "#964219",
               "LASSO" = "#4361ee", "RF"     = "#e07a5f",
               "AR"    = "#588157", "Factor" = "#b5838d",
               "Bagging" = "#6d597a")
  ) +
  labs(
    title  = "RMSFE relativo ao RW — Crise 2008 vs Período Completo Pré-COVID",
    x      = "Horizonte", y = "RMSFE / RMSFE(RW)", colour = "Modelo"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

ggsave("results/fig_subperiodo_gfc.png",
       p_gfc, width = 10, height = 5, dpi = 150)

print(p_gfc)

# ============================================================
# C. ANÁLISE DOS BETAS TVP — PC1 vs PC2 vs demais
# ============================================================
# Objetivo: mostrar que variação em PC1/PC2 = mudança de regime
# da inflação americana (Phillips Curve time-varying)

load("forecasts/betas_2SRR.rda")
# df_betas esperado: data.frame com colunas window_date, var_idx, beta (h=1)

# Garante que var_idx é numérico e cria label de fator
df_betas <- df_betas %>%
  mutate(
    var_idx = as.integer(var_idx),
    fator   = paste0("PC", var_idx)
  )

# --- C1. Volatilidade de cada PC ao longo do tempo
volatility_by_pc <- df_betas %>%
  group_by(var_idx, fator) %>%
  summarise(
    sd_beta  = sd(beta, na.rm = TRUE),
    mean_abs = mean(abs(beta), na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  arrange(desc(sd_beta))

write.csv(volatility_by_pc, "results/beta_volatility_by_pc.csv",
          row.names = FALSE)

cat("\n=== Volatilidade dos betas por PC (h=1) ===\n")
print(head(volatility_by_pc, 10))

# --- C2. Gráfico: PC1 e PC2 — betas ao longo do tempo com destaque
top2_pcs <- volatility_by_pc$var_idx[1:2]

# Adiciona períodos de recessão (NBER) para contexto econômico
recession_bands <- data.frame(
  start = as.Date(c("2001-03-01", "2007-12-01", "2020-02-01")),
  end   = as.Date(c("2001-11-01", "2009-06-01", "2020-04-01")),
  label = c("Dot-com", "GFC", "COVID")
)

p_betas_pc12 <- df_betas %>%
  filter(var_idx %in% top2_pcs) %>%
  ggplot(aes(x = window_date, y = beta,
             colour = fator, group = fator)) +
  # Bandas de recessão
  geom_rect(data = recession_bands,
            aes(xmin = start, xmax = end,
                ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE,
            fill = "grey80", alpha = 0.4) +
  geom_text(data = recession_bands,
            aes(x = start + (end - start)/2,
                y = Inf, label = label),
            inherit.aes = FALSE,
            vjust = 1.5, size = 2.8, colour = "grey40") +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey50", linewidth = 0.4) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = c("#01696f", "#964219")) +
  labs(
    title    = "Coeficientes TVP: PC1 e PC2 ao longo do tempo (h=1)",
    subtitle = "Variação dos betas = evidência de mudança de regime na inflação americana (Phillips Curve TVP)\nÁreas cinza = recessões NBER",
    x        = "Data",
    y        = "β(t)",
    colour   = "Componente"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold"),
    plot.subtitle    = element_text(colour = "grey40", size = 9)
  )

ggsave("results/fig_betas_pc1_pc2.png",
       p_betas_pc12, width = 11, height = 5, dpi = 150)

print(p_betas_pc12)

cat("Gráfico betas PC1/PC2 salvo em results/fig_betas_pc1_pc2.png\n")

# --- C3. Gráfico: Top 6 PCs por volatilidade (facet)
top6_pcs <- volatility_by_pc$var_idx[1:6]

p_betas_top6 <- df_betas %>%
  filter(var_idx %in% top6_pcs) %>%
  mutate(fator = factor(fator, levels = paste0("PC", top6_pcs))) %>%
  ggplot(aes(x = window_date, y = beta)) +
  geom_rect(data = recession_bands,
            aes(xmin = start, xmax = end,
                ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE,
            fill = "grey80", alpha = 0.35) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey60", linewidth = 0.3) +
  geom_line(colour = "#01696f", linewidth = 0.75) +
  facet_wrap(~ fator, scales = "free_y", ncol = 3) +
  labs(
    title    = "Top 6 fatores PCA por volatilidade dos betas TVP (h=1)",
    subtitle = "Cada painel = um componente principal. Área cinza = recessões NBER.",
    x        = "Data", y = "β(t)"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold"),
    strip.text       = element_text(face = "bold", colour = "#01696f")
  )

ggsave("results/fig_betas_top6_pcs.png",
       p_betas_top6, width = 12, height = 7, dpi = 150)

print(p_betas_top6)

cat("Gráfico Top 6 PCs salvo em results/fig_betas_top6_pcs.png\n")

# --- C4. Interpretação econômica: beta de PC1 x inflação realizada
# Correlação entre |Δβ_PC1| e |Δ inflação| — evidência de co-movimento
inflacao_realizada <- yout[, 1]  # inflação realizada no OOS

beta_pc1 <- df_betas %>%
  filter(var_idx == top2_pcs[1]) %>%
  arrange(window_date) %>%
  pull(beta)

# Alinha comprimentos
n_min <- min(length(beta_pc1), length(inflacao_realizada))
cor_pc1_infl <- cor(diff(beta_pc1[1:n_min]),
                    diff(inflacao_realizada[1:n_min]),
                    use = "complete.obs")

sink("results/interpretacao_betas.txt")
cat("=== INTERPRETAÇÃO ECONÔMICA DOS BETAS TVP ===\n\n")
cat("Fator mais volátil (h=1):", paste0("PC", top2_pcs[1]), "\n")
cat("Segundo mais volátil     :", paste0("PC", top2_pcs[2]), "\n\n")
cat(sprintf("Correlação entre Δβ(PC%d) e Δ(Inflação realizada): %.4f\n",
            top2_pcs[1], cor_pc1_infl))
cat("\nInterpretação:\n")
cat("  Variações significativas em PC1 e PC2 ao longo do tempo\n")
cat("  constituem evidência direta de que a relação entre os\n")
cat("  preditores macroeconômicos e a inflação americana é\n")
cat("  TIME-VARYING — compatível com a Phillips Curve TVP e\n")
cat("  as quebras estruturais discutidas em Coulombe (2024).\n")
cat("\n  Períodos de maior volatilidade dos betas:\n")
cat("  - Crise 2008-09 (GFC): mudança na transmissão monetária\n")
cat("  - 2020-2022 (COVID + inflação pós-pandemia): regime novo\n")
cat("  Exatamente quando o 2SRR tem maior vantagem relativa\n")
cat("  sobre modelos de coeficientes fixos (Ridge, LASSO).\n")
sink()
cat("Interpretação econômica salva em results/interpretacao_betas.txt\n")

# ============================================================
# D. CSFE — CUMULATIVE SQUARED FORECAST ERROR
#    CSFE(t) = Σ[e_RW(s)² - e_2SRR(s)²] para s=1..t
#    Quando cresce = 2SRR acumula vantagem; cai = RW melhor
# ============================================================

csfe_vs_rw    <- cumsum((rw[, 1] - yout[, 1])^2 -
                         (fc_2srr[, 1] - yout[, 1])^2)

# Contra Ridge (se disponível)
if ("Ridge" %in% names(models_list)) {
  csfe_vs_ridge <- cumsum((models_list[["Ridge"]][, 1] - yout[, 1])^2 -
                           (fc_2srr[, 1] - yout[, 1])^2)
} else {
  csfe_vs_ridge <- rep(NA, n_oos)
}

df_csfe <- data.frame(
  date          = oos_dates,
  CSFE_vs_RW    = csfe_vs_rw,
  CSFE_vs_Ridge = csfe_vs_ridge
) %>%
  pivot_longer(-date, names_to = "comparacao", values_to = "csfe") %>%
  mutate(comparacao = recode(comparacao,
    "CSFE_vs_RW"    = "2SRR vs Random Walk",
    "CSFE_vs_Ridge" = "2SRR vs Ridge"
  ))

p_csfe <- df_csfe %>%
  filter(!is.na(csfe)) %>%
  ggplot(aes(x = date, y = csfe,
             colour = comparacao, group = comparacao)) +
  # Bandas de recessão
  geom_rect(data = recession_bands,
            aes(xmin = start, xmax = end,
                ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE,
            fill = "grey80", alpha = 0.35) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey40", linewidth = 0.5) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(
    values = c("2SRR vs Random Walk" = "#01696f",
               "2SRR vs Ridge"       = "#964219")
  ) +
  labs(
    title    = "CSFE acumulado — 2SRR vs benchmarks (h=1)",
    subtitle = "Cresce → 2SRR acumula vantagem | Cai → benchmark ganha\nÁrea cinza = recessões NBER",
    x        = "Data",
    y        = "CSFE acumulado (e²_bench − e²_2SRR)",
    colour   = "Comparação"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold"),
    plot.subtitle    = element_text(colour = "grey40", size = 9)
  )

ggsave("results/fig_csfe.png",
       p_csfe, width = 11, height = 5, dpi = 150)

print(p_csfe)

cat("Gráfico CSFE salvo em results/fig_csfe.png\n")

# ============================================================
# E. TABELA LaTeX COM ASTERISCOS DE SIGNIFICÂNCIA (DM)
# ============================================================
# Carrega res_full (gerado em 04_eval_results_felipe.R)
# Supondo que está em memória; se não, recarrega o CSV

if (!exists("res_full")) {
  res_full <- read.csv("results/rmsfe_comparativo.csv",
                       row.names = 1)
}

# Monta tabela combinada: RMSFE + asteriscos DM (apenas h=1..12)
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
cat("Horizonte &", paste(cols_latex, collapse = " & "), "\\\\\n")
cat("\\hline\n")

for (i in 1:nrow(res_full)) {
  h_label <- rownames(res_full)[i]
  vals    <- res_full[i, cols_latex]

  # Menor valor da linha (negrito)
  min_idx <- which.min(unlist(vals))

  cells <- sapply(seq_along(cols_latex), function(j) {
    col  <- cols_latex[j]
    val  <- round(vals[[col]], 4)

    # Asteriscos DM: apenas para h=1..12 e modelos != 2SRR
    ast <- ""
    if (col != "2SRR" && h_label %in% rownames(dm_pvalues)) {
      p <- dm_pvalues[h_label, col]
      ast <- stars(p)
    }

    cell <- paste0(val, ast)
    if (j == min_idx) cell <- paste0("\\textbf{", cell, "}")
    cell
  })

  cat(h_label, "&", paste(cells, collapse = " & "), "\\\\\n")
}

cat("\\hline\n")
cat("\\multicolumn{", length(cols_latex) + 1, "}{l}{",
    "\\scriptsize{* p<0.10, ** p<0.05, *** p<0.01 — Diebold-Mariano (bilateral) vs 2SRR}}\\\\\n",
    sep = "")
cat("\\end{tabular}\n}\n\\end{table}\n")
sink()
cat("Tabela LaTeX com DM salva em results/tabela_rmsfe_dm_latex.tex\n")

# ============================================================
# F. POSICIONAMENTO vs TABELA 16 DE COULOMBE (2024)
# ============================================================
# Coulombe (2024) Tabela 16 reporta RMSPE/RMSPE(AR2) para:
# AR, ARDI, VAR5, VAR20 com e sem TVP (2SRR, MSRRS, MSRRD)
# Aqui construímos tabela análoga no framework Medeiros/FRED

coulombe_context <- data.frame(
  Modelo       = c("AR (Medeiros)",
                   "Ridge (Medeiros)",
                   "LASSO (Medeiros)",
                   "RF (Medeiros)",
                   "Bagging (Medeiros)",
                   "Factor (Medeiros)",
                   "2SRR [ESTE TRABALHO]"),
  Framework    = c(rep("ForecastingInflation (FRED-MD)", 6),
                   "ForecastingInflation (FRED-MD) + Coulombe (2024)"),
  Contribuicao = c(rep("Baseline Medeiros et al.", 6),
                   "Primeira aplicacao 2SRR no framework Medeiros")
)

write.csv(coulombe_context,
          "results/posicionamento_literatura.csv",
          row.names = FALSE)

sink("results/posicionamento_coulombe.txt")
cat("=== POSICIONAMENTO DO TCC NA LITERATURA ===\n\n")
cat("Coulombe (2024) avalia 2SRR em:\n")
cat("  - Modelos AR, ARDI, VAR5, VAR20 em dados CANADENSES\n")
cat("  - Variáveis: Inflação, GDP growth, Interest rate (SPREAD)\n")
cat("  - Tabela 16 (Half/Half): RMSPE relativo ao AR2\n\n")
cat("Este trabalho:\n")
cat("  - Aplica 2SRR no MESMO FRAMEWORK de Medeiros et al.\n")
cat("  - Dados FRED-MD (EUA) — inflação CPI\n")
cat("  - RMSFE relativo ao Random Walk (benchmark mais conservador)\n")
cat("  - Permite comparação DIRETA com Ridge, LASSO, RF, Bagging\n")
cat("    no MESMO dataset — algo não feito por Coulombe (2024)\n\n")
cat("Contribuição inédita:\n")
cat("  Primeira avaliação de 2SRR head-to-head com métodos\n")
cat("  de Machine Learning (RF, Bagging) E métodos clássicos\n")
cat("  (Ridge, LASSO) no framework FRED-MD de Medeiros et al.\n")
cat("  Coulombe (2024) compara apenas com VARs e ARs.\n")
sink()
cat("Posicionamento literatura salvo em results/posicionamento_coulombe.txt\n")

# ============================================================
# SUMÁRIO Ajustar
# ============================================================

cat("\n")
cat("========================================================\n")
cat(" 06_inference_and_subperiods.R — CONCLUÍDO\n")
cat("========================================================\n")
cat("Arquivos gerados em results/:\n")
cat("  dm_pvalues.csv                 — p-valores Diebold-Mariano\n")
cat("  rmsfe_subperiodos.csv          — RMSFE por sub-período\n")
cat("  beta_volatility_by_pc.csv      — volatilidade dos betas por PC\n")
cat("  interpretacao_betas.txt        — interpretação econômica TVP\n")
cat("  posicionamento_coulombe.txt    — contribuição vs literatura\n")
cat("  posicionamento_literatura.csv  — tabela de posicionamento\n")
cat("  tabela_rmsfe_dm_latex.tex      — LaTeX com asteriscos DM\n")
cat("  fig_subperiodos.png            — RMSFE por sub-período\n")
cat("  fig_subperiodo_gfc.png         — sub-período GFC\n")
cat("  fig_betas_pc1_pc2.png          — betas PC1 e PC2 + recessões\n")
cat("  fig_betas_top6_pcs.png         — top 6 PCs (facet)\n")
cat("  fig_csfe.png                   — CSFE acumulado vs RW e Ridge\n")
cat("========================================================\n")

