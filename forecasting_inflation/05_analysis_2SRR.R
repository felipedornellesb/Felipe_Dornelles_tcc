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