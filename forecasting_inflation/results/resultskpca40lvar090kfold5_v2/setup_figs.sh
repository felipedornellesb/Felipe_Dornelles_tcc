#!/usr/bin/env bash
# setup_figs.sh
# Baixa todas as figuras do repositório GitHub para a pasta local.
# Execute UMA VEZ antes de compilar o .tex:
#   chmod +x setup_figs.sh && ./setup_figs.sh

BASE="https://raw.githubusercontent.com/felipedornellesb/Felipe_Dornelles_tcc/main/forecasting_inflation/results"

FIGS=(
  "fig1_rmsfe_todos.png"
  "fig2_2srr_vs_selecionados.png"
  "fig_2srr_vs_ridge.png"
  "fig_decomposicao_periodos.png"
  "fig_borda_ranking.png"
  "fig_csfe.png"
  "fig_betas_2SRR_tempo.png"
  "fig_betas_top2_volateis.png"
  "fig_betas_pc1_pc2_pc3.png"
  "fig_betas_top6_pcs.png"
  "fig_contrib_pcs.png"
)

for fig in "${FIGS[@]}"; do
  if [ ! -f "$fig" ]; then
    echo "Baixando $fig..."
    curl -sL "$BASE/$fig" -o "$fig"
  else
    echo "Já existe: $fig (pulando)"
  fi
done

echo "Concluído. Compile com: pdflatex relatorio_2SRR.tex"