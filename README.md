# Previsão de Inflação com Parâmetros Variando no Tempo (TVP)

Implementação do modelo **Two-Step Ridge Regression (2SRR)** de Coulombe (2025) aplicado à base de dados de Medeiros (2024) para previsão de inflação nos EUA. O projeto avalia o ganho preditivo de parâmetros variantes no tempo frente a modelos estáticos (Ridge, LASSO, Random Forest etc.) em horizontes de 1 a 12 meses, com validação econométrica completa fora-da-amostra (OOS).

---

## Referências Principais

- **Coulombe, P. G. (2025)**. Time-Varying Parameters as Ridge Regressions. *International Journal of Forecasting*, 41(3).
  - DOI: [10.1016/j.ijforecast.2024.08.006](https://doi.org/10.1016/j.ijforecast.2024.08.006)
  - Replicação: [github.com/hugocout/Replication-codes-for-Time-Varying-Parameters-as-Ridge-Regressions](https://github.com/hugocout/Replication-codes-for-Time-Varying-Parameters-as-Ridge-Regressions)

- **Medeiros, M. C. et al.** Base de dados de preditores macroeconômicos.
  - Código original: [github.com/gabrielrvsc/ForecastingInflation](https://github.com/gabrielrvsc/ForecastingInflation)

- **Nathalia, O. (2024)**. Previsão de Inflação no Brasil com Machine Learning.
  - Tese: [lume.ufrgs.br/handle/10183/272999](https://lume.ufrgs.br/handle/10183/272999)
  - Código: [github.com/nathaliaoreda/thesis_UFRGS](https://github.com/nathaliaoreda/thesis_UFRGS)

---

## Fluxo de Execução

```
01 → 02 → 03 → 03_felipe → 04_felipe → 05 → 06 → 07 → 08/09
```

Execute os scripts na ordem abaixo a partir de `forecasting_inflation/`.

---

## Scripts

### Preparação de Dados

| Script | Descrição | Output principal |
|--------|-----------|-----------------|
| `01_get_fred_data.R` | Baixa/atualiza dados do FRED via API | `data/2026-03-MD.csv` |
| `02_random_walk_oos_y.R` | Gera benchmark Random Walk e variável OOS | `forecasts/rw.rda`, `forecasts/yout.rda` |

### Modelos Baseline (Medeiros)

| Script | Descrição | Output principal |
|--------|-----------|-----------------|
| `03_call_model.R` | Estima Ridge, LASSO, Elastic Net, Random Forest, AdaLASSO, AdaElNET | `forecasts/*.rda` por modelo |

### Modelo Principal (2SRR)

| Script | Descrição | Output principal |
|--------|-----------|-----------------|
| `03_call_model_felipe.R` | Versão inicial do 2SRR (legada) | `forecasts/2SRR.rda`, `forecasts/betas_2SRR.rda` |
| `06_coulombe_2SRR_pipeline.R` | **Pipeline principal v8.0** — baixa funções Coulombe, executa loop POOS com 312 janelas, estima Ridge e 2SRR para h=1,3,6,12, salva betas para todos os horizontes | `forecasts/coulombe_fc_h{01/03/06/12}.csv`, `forecasts/coulombe_betas_2SRR.rda`, `forecasts/coulombe_betas_ridge.rda` |

>  O `06` demanda ~30 a 40h. Requer `fGarch` instalado (`install.packages("fGarch")`).

### Avaliação e Análise

| Script | Descrição | Outputs principais |
|--------|-----------|-------------------|
| `04_eval_results.R` | Avaliação básica dos modelos Medeiros (RMSFE relativo ao RW) | Tabelas e figuras em `results/` |
| `04_eval_results_felipe.R` | Avaliação completa: RMSFE de todos os modelos, gráfico comparativo, betas TVP do 2SRR (h=1) e comparação direta 2SRR vs Ridge | `results/rmsfe_comparativo.csv`, `results/fig_rmsfe_comparativo.png`, `results/fig_betas_2SRR_tempo.png`, `results/fig_2srr_vs_ridge.png` |
| `05_v2_analysis_2SRR.R` | Análise aprofundada dos parâmetros variantes: trajetórias de betas, volatilidade por componente, ranking de instabilidade | `results/beta_traj_h**.csv`, `results/beta_var_h**.csv` |
| `07_validacao_econometrica.R` | **Validação completa v3.0** — 13 partes: métricas (RMSE, MAE, MAPE), Diebold-Mariano pairwise, Clark-West, Mincer-Zarnowitz, Forecast Encompassing, Giacomini-White, CSSED + Fluctuation Test, Model Confidence Set, análise TVP dos betas, sub-amostras (pré/pós COVID), gráficos e tabelas LaTeX | Ver seção *Arquivos de Resultado* |
| `08_tvp_comparison.R` | Compara desempenho preditivo entre três especificações TVP (AR, Factor, FAVAR) | `results/tvp_comparison_table.csv` |
| `09_tvp_analysis.R` | Análise adicional das especificações TVP: tabela RMSE por caso, DM pairwise entre 2SRR de cada especificação, gráficos comparativos | `results/figures/tvp_comparison_rmse.pdf`, `results/figures/tvp_comparison_ratio.pdf` |

---

## Especificações TVP Comparadas (`08` e `09`)

| Caso | Especificação | Descrição |
|------|--------------|-----------|
| 1 | **TVP-AR** | Somente lags de y (univariado) |
| 2 | **TVP-Factor** | Somente fatores PCA (sem lags) |
| 3 | **TVP-FAVAR** | Fatores PCA + lags de y — modelo completo (rodado em `06`) |

---

## Arquivos de Resultado (`results/`)

| Arquivo | O que é |
|---------|---------|
| `tabela_completa.csv` | RMSE, MAE, rank e ratio vs 2SRR de todos os modelos por horizonte |
| `dm_results.csv` | Testes Diebold-Mariano (2SRR vs Ridge): estatística DM e p-valor por horizonte |
| `tabela_forecast.tex` / `tabela_dm.tex` | Tabelas prontas para LaTeX |
| `cssed_h{01/03/06/12}.csv` | Diferença acumulada de erro quadrático (2SRR vs Ridge) ao longo do tempo — curva caindo = 2SRR ganha |
| `beta_traj_h{01/03/06/12}.csv` | Trajetória temporal dos coeficientes do 2SRR por janela OOS — evolução dos parâmetros ao longo do tempo |
| `beta_var_h{01/03/06/12}.csv` | Variância e coeficiente de variação (CV) por beta — resume quais componentes mais variam em cada horizonte |
| `betas_2srr_h{01/03/06/12}.csv` | Coeficientes estimados pelo 2SRR em cada janela OOS (escala original) |
| `betas_ridge_h{01/03/06/12}.csv` | Coeficientes estimados pelo Ridge estático em cada janela — contraste direto com o 2SRR |
| `tvp_comparison_table.csv` | RMSE e ratio do 2SRR vs Ridge para cada especificação TVP (AR, Factor, FAVAR) e horizonte |
| `figures/` | Gráficos de forecast, CSSED, betas TVP e comparação de especificações (`.pdf`) |

---

## Estrutura de Dados

```
data/
├── 2026-03-MD.csv          ← dados brutos do FRED (atualizado)
├── data.rda                ← dados compilados em formato R
└── medeiros/               ← base de Medeiros

forecasts/
├── coulombe_fc_h{01/03/06/12}.csv   ← previsões Ridge e 2SRR por horizonte
├── coulombe_betas_2SRR.rda          ← betas 2SRR por janela (todos horizontes)
├── coulombe_betas_ridge.rda         ← betas Ridge por janela (todos horizontes)
├── tvp_TVP_AR_h{**}.csv             ← previsões TVP-AR
├── tvp_TVP_Factor_h{**}.csv         ← previsões TVP-Factor
├── rw.rda                           ← benchmark Random Walk
├── yout.rda                         ← valores realizados OOS
└── *.rda                            ← modelos Medeiros (LASSO, RF, etc.)
```

---

## Observações

- Validação cruzada **k-fold** (k=5) em janelas móveis OOS com **312 janelas**.
- Hiperparâmetros do 2SRR: K_max = 40 fatores PCA, variância explicada ≥ 90%, lambda via CV.
- Checkpoints automáticos a cada 50 janelas em `checkpoints/` para retomada de execução.
- Pacotes necessários: `glmnet`, `pracma`, `fGarch`, `ggplot2`, `forecast`, `lmtest`, `sandwich`, `car`, `xtable`, `MCS`.
