# Previsão de Inflação com Parâmetros Variando no Tempo

## Descrição do Projeto

Implementação de modelos de previsão de inflação com ênfase na metodologia de **Two-Step Ridge Regression (2SRR)** com parâmetros variantes no tempo. O projeto adapta a implementação original de Coulombe (2025) para dados do FRED e da base de Medeiros (2024), permitindo análise de estabilidade estrutural em relações inflacionárias.

**Nota**: Inclui tentativas exploratórias de integração com dados brasileiros (Nathalia, 2024), sem resultados conclusivos.

---

## Referências Principais

- **Coulombe, P. G. (2025)**. Time-Varying Parameters as Ridge Regressions. *International Journal of Forecasting*, 41(3).
  - DOI: [10.1016/j.ijforecast.2024.08.006](https://doi.org/10.1016/j.ijforecast.2024.08.006)
  - Replicação: [github.com/hugocout/Replication-codes-for-Time-Varying-Parameters-as-Ridge-Regressions](https://github.com/hugocout/Replication-codes-for-Time-Varying-Parameters-as-Ridge-Regressions)

- **Medeiros, M. C. et al. (2024)???**. Não encontrei de onde o professor encontrou tal base, a definir a referência correta.
  - Base de dados: [github.com/gabrielrvsc/ForecastingInflation](https://github.com/gabrielrvsc/ForecastingInflation)

- **Vasconcellos, G. (2024)**. Não encontrei de onde o professor encontrou tal base, a definir a referência correta.
  - Código: [github.com/gabrielrvsc/ForecastingInflation](https://github.com/gabrielrvsc/ForecastingInflation)

- **Nathalia, O. (2024)**. Previsão de Inflação no Brasil com Machine Learning.
  - Tese: [lume.ufrgs.br/handle/10183/272999](https://lume.ufrgs.br/handle/10183/272999)
  - Código: [github.com/nathaliaoreda/thesis_UFRGS](https://github.com/nathaliaoreda/thesis_UFRGS)

---

## Estrutura dos Scripts

| Script | Descrição |
|--------|-----------|
| `01_get_fred_data.R` | Baixa/atualiza dados do FRED (executar conforme necessário) |
| `02_random_walk_oos_y.R` | Gera benchmarks: `rw.rda` (random walk) e `yout.rda` (dados OOS) |
| `03_call_model.R` | Executa modelos baseline (Ridge, Random Forest, LASSO, Elastic Net, etc.) |
| `03_call_model_felipe.R` | **Principal**: executa 2SRR, salva `2SRR.rda` e `betas_2SRR.rda` |
| `04_eval_results_felipe.R` | Gera tabelas e gráficos de desempenho comparativo |
| `05_v2_analysis_2SRR.R` | Análise complementar de parâmetros variantes no tempo |

---

## Estrutura de Dados

```
data/
├── 2026-03-MD.csv          ← dados brutos do FRED (atualizado)
├── data.rda                ← dados compilados em formato R
├── fred/                   ← série histórica do FRED (se houver)
├── medeiros/               ← base de Medeiros (2024)
└── [outras subpastas]      ← dados auxiliares e intermediários
```

---

## Estrutura de Resultados

```
results/
├── Arquivos de resumo
│   ├── resumo_2SRR.txt
│   ├── posicionamento_coulombe.txt
│   └── ranking_borda.csv
├── Tabelas comparativas
│   ├── rmsfe_comparativo.csv
│   ├── rmsfe_subperiodos.csv
│   ├── dm_pvalues.csv              ← testes Diebold-Mariano
│   ├── contrib_media_pcs.csv       ← contribuição por componente
│   └── tabela_rmsfe_dm_latex.tex   ← tabela em formato LaTeX
├── Betas (parâmetros variantes)
│   ├── betas_2SRR.rda              ← série completa de coeficientes
│   ├── beta_volatility_by_pc.csv   ← volatilidade por horizonte
│   ├── betas_2SRR_h1.csv a h12.csv ← betas por horizonte individual
│   └── [subpastas com configs diferentes: k, lvar, kfold]
├── Previsões
│   ├── 2SRR.rda
│   ├── rw.rda              ← benchmark random walk
│   ├── yout.rda            ← valores observados OOS
│   ├── AdaElNET.rda, AdaLASSO.rda, ...
│   └── [arquivos .rda para cada modelo]
├── Documentação
│   ├── relatorio_2SRR.tex
│   ├── relatorio_2SRR.toc
│   ├── relatorio_2SRR.aux
│   ├── relatorio_2SRR.fdb_latexmk
│   └── [arquivos LaTeX auxiliares]
└── resultsk*/ (subpastas com configs)
    ├── resultsk20lvar085kfold3/
    ├── resultskpca40lvar090kfold5/
    └── [outras combinações de hiperparâmetros]
```

---

## Fluxo de Execução Recomendado

1. **Preparação de dados**: `01_get_fred_data.R` (se necessário)
2. **Benchmarks**: `02_random_walk_oos_y.R`
3. **Modelos baseline**: `03_call_model.R`
4. **Modelo principal**: `03_call_model_felipe.R` (2SRR)
5. **Análise e visualização**: `04_eval_results_felipe.R`
6. **Análise aprofundada**: `05_v2_analysis_2SRR.R`

---

## Observações

- Os modelos utilizam validação cruzada **k-fold** em janelas móveis fora-da-amostra (OOS).
- Hiperparâmetros do 2SRR: $K_{\max} = 40$, variância explicada ≥ 90%, $k = 5$ dobras.
- Tempos de execução variam: ~8-10 horas para o 2SRR com 312 janelas OOS.
- Resultados organizados em subpastas para diferentes combinações de hiperparâmetros.