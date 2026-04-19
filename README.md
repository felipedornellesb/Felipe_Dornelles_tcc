Tentando modificar o código de 2SRR do coulombe, para utilizar com a base do Medeiros. Contem também tentativa de analisar dados de Brasil com os da Nathalia, porem sem resultados satisfatórios.

https://github.com/gabrielrvsc/ForecastingInflation

https://github.com/nathaliaoreda/thesis_UFRGS
https://github.com/hugocout/Replication-codes-for-Time-Varying-Parameters-as-Ridge-Regressions/

https://doi.org/10.1016/j.ijforecast.2024.08.006 (Coulombe - 2025)
https://lume.ufrgs.br/handle/10183/272999# (Nathalia - 2024)

## Organização dos scripts:
01_get_fred_data.R          ← só se quiser atualizar os dados
02_random_walk_oos_y.R      ← gera rw.rda e yout.rda
03_call_model.R             ← roda os modelos do Medeiros (Ridge, RF, etc.)
03_call_model_felipe.R      ← roda o 2SRR, salva 2SRR.rda e betas_2SRR.rda
04_eval_results_felipe.R    ← gera tabelas e gráficos em results/

## Organização dos dados:
data/                      ← pasta com os dados brutos
data/fred/                 ← dados do FRED
data/medeiros/              ← dados do Medeiros (2024)
data/2SRR/                 ← dados intermediários do 2SRR (resíduos, etc.)
data/results/              ← resultados finais (previsões, betas, etc.)

## Organização dos resultados:
results/                    ← pasta com os resultados finais
results/plots/              ← gráficos de desempenho dos modelos
results/tables/             ← tabelas de desempenho dos modelos
results/2SRR/              ← resultados específicos do 2SRR (betas, etc.)