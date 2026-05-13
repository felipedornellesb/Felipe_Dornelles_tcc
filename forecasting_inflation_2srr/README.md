# Inflation Forecasting with 2SRR

Two-Step Ridge Regression (Coulombe, IJF 2025) applied to US CPI inflation forecasting, benchmarked against Ridge, LASSO, Random Forest and other models from the Medeiros framework.

## How to Run

```r
# Set working directory to project root, then:
source("00_prog/00_setup.R")            # install packages, download functions
source("00_prog/01_data_prep.R")        # generate yout.rda and rw.rda
source("00_prog/02_forecast_medeiros.R") # run Medeiros baseline models
source("00_prog/03_forecast_2srr.R")    # run 2SRR
source("00_prog/04_analysis.R")         # RMSFE, tests, betas, figures
```

Place `data.rda` in `10_data/` before running. This file comes from the Medeiros repository or from the advisor.

## Structure

```
00_prog/               Scripts (run in order)
10_data/               data.rda (FRED-MD)
20_tools/
  21_coulombe/         Coulombe functions (10 files, auto-downloaded)
  22_medeiros/         Medeiros functions (auto-downloaded)
  23_adapted/          tvp_functions.R (standalone 2SRR + run2srr wrapper)
30_output/             forecasts/, betas/, checkpoints/
40_results/            tables/, figures/
```

## Coulombe Functions (21_coulombe/)

| File | Key function | Role |
|------|-------------|------|
| dualGRRmdA_v190215.R | `dualGRR()` | Dual ridge solver (Eq. 9, 11) |
| TVPRR_v181111.R | `TVPRR()` | TVP ridge core |
| TVPRRcosso_v181120.R | `TVPRR_cosso()` | Algorithm 1 orchestrator |
| CVGSBHK_v181127.R | `cvgs.bhk2015()` | Cross-validation |
| CVKFMV_v190214.R | (CV helpers) | Additional CV routines |
| zfun_v190304.R | `Zfun()`, `make_reg_matrix()` | Z matrix construction |
| fastZrot_v181125.R | `fastZrot()` | Fast Z rotation |
| EM_sw.R | `EM_sw()` | Stock-Watson EM imputation |
| ICp2.R | `ICp2()` | Information criteria |
| Xgenerators_v190127.R | (data generators) | Simulation support |

The `factor()` PCA function is defined in `00_setup.R` before sourcing, because `EM_sw()` calls `factor(X, n_fac=n)` which conflicts with `base::factor()`.

## Required Packages

`glmnet`, `pracma`, `randomForest`, `forecast`, `lmtest`, `sandwich`, `ggplot2`, `reshape2`, `xtable`, `HDeconometrics` (from GitHub), `rugarch` or `fGarch` (for GARCH in Step 2).

## References

- Coulombe, P. G. (2025). Time-Varying Parameters as Ridge Regressions. IJF 41(3), 982-1002.
- Medeiros et al. github.com/gabrielrvsc/ForecastingInflation

Felipe Dornelles, UFRGS 2026
