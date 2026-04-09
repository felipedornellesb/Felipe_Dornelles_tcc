# ============================================================
# 02_FORECAST.R
#
# Roda o modelo TVP-2SRR (Coulombe 2022) em janela expandida
# para todas as combinações variable × horizon definidas no
# grid all_options carregado do 00_data_download.R.
#
# Inputs  (de 10_data/data_MM_DD_YYYY/):
#   df_transf.rda   → painel transformado e estacionário
#   targets_br.rda  → lista V1..V5 com nomes das variáveis
#   all_options.rda → grid M × V × H (60 combinações)
#
# Outputs (em 30_output/):
#   coulombe_2srr.rds   → lista com previsões por var × h
#   actual_values.rds   → matriz de valores realizados
#   results_summary.rds → MSFE, RMSFE, MAE, DM vs RW
# ============================================================

rm(list = ls())

# ============================================================
# SETUP
# ============================================================

wd <- "C:/Users/felip/OneDrive/Documentos/tcc/Felipe_Dornelles_tcc/"
setwd(wd)

paths <- list(
  program   = "00_program",
  data      = "10_data",
  tools     = "20_tools",
  functions = "20_tools/functions",
  output    = "30_output",
  results   = "40_results"
)

# Pasta datada dentro de 30_output/ — ex.: 30_output/outputs_04_08_2026/
run_date    <- format(Sys.Date(), "%m_%d_%Y")
output_run  <- file.path(paths$output, paste0("outputs_", run_date))

for (p in c(paths$output, output_run, paths$results)) {
  if (!dir.exists(p)) dir.create(p, recursive = TRUE)
}

cat(sprintf("Output folder: %s\n", output_run))

# ============================================================
# PACKAGES
# ============================================================

myPKGs <- c("dplyr", "randomForest", "mboost", "e1071", "readr",
            "pracma", "glmnet", "fGarch", "matrixcalc", "forecast")

InstalledPKGs    <- names(installed.packages()[, "Package"])
InstallThesePKGs <- myPKGs[!myPKGs %in% InstalledPKGs]
if (length(InstallThesePKGs) > 0)
  install.packages(InstallThesePKGs, repos = "http://cran.us.r-project.org")
invisible(lapply(myPKGs, library, character.only = TRUE))

# ============================================================
# LOAD TOOLS & FUNCTIONS
# ============================================================

source(file.path(paths$functions, "00_Nathalia_functions.R"))
source(file.path(paths$functions, "MV2SRR_v221103.R"))

# ============================================================
# CARREGA DADOS DO RUN MAIS RECENTE EM 10_data/
# ------------------------------------------------------------
# Sempre usa a pasta data_MM_DD_YYYY mais recente disponível,
# evitando hardcodar a data. Se quiser fixar uma data
# específica, substitua por:
#   run_folder <- "10_data/data_04_08_2026"
# ============================================================

data_subfolders <- list.dirs(paths$data, recursive = FALSE, full.names = TRUE)
data_subfolders <- data_subfolders[grepl("data_\\d{2}_\\d{2}_\\d{4}$",
                                         basename(data_subfolders))]

if (length(data_subfolders) == 0)
  stop("Nenhuma pasta data_MM_DD_YYYY encontrada em 10_data/. Rode o 00_data_download.R primeiro.")

run_folder <- data_subfolders[length(data_subfolders)]   # pega a mais recente
cat(sprintf("Loading data from: %s\n", run_folder))

load(file.path(run_folder, "df_transf.rda"))    # df_transf
load(file.path(run_folder, "targets_br.rda"))   # targets_br (V1..V5)
load(file.path(run_folder, "all_options.rda"))  # all_options (grid M x V x H)

# ============================================================
# DEFINE VARIÁVEIS-ALVO E HORIZONTE
# ------------------------------------------------------------
# targets_br vem do 00_data_download.R com os nomes limpos
# (PIB, IPCA, SELIC, CAMBIO, DESEMPREGO). Usa esses nomes
# para indexar df_transf — sem depender dos códigos IPEA.
# ============================================================

variable_list <- unlist(targets_br)   # V1="PIB", V2="IPCA", ...
horizon_list  <- c(1, 3, 6, 12)
lag_orders    <- 1:6
model_name    <- "coulombe_2srr"
window_type   <- "expanding"          # "expanding" ou "rolling"

# Confirma que todas as variáveis-alvo estão no painel
ausentes <- variable_list[!variable_list %in% colnames(df_transf)]
if (length(ausentes) > 0)
  stop(paste("Variaveis ausentes em df_transf:", paste(ausentes, collapse = ", ")))

# Remove coluna de data antes de escalar
df_num <- df_transf |> select(-date)

# Normaliza o painel inteiro (z-score) — necessário para o TVP-ridge
df_scaled <- as.data.frame(scale(df_num))

# ============================================================
# JANELA DE PREVISÃO
# ------------------------------------------------------------
# Período in-sample encerra em Dez/2008.
# df_transf começa em Fev/1996 (Jan/1996 perdido na diff).
# Dez/2008 = 155 meses desde Fev/1996 → start_window = 156
# (o índice 156 é a primeira previsão OOS, Jan/2009).
# Ajuste start_window se seu df_transf começar em outra data.
# ============================================================

n_total      <- nrow(df_scaled)
start_window <- 156        # primeiro índice OOS = Jan/2009
end_window   <- n_total
nwindows     <- end_window - start_window + 1

cat(sprintf("Panel: %d obs x %d vars | OOS windows: %d (index %d to %d)\n",
            n_total, ncol(df_scaled), nwindows, start_window, end_window))

stopifnot("start_window deve ser >= 50" = start_window >= 50)
stopifnot("start_window deve ser < n_total" = start_window <= n_total)

actual_values <- as.matrix(df_scaled[start_window:end_window, variable_list,
                                     drop = FALSE])

# ============================================================
# dataprep_generic
# ------------------------------------------------------------
# Prepara Xin, Xout, yin para qualquer janela de df já
# transformado. Usa os K primeiros fatores principais de X
# como regressores, concatenados com lags de y — exatamente
# como dataprep(..., dataset="B0ARDI") da Nathalia, mas sem
# depender do nome do dataset.
# ============================================================

dataprep_generic <- function(df, horizon, variable, lag = 12, K = 6) {

  df    <- as.data.frame(df)
  n     <- nrow(df)
  y_all <- df[[variable]]

  # Regressores: todas as colunas exceto a variável-alvo
  x_mat <- as.matrix(df[, setdiff(colnames(df), variable), drop = FALSE])

  # 1. Remove colunas com NA ou Inf
  cols_ok <- apply(x_mat, 2, function(col) all(is.finite(col)))
  x_mat   <- x_mat[, cols_ok, drop = FALSE]

  # 2. Remove colunas com variância zero (constantes dentro da janela)
  cols_var <- apply(x_mat, 2, function(col) var(col) > 0)
  x_mat    <- x_mat[, cols_var, drop = FALSE]

  # 3. Garante que y também está limpo
  if (any(!is.finite(y_all)))
    stop(sprintf("dataprep_generic: y contém NA/Inf (var=%s)", variable))

  if (ncol(x_mat) == 0L)
    stop(sprintf("dataprep_generic: nenhuma coluna válida em x (n=%d, h=%d)",
                 n, horizon))

  K_actual <- min(K, ncol(x_mat), n - 1L)

  # Fatores principais
  factors  <- prcomp(x_mat, scale. = FALSE)$x[, 1:K_actual, drop = FALSE]

  # X = [y_lags | fatores_lags]
  x_full  <- cbind(y_all, factors)
  X_embed <- embed(as.matrix(x_full), lag)

  T_embed <- nrow(X_embed)
  y_start <- lag + horizon
  y_end   <- n
  n_obs   <- y_end - y_start + 1L

  if (n_obs <= 0L)
    stop(sprintf("dataprep_generic: window too small (n=%d, lag=%d, h=%d)",
                 n, lag, horizon))

  Xin  <- X_embed[1:n_obs,  , drop = FALSE]
  yin  <- y_all[y_start:y_end]
  Xout <- X_embed[T_embed,   , drop = FALSE]

  stopifnot(nrow(Xin) == length(yin))

  list(Xin = Xin, Xout = as.vector(Xout), yin = yin)
}

# ============================================================
# func_coulombe_2srr
# ------------------------------------------------------------
# Wrapper para tvp.ridge() (MV2SRR_v221103.R).
# Re-otimiza lambda via CV quando reoptimize_hyperparameters
# = TRUE; nas demais janelas reutiliza best$lambda.
# ============================================================

func_coulombe_2srr <- function(df, horizon, variable, lag_orders,
                               reoptimize_hyperparameters = FALSE,
                               best = list()) {

  lag  <- if (!is.null(best$lag)) best$lag else max(lag_orders)
  prep <- dataprep_generic(df, horizon = horizon, variable = variable,
                           lag = lag, K = 6L)

  # Grid de lambdas do paper (Coulombe 2022)
  lambdavec <- exp(pracma::linspace(-6, 20, n = 15))

  # Na primeira janela nunca há lambda salvo — força re-otimização
  if (is.null(best$lambda)) reoptimize_hyperparameters <- TRUE

  result <- tryCatch({

    lambda_use <- if (reoptimize_hyperparameters) lambdavec else best$lambda

    out <- tvp.ridge(
      X                 = prep$Xin,
      Y                 = matrix(prep$yin, ncol = 1),
      lambda.candidates = lambda_use,
      oosX              = prep$Xout,
      kfold             = 5L,
      CV.2SRR           = TRUE,
      CV.plot           = FALSE,
      sig.eps.param     = 0.75,
      sig.u.param       = 0.75
    )

    # Atualiza best apenas quando re-otimizou
    if (reoptimize_hyperparameters)
      best <- list(lag = lag, lambda = out$lambdas)

    list(pred = as.numeric(out$forecast), best = best)

  }, error = function(e) {
    message(sprintf("  [ERROR tvp.ridge] var=%s h=%d: %s",
                    variable, horizon, e$message))
    list(pred = NA_real_, best = best)
  })

  result
}

# ============================================================
# FORECASTING LOOP
# ------------------------------------------------------------
# Loop externo: variáveis × horizontes
# Loop interno: janelas OOS (expanding window)
# Lambda re-otimizado a cada 12 janelas e sempre na primeira.
# ============================================================

set.seed(1234)
forecasts_list <- list()

for (v in variable_list) {

  for (h in horizon_list) {

    cat(sprintf("\n[%s] var=%-12s | h=%2d\n", model_name, v, h))

    best       <- list()
    model_list <- vector("list", nwindows)

    for (i in seq_len(nwindows)) {

      # Índice da última observação de treino nesta janela
      train_end <- start_window - h - 1L + i

      if (train_end < 30L) {
        model_list[[i]] <- list(pred = NA_real_)
        next
      }

      Df <- if (window_type == "expanding") {
        df_scaled[1:train_end, , drop = FALSE]
      } else {
        df_scaled[max(1L, train_end - 119L):train_end, , drop = FALSE]
      }

      # Re-otimiza lambda na primeira janela e a cada 12
      reopt <- (i == 1L) || (i %% 12L == 1L)
      if (reopt) cat(sprintf("  re-optimizing lambda (window %d/%d)\n",
                             i, nwindows))

      model <- func_coulombe_2srr(
        df                         = Df,
        horizon                    = h,
        variable                   = v,
        lag_orders                 = lag_orders,
        reoptimize_hyperparameters = reopt,
        best                       = best
      )

      # Fallback: replica previsão anterior se der NA
      if (is.na(model$pred)) {
        model$pred <- if (i > 1L && !is.na(model_list[[i - 1L]]$pred))
          model_list[[i - 1L]]$pred
        else 0
        cat(sprintf("  [fallback] window %d uses previous pred\n", i))
      }

      model_list[[i]] <- list(pred = model$pred)
      best            <- model$best
    }

    preds <- vapply(model_list, function(x) x$pred, numeric(1L))

    forecasts_list[[length(forecasts_list) + 1L]] <- list(
      variable = v,
      horizon  = h,
      model    = model_name,
      pred     = preds
    )

    cat(sprintf("  -> %d forecasts | NAs: %d\n",
                length(preds), sum(is.na(preds))))
  }
}

# ============================================================
# SALVA PREVISÕES
# ============================================================

saveRDS(forecasts_list,
        file = file.path(paths$output, paste0(model_name, ".rds")))
saveRDS(actual_values,
        file = file.path(paths$output, "actual_values.rds"))

cat(sprintf("\nForecast completo. %d combinacoes (var x horizonte) salvas em %s/\n",
            length(forecasts_list), paths$output))

# ============================================================
# CHECKUPS — MSFE, RMSFE, MAE
# ------------------------------------------------------------
# err[i] = pred[i] - actual[i]: alinha diretamente pois
# pred[i] já é a previsão h passos à frente do passo i.
# ============================================================

forecasts <- readRDS(file.path(paths$output, "coulombe_2srr.rds"))
actual    <- readRDS(file.path(paths$output, "actual_values.rds"))

results <- do.call(rbind, lapply(forecasts, function(f) {
  act <- actual[, f$variable]
  n   <- min(length(f$pred), length(act))
  err <- f$pred[1:n] - act[1:n]
  data.frame(
    variable = f$variable,
    horizon  = f$horizon,
    MSFE     = mean(err^2,    na.rm = TRUE),
    RMSFE    = sqrt(mean(err^2, na.rm = TRUE)),
    MAE      = mean(abs(err), na.rm = TRUE)
  )
}))

cat("\n=== In-sample fit stats ===\n")
print(results, digits = 4)

# ============================================================
# BENCHMARK — RANDOM WALK
# ------------------------------------------------------------
# RW: pred[t+h] = act[t]  (sem deriva)
# MSFE_rel < 1 = 2SRR bate o RW
# ============================================================

rw_results <- do.call(rbind, lapply(forecasts, function(f) {
  act <- actual[, f$variable]
  h   <- f$horizon
  n   <- length(act)

  if (n <= h) return(NULL)

  rw_pred <- act[1:(n - h)]
  rw_act  <- act[(1 + h):n]
  rw_err  <- rw_pred - rw_act

  data.frame(
    variable = f$variable,
    horizon  = h,
    MSFE_RW  = mean(rw_err^2, na.rm = TRUE)
  )
}))

# MSFE relativo ao RW
results_rel <- merge(results, rw_results, by = c("variable", "horizon"))
results_rel$MSFE_rel <- results_rel$MSFE / results_rel$MSFE_RW

cat("\n=== MSFE relativo ao Random Walk (< 1 = 2SRR superior) ===\n")
print(results_rel[, c("variable", "horizon", "MSFE", "MSFE_RW",
                      "MSFE_rel", "RMSFE")],
      digits = 4, row.names = FALSE)

# ============================================================
# TESTE DIEBOLD-MARIANO — 2SRR vs Random Walk
# ------------------------------------------------------------
# H0: MSFE iguais | H1 (less): MSFE do 2SRR < MSFE do RW
# p < 0.05 → 2SRR significativamente superior ao RW
# ============================================================

dm_results <- do.call(rbind, lapply(forecasts, function(f) {
  act <- actual[, f$variable]
  h   <- f$horizon
  n   <- min(length(f$pred), length(act))

  e1 <- f$pred[1:n] - act[1:n]                    # erro 2SRR

  # Erro RW: alinhado com o mesmo índice de e1
  # e2[i] = act[i] - act[i - h]  (pred RW h passos antes)
  e2 <- rep(NA_real_, n)
  for (i in (h + 1L):n) e2[i] <- act[i - h] - act[i]

  # Remove NAs para o teste
  keep <- !is.na(e1) & !is.na(e2)
  if (sum(keep) < 10L){
    return(data.frame(variable = f$variable,
                      horizon  = h,
                      DM_stat  = NA_real_,
                      p_value  = NA_real_,
                      MSFE_rel = NA_real_))
  }

  dm_test <- tryCatch(
    forecast::dm.test(e1[keep], e2[keep],
                      alternative = "less",
                      h           = h,
                      power       = 2),
    error = function(e) NULL
  )

  p_val <- if (!is.null(dm_test)) dm_test$p.value else NA_real_
  stat  <- if (!is.null(dm_test)) dm_test$statistic else NA_real_

  data.frame(
    variable = f$variable,
    horizon  = h,
    DM_stat  = as.numeric(stat),
    p_value  = p_val,
    MSFE_rel = mean(e1[keep]^2, na.rm = TRUE) /
               mean(e2[keep]^2, na.rm = TRUE)
  )
}))

cat("\n=== Teste Diebold-Mariano: 2SRR vs Random Walk ===\n")
cat("H0: MSFE iguais | H1 (less): MSFE do 2SRR < MSFE do RW\n\n")
print(dm_results[, c("variable", "horizon", "DM_stat",
                      "p_value", "MSFE_rel")],
      digits = 4, row.names = FALSE)

# ============================================================
# CONSOLIDA E SALVA RESULTADOS FINAIS
# ============================================================

results_summary <- list(
  fit_stats  = results,
  msfe_rel   = results_rel,
  dm_vs_rw   = dm_results
)

saveRDS(results_summary,
        file = file.path(paths$results, "results_summary.rds"))

# Exporta também em CSV para leitura rápida
write.csv(results,
          file      = file.path(paths$results, "TAB01_fit_stats.csv"),
          row.names = FALSE)

write.csv(results_rel,
          file      = file.path(paths$results, "TAB02_msfe_relativo.csv"),
          row.names = FALSE)

write.csv(dm_results,
          file      = file.path(paths$results, "TAB03_diebold_mariano.csv"),
          row.names = FALSE)

cat(sprintf("\n✅ Resultados salvos em %s/\n", paths$results))
cat("   TAB01_fit_stats.csv\n")
cat("   TAB02_msfe_relativo.csv\n")
cat("   TAB03_diebold_mariano.csv\n")
cat("   results_summary.rds\n")

# ============================================================
# RESUMO FINAL NO CONSOLE
# ============================================================

cat("\n=== RESUMO GERAL ===\n")
cat(sprintf("Modelo      : %s\n",   model_name))
cat(sprintf("Janela      : %s\n",   window_type))
cat(sprintf("Variáveis   : %s\n",   paste(variable_list, collapse = ", ")))
cat(sprintf("Horizontes  : %s\n",   paste(horizon_list,  collapse = ", ")))
cat(sprintf("Combinações : %d\n",   length(forecasts_list)))
cat(sprintf("OOS windows : %d\n",   nwindows))
cat(sprintf("Período OOS : índices %d a %d\n", start_window, end_window))

cat("\n--- MSFE relativo ao RW (< 1 = 2SRR superior) ---\n")
tab_print <- results_rel[, c("variable", "horizon", "MSFE_rel", "RMSFE")]
tab_print  <- tab_print[order(tab_print$variable, tab_print$horizon), ]
print(tab_print, digits = 3, row.names = FALSE)

cat("\n--- DM p-valores (< 0.05 = 2SRR sign. superior ao RW) ---\n")
dm_print <- dm_results[, c("variable", "horizon", "p_value", "DM_stat")]
dm_print  <- dm_print[order(dm_print$variable, dm_print$horizon), ]
print(dm_print, digits = 3, row.names = FALSE)

cat("\n 02_FORECAST.R finalizado.\n")