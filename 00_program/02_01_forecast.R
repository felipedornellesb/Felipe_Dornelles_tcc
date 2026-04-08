# =============================================================================================================
# 02_FORECAST.R
# =============================================================================================================

rm(list = ls())

wd <- 'C:/Users/felip/OneDrive/Documentos/tcc/Felipe_Dornelles_tcc/'
setwd(wd)

paths <- list(
  program   = "00_program",
  data      = "10_data",
  tools     = "20_tools",
  functions = "20_tools/functions",
  output    = "30_output",
  results   = "40_results"
)

# =============================================================================================================
# PACKAGES
# =============================================================================================================

myPKGs <- c('dplyr', 'randomForest', 'mboost', 'e1071', 'readr',
            'pracma', 'glmnet', 'fGarch', 'matrixcalc')

InstalledPKGs    <- names(installed.packages()[, 'Package'])
InstallThesePKGs <- myPKGs[!myPKGs %in% InstalledPKGs]
if (length(InstallThesePKGs) > 0)
  install.packages(InstallThesePKGs, repos = "http://cran.us.r-project.org")
invisible(lapply(myPKGs, library, character.only = TRUE))

# =============================================================================================================
# LOAD TOOLS & FUNCTIONS
# =============================================================================================================

source(paste(paths$functions, '00_Nathalia_functions.R', sep = '/'))
source(paste(paths$functions,     'MV2SRR_v221103.R',        sep = '/'))

# =============================================================================================================
# dataprep_generic
# Prepara Xin, Xout, yin para qualquer janela de df já transformado.
# Usa os K primeiros fatores principais de X como regressores,
# concatenados com lags de y — exatamente como dataprep(..., dataset="B0ARDI")
# da Nathalia, mas sem depender do nome do dataset.
# =============================================================================================================

dataprep_generic <- function(df, horizon, variable, lag = 12, K = 6) {
  
  df    <- as.data.frame(df)
  n     <- nrow(df)
  y_all <- df[[variable]]
  
  # regressores: todas as colunas exceto a variável-alvo
  x_mat    <- as.matrix(df[, setdiff(colnames(df), variable), drop = FALSE])
  K_actual <- min(K, ncol(x_mat), n - 1)
  
  # fatores principais (df já está scaled no loop principal)
  factors  <- prcomp(x_mat, scale. = FALSE)$x[, 1:K_actual, drop = FALSE]
  
  # X = [y | fatores]
  x_full   <- cbind(y_all, factors)
  X_embed  <- embed(as.matrix(x_full), lag)   # (n - lag + 1) x (lag * (K+1))
  
  T_embed  <- nrow(X_embed)                   # = n - lag + 1
  # yin[t] = y em t + lag - 1 + horizon   (h passos a frente)
  y_start  <- lag + horizon                   # primeiro índice válido de y
  y_end    <- n                               # último índice de y
  n_obs    <- y_end - y_start + 1            # número de pares (X, y)
  
  if (n_obs <= 0)
    stop(sprintf("dataprep_generic: janela muito pequena (n=%d, lag=%d, h=%d)", n, lag, horizon))
  
  Xin  <- X_embed[1:n_obs,       , drop = FALSE]
  yin  <- y_all[(y_start):(y_end)]
  Xout <- X_embed[T_embed,        , drop = FALSE]   # última linha: previsão OOS
  
  stopifnot(nrow(Xin) == length(yin))
  
  return(list(Xin = Xin, Xout = as.vector(Xout), yin = yin))
}

# =============================================================================================================
# func_coulombe_2srr
# Wrapper para tvp.ridge() (MV2SRR_v221103.R).
# Re-otimiza lambda via CV quando reoptimize_hyperparameters = TRUE;
# nas demais janelas reutiliza o lambda salvo em best$lambda.
# =============================================================================================================

func_coulombe_2srr <- function(df, horizon, variable, lag_orders,
                               reoptimize_hyperparameters = FALSE,
                               type = "cv", best = list(), ...) {
  
  lag <- if (!is.null(best$lag)) best$lag else max(lag_orders)
  
  prep  <- dataprep_generic(df, horizon = horizon, variable = variable,
                            lag = lag, K = 6)
  Xin   <- prep$Xin
  yin   <- prep$yin
  Xout  <- prep$Xout
  
  # grid de lambdas do paper (Coulombe 2022, eq. CV)
  lambdavec <- exp(pracma::linspace(-6, 20, n = 15))
  
  # na primeira janela nunca há lambda salvo — força re-otimização
  if (is.null(best$lambda)) reoptimize_hyperparameters <- TRUE
  
  result <- tryCatch({
    
    lambda_use <- if (reoptimize_hyperparameters) lambdavec else best$lambda
    
    out <- tvp.ridge(
      X                 = Xin,
      Y                 = matrix(yin, ncol = 1),
      lambda.candidates = lambda_use,
      oosX              = Xout,
      kfold             = 5,
      CV.2SRR           = TRUE,
      CV.plot           = FALSE,
      sig.eps.param     = 0.75,
      sig.u.param       = 0.75
    )
    
    # atualiza best apenas quando re-otimizou
    if (reoptimize_hyperparameters) {
      best <- list(lag = lag, lambda = out$lambdas)
    }
    
    list(pred = as.numeric(out$forecast), best = best)
    
  }, error = function(e) {
    message(sprintf("  [ERRO tvp.ridge] var=%s h=%d: %s", variable, horizon, e$message))
    list(pred = NA_real_, best = best)
  })
  
  return(result)
}

# =============================================================================================================
# FORECASTING LOOP
# =============================================================================================================

load(paste(paths$data, "df.rda", sep = '/'))

# Confirma variáveis-alvo presentes
IPCA   <- "PRECOS12_IPCA12"
UNRATE <- "SEADE12_TDAGSP12"
SPREAD <- "JPM366_EMBI366"
variable_list <- c(IPCA, UNRATE, SPREAD)

ausentes <- variable_list[!variable_list %in% colnames(df)]
if (length(ausentes) > 0)
  stop(paste("Variáveis não encontradas no df.rda:", paste(ausentes, collapse = ", ")))

# Normaliza DEPOIS de confirmar colunas
df_scaled <- as.data.frame(scale(df))

# ---- Janela de previsão ----
# df.rda começa em Jan/1996 (linha 1). Após burn-in de 2 obs no 01_data_prep,
# linha 1 corresponde a Mar/1996.
# Treino encerra em Dez/2008:
#   Dec/2008 = linha (279 - (Mai/2019 - Dez/2008 em meses))
#            = 279 - 125 = 154  =>  start_window = 155
# Ajuste se o número de linhas do seu df.rda for diferente de 279.
n_total      <- nrow(df_scaled)
start_window <- n_total - 125        # primeiro passo OOS = Jan/2009
end_window   <- n_total
nwindows     <- end_window - start_window + 1

cat(sprintf("df: %d obs x %d vars | janelas OOS: %d (%d a %d)\n",
            n_total, ncol(df_scaled), nwindows, start_window, end_window))

stopifnot("start_window deve ser >= 50" = start_window >= 50)

actual_values <- as.matrix(df_scaled[start_window:end_window, variable_list])
horizon_list  <- c(1, 3, 6, 12)
lag_orders    <- 1:6
model_name    <- "coulombe_2srr"
window_type   <- "expanding"   # "expanding" ou "rolling"

set.seed(1234)
forecasts_list <- list()

for (v in variable_list) {
  
  for (h in horizon_list) {
    
    cat(sprintf("\n[%s] var=%s | h=%d\n", model_name, v, h))
    
    best       <- list()
    model_list <- vector("list", nwindows)
    
    for (i in seq_len(nwindows)) {
      
      # índice da última obs de treino desta janela
      train_end <- start_window - h - 1 + i
      
      if (train_end < 30) {
        model_list[[i]] <- list(pred = NA_real_)
        next
      }
      
      Df <- if (window_type == "expanding") {
        df_scaled[1:train_end, , drop = FALSE]
      } else {
        df_scaled[max(1, train_end - 119):train_end, , drop = FALSE]
      }
      
      # re-otimiza lambda a cada 12 janelas (e sempre na primeira)
      reopt <- (i == 1L) || (i %% 12L == 1L)
      if (reopt) cat(sprintf("  re-otimizando lambda (janela %d)\n", i))
      
      model <- func_coulombe_2srr(
        df                         = Df,
        horizon                    = h,
        variable                   = v,
        lag_orders                 = lag_orders,
        reoptimize_hyperparameters = reopt,
        type                       = "cv",
        best                       = best
      )
      
      # fallback: replica previsão anterior se der NA
      if (is.na(model$pred)) {
        model$pred <- if (i > 1L && !is.na(model_list[[i - 1L]]$pred))
          model_list[[i - 1L]]$pred
        else 0
        cat(sprintf("  [fallback] janela %d usa pred anterior\n", i))
      }
      
      model_list[[i]] <- list(pred = model$pred)
      best            <- model$best
    }
    
    preds <- vapply(model_list, function(x) x$pred, numeric(1))
    
    forecasts_list[[length(forecasts_list) + 1]] <- list(
      variable = v,
      horizon  = h,
      model    = model_name,
      pred     = preds
    )
    
    cat(sprintf("  -> %d previsoes | NAs: %d\n", length(preds), sum(is.na(preds))))
  }
}

# =============================================================================================================
# SAVE
# =============================================================================================================

saveRDS(forecasts_list,
        file = paste(paths$output, paste0(model_name, ".rds"), sep = "/"))
saveRDS(actual_values,
        file = paste(paths$output, "actual_values.rds", sep = "/"))

cat(sprintf("\nForecast completo. %d combinacoes (var x horizonte) salvas em 30_output/\n",
            length(forecasts_list)))

# =============================================================================================================
# CHECKUPS
# =============================================================================================================

# MSFE, RMSFE e MAE para cada combinação (variável x horizonte) para tabela base

forecasts <- readRDS(paste(paths$output, "coulombe_2srr.rds", sep = "/"))
actual    <- readRDS(paste(paths$output, "actual_values.rds",  sep = "/"))

results <- do.call(rbind, lapply(forecasts, function(f) {
  act  <- actual[, f$variable]
  # alinha: pred[i] é previsão h passos à frente do passo i
  n    <- min(length(f$pred), length(act))
  err  <- f$pred[1:n] - act[1:n]
  data.frame(
    variable = f$variable,
    horizon  = f$horizon,
    MSFE     = mean(err^2, na.rm = TRUE),
    RMSFE    = sqrt(mean(err^2, na.rm = TRUE)),
    MAE      = mean(abs(err), na.rm = TRUE)
  )
}))

print(results, digits = 4)

# Benchmark: Random Walk

rw_results <- do.call(rbind, lapply(forecasts, function(f) {
  act <- actual[, f$variable]
  h   <- f$horizon
  n   <- length(act)
  
  # random walk: previsão é o último valor observado antes da janela OOS
  # act está em escala padronizada, então RW pred = act[t], target = act[t+h]
  rw_pred <- head(act, n - h)
  rw_act  <- tail(act, n - h)
  rw_err  <- rw_pred - rw_act
  
  data.frame(
    variable = f$variable,
    horizon  = f$horizon,
    MSFE_RW  = mean(rw_err^2, na.rm = TRUE)
  )
}))

# MSFE relativo ao RW (< 1 = bate o benchmark)
results_rel <- merge(results, rw_results, by = c("variable", "horizon"))
results_rel$MSFE_rel <- results_rel$MSFE / results_rel$MSFE_RW

print(results_rel[, c("variable","horizon","MSFE","MSFE_RW","MSFE_rel","RMSFE")],
      digits = 4)

# MSFE_rel < 1 significa que o 2SRR bate o random walk — que é o resultado esperado para IPCA e EMBI em horizontes curtos (h=1,3) e mais incerto para h=12. Cole o output aqui para interpretar os resultados finais.

# Teste Diebold-Mariano para comparar MSFE do 2SRR com o RW

library(forecast)

dm_results <- do.call(rbind, lapply(forecasts, function(f) {
  act  <- actual[, f$variable]
  h    <- f$horizon
  n    <- min(length(f$pred), length(act))

  e1   <- f$pred[1:n] - act[1:n]          # erro 2SRR
  rw   <- head(act, n - h + 0)            # RW: pred = valor anterior
  # alinha RW com o mesmo n
  e2   <- head(act, n) - c(NA, head(act, n-1))
  e2[1] <- 0

  dm <- tryCatch(
    dm.test(e1, e2, alternative = "less", h = h, power = 2),
    error = function(e) NULL
  )

  data.frame(
    variable  = f$variable,
    horizon   = h,
    MSFE_rel  = mean(e1^2, na.rm=TRUE) / mean(e2^2, na.rm=TRUE),
    DM_stat   = if (!is.null(dm)) dm$statistic else NA,
    p_value   = if (!is.null(dm)) dm$p.value   else NA
  )
}))

print(dm_results, digits = 3)