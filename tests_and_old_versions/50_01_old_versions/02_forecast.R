# =============================================================================================================
# INITIALIZATION & DIRECTORY SETUP
# =============================================================================================================

rm(list = ls()) # Clear environment

# INSERT YOUR PATH HERE
wd = 'C:/Users/felip/OneDrive/Documentos/tcc/Felipe_Dornelles_tcc'
setwd(wd)

# Define paths relative to the working directory
paths <- list(program = "00_program",
              data = "10_data",
              tools = "20_tools",
              functions = "20_tools/functions",
              output = "30_output",
              results = "40_results")

# =============================================================================================================
# PACKAGE MANAGEMENT
# =============================================================================================================

myPKGs <- c('dplyr', 'randomForest', 'mboost', 'e1071', 'readr', 'GA', 'pracma', 
            'doParallel', 'foreach', 'glmnet', 'timeSeries', 'fGarch', 'matrixcalc')

InstalledPKGs <- names(installed.packages()[,'Package'])
InstallThesePKGs <- myPKGs[!myPKGs %in% InstalledPKGs]
if (length(InstallThesePKGs) > 0) install.packages(InstallThesePKGs, repos = "http://cran.us.r-project.org")

invisible(lapply(myPKGs, library, character.only = TRUE))

# =============================================================================================================
# LOAD TOOLS & FUNCTIONS
# =============================================================================================================

source(paste(paths$functions, '00_Nathalia_functions.R', sep='/'))
source(paste(paths$tools, 'EM_sw.R', sep='/'))
source(paste(paths$tools, 'ICp2.R', sep='/'))
source(paste(paths$tools, 'Xgenerators_v190127.R', sep='/'))
source(paste(paths$functions, 'dualGRRmdA_v190215.R', sep='/'))
source(paste(paths$functions, 'CVGSBHK_v181127.R', sep='/'))
source(paste(paths$functions, 'zfun_v190304.R', sep='/'))
source(paste(paths$functions, 'factor.R', sep='/'))
source(paste(paths$functions, 'TVPRRcosso_v181120.R', sep='/'))
source(paste(paths$functions, 'TVPRRcossoF_v190125.R', sep='/'))
source(paste(paths$functions, 'TVPRR_v181111.R', sep='/')) 
source(paste(paths$functions, 'fastZrot_v181125.R', sep='/'))
source(paste(paths$functions, 'CVKFMV_v190214.R', sep='/'))
source(paste(paths$functions, 'TVPRR_VARF_v190304.R', sep='/'))

# =============================================================================================================
# WRAPPER FUNCTIONS 
# =============================================================================================================

func_coulombe_2srr = function(df, horizon, variable, lag_orders, reoptimize_hyperparameters=FALSE, type, best, ...) {
  lag = ifelse(is.null(best$lag), max(lag_orders), best$lag)
  
  # Reduzido para K=6
  prep_data = dataprep(df, horizon, variable, lag=lag, K=6, dataset="B0_ARDI") 
  Xin = prep_data$Xin
  yin = prep_data$yin
  Xout = prep_data$Xout
  
  lambdavec = exp(pracma::linspace(-6, 20, n=15))
  
  # CORREÇÃO EFEITO CASCATA: Se o best estiver vazio, força a otimização
  if(is.null(best$lambda1)) {
    reoptimize_hyperparameters <- TRUE
  }
  
  result <- tryCatch({
    if (reoptimize_hyperparameters == TRUE) {
      out = TVPRR_cosso(X=Xin, y=yin, type=2, lambdavec=lambdavec, oosX=Xout, kfold=5, silent=1)
      best = list()
      best$lag = lag
      best$lambda1 = out$lambda1
      best$lambda2 = 0.1 
      pred = out$fcast
    } else {
      out = TVPRR_cosso(X=Xin, y=yin, type=2, lambdavec=best$lambda1, lambda2=best$lambda2, oosX=Xout, kfold=5, silent=1)
      pred = out$fcast
    }
    list(pred=as.numeric(pred), best=best)
    
  }, error = function(e) {
    message(paste(" -> ERRO IGNORADO VIA TRYCATCH:", e$message))
    return(list(pred=NA, best=best))
  })
  
  return(result)
}

func_coulombe_msrrs = function(df, horizon, variable, lag_orders, reoptimize_hyperparameters=FALSE, type, best, ...) {
  lag = ifelse(is.null(best$lag), max(lag_orders), best$lag)
  prep_data = dataprep(df, horizon, variable, lag=lag, K=6, dataset="B0_ARDI")
  Xin = prep_data$Xin
  yin = prep_data$yin
  Xout = prep_data$Xout
  lambdavec = exp(pracma::linspace(-6, 20, n=15))
  
  if(is.null(best$lambda1)) {
    reoptimize_hyperparameters <- TRUE
  }
  
  result <- tryCatch({
    if (reoptimize_hyperparameters == TRUE) {
      out = TVPRR(X=Xin, y=yin, type=3, lambdavec=lambdavec, oosX=Xout, kfold=5, silent=1)
      best = list()
      best$lag = lag
      best$lambda1 = out$lambda1
      pred = out$fcast
    } else {
      out = TVPRR(X=Xin, y=yin, type=3, lambdavec=best$lambda1, oosX=Xout, kfold=5, silent=1)
      pred = out$fcast
    }
    list(pred=as.numeric(pred), best=best)
    
  }, error = function(e) {
    message(paste(" -> ERRO IGNORADO VIA TRYCATCH:", e$message))
    return(list(pred=NA, best=best))
  })
  return(result)
}

# =============================================================================================================
# FORECASTING LOOP 
# =============================================================================================================

load(paste(paths$data, "df.rda", sep='/'))
df = as.matrix(scale(df))

UNRATE <- "SEADE12_TDTGSP12"
IPCA <- "PRECOS12_IPCA12"
SPREAD <- "JPM366_EMBI366"

variable <-  c(IPCA, UNRATE, SPREAD)
horizon_list <- c(1, 3, 6, 12)
lag_orders <- c(1:6) 

set.seed(1234)

model_name = "coulombe_2srr"
model_function = func_coulombe_2srr 
forecasts_list = list()

start_window = 155 
end_window = nrow(df) 
nwindows = end_window - start_window + 1 
actual_values = df[start_window:end_window, ]
window_type = "expanding"

for(v in variable){
  
  # PROTEÇÃO CONTRA VARIÁVEIS DELETADAS NO PREP
  if(!v %in% colnames(df)) {
    cat(sprintf("\n[ALERTA FATAL] Variável %s não existe no df.rda! Pulando...\n", v))
    next
  }
  
  for(h in horizon_list){
    window_sup <- start_window - h - 1
    model_list = list()
    best = list() 
    
    cat(sprintf("\nEstimating %s for horizon %d...\n", v, h))
    
    for (i in 1:nwindows){
      reoptimize_hyperparameters = FALSE
      
      if (window_type=="expanding"){
        Df = df[1:(window_sup+i),]
      } else { 
        Df = df[i:(window_sup+i),]
      }
      
      if (i %% 12 == 1) {
        reoptimize_hyperparameters = TRUE
        cat(sprintf("  -> Re-optimizing at step %d\n", i))
      }
      
      model = model_function(Df, horizon=h, variable=v, lag_orders=lag_orders, 
                             reoptimize_hyperparameters=reoptimize_hyperparameters, type="cv", best=best)
      
      if(is.na(model$pred)) {
        if(i > 1) {
          model$pred <- model_list[[i-1]]$pred
        } else {
          model$pred <- 0 
        }
      }
      
      model_list[[i]] <- list(pred=model$pred)
      best <- model$best 
    }
    
    forecast = head(unlist(lapply(model_list, function(x) x$pred)), nwindows)
    forecasts_list[[length(forecasts_list)+1]] <- list(variable=v, horizon=h, model=model_name, pred=forecast)
  }
}

saveRDS(forecasts_list, file = paste(paths$output, paste0(model_name, ".rda"), sep = "/"))
saveRDS(actual_values, file = paste(paths$output, "actual_values.rda", sep = "/"))

cat("\nForecasting complete! Results saved to 30_output folder.\n")
