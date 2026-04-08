# =============================================================================================================
# INITIALIZATION & DIRECTORY SETUP
# =============================================================================================================

rm(list = ls()) # Clear environment

# INSERT YOUR PATH HERE
wd = 'C:/Users/felip/OneDrive/Documentos/tcc/Felipe_Dornelles_tcc/'
setwd(wd)

# Define paths relative to the working directory
paths <- list(program = "00_program",
              data = "10_data",
              tools = "20_tools",
              functions = "20_tools/functions",
              output = "30_output",
              results = "40_results",
              tests = "50_tests") # <-- Nova pasta de testes adicionada

# =============================================================================================================
# PACKAGE MANAGEMENT
# =============================================================================================================

myPKGs <- c('dplyr', 'randomForest', 'mboost', 'e1071', 'readr', 'GA', 'pracma', 
            'doParallel', 'foreach', 'glmnet', 'timeSeries', 'fGarch', 'matrixcalc')

InstalledPKGs <- names(installed.packages()[,'Package'])
InstallThesePKGs <- myPKGs[!myPKGs %in% InstalledPKGs]
if (length(InstallThesePKGs) > 0) install.packages(InstallThesePKGs, repos = "http://cran.us.r-project.org")

# Load libraries
invisible(lapply(myPKGs, library, character.only = TRUE))

# =============================================================================================================
# LOAD TOOLS & FUNCTIONS
# =============================================================================================================

# Nathalia's data functions
source(paste(paths$functions, '00_Nathalia_functions.R', sep='/'))

# Coulombe's Tools
source(paste(paths$tools, 'EM_sw.R', sep='/'))
source(paste(paths$tools, 'ICp2.R', sep='/'))
source(paste(paths$tools, 'Xgenerators_v190127.R', sep='/'))

# Coulombe's TVP Functions
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
# WRAPPER FUNCTIONS (Bridging Nathalia's loop with Coulombe's models)
# =============================================================================================================

# Wrapper for 2-Step Ridge Regression (2SRR) - Coulombe's Type 2
func_coulombe_2srr = function(df, horizon, variable, lag_orders, reoptimize_hyperparameters=FALSE, type, best, ...) {
  
  # For Coulombe's models, we default to the maximum lag evaluated or the optimized one
  lag = ifelse(is.null(best$lag), max(lag_orders), best$lag)
  
  # Prepare data using Nathalia's Data-Rich ARDI setup
  prep_data = dataprep(df, horizon, variable, lag=lag, K=10, dataset="B0_ARDI")
  Xin = prep_data$Xin
  yin = prep_data$yin
  Xout = prep_data$Xout
  
  # Coulombe's lambda candidates
  lambdavec = exp(pracma::linspace(-6, 20, n=15))
  
  if (reoptimize_hyperparameters == TRUE) {
    best = list()
    best$lag = lag
    
    # TVPRR_cosso Type 2 is the 2SRR model
    out = TVPRR_cosso(X=Xin, y=yin, type=2, lambdavec=lambdavec, 
                      oosX=Xout, kfold=5, silent=1)
    
    best$lambda1 = out$lambda1
    best$lambda2 = 0.1 # default from Coulombe
    pred = out$fcast
  } else {
    # Use saved hyperparameters to save time
    out = TVPRR_cosso(X=Xin, y=yin, type=2, lambdavec=best$lambda1, lambda2=best$lambda2,
                      oosX=Xout, kfold=5, silent=1)
    pred = out$fcast
  }
  
  return(list(pred=as.numeric(pred), best=best))
}

# Wrapper for Multistep Ridge Regression, Sparse TVPs (MSRRs) - Coulombe's Type 3
func_coulombe_msrrs = function(df, horizon, variable, lag_orders, reoptimize_hyperparameters=FALSE, type, best, ...) {
  lag = ifelse(is.null(best$lag), max(lag_orders), best$lag)
  prep_data = dataprep(df, horizon, variable, lag=lag, K=10, dataset="B0_ARDI")
  Xin = prep_data$Xin
  yin = prep_data$yin
  Xout = prep_data$Xout
  lambdavec = exp(pracma::linspace(-6, 20, n=15))
  
  if (reoptimize_hyperparameters == TRUE) {
    best = list()
    best$lag = lag
    out = TVPRR(X=Xin, y=yin, type=3, lambdavec=lambdavec, oosX=Xout, kfold=5, silent=1)
    best$lambda1 = out$lambda1
    pred = out$fcast
  } else {
    out = TVPRR(X=Xin, y=yin, type=3, lambdavec=best$lambda1, oosX=Xout, kfold=5, silent=1)
    pred = out$fcast
  }
  return(list(pred=as.numeric(pred), best=best))
}

# =============================================================================================================
# FORECASTING LOOP (SMOKE TEST MODE)
# =============================================================================================================

# Load the cleaned dataset generated by 01_data_prep.R
load(paste(paths$data, "df.rda", sep='/'))
df = as.matrix(scale(df))

# Variables of interest
UNRATE <- "SEADE12_TDTGSP12"
IPCA <- "PRECOS12_IPCA12"
SPREAD <- "JPM366_EMBI366"

# <--- MODIFICAÇÕES DE TESTE AQUI --->
variable <-  c(IPCA)         # Testando apenas 1 variável
horizon_list <- c(1)         # Testando apenas 1 horizonte
lag_orders <- c(1:6)         

set.seed(1234)

# Model Selection
model_name = "test_coulombe_2srr" # Nome alterado para evitar confusão
model_function = func_coulombe_2srr 
forecasts_list = list()

# Window Setup
start_window = 155 
end_window = 156             # <--- TRAVADO EM 156 PARA RODAR APENAS 2 PASSOS --->
nwindows = end_window - start_window + 1 
actual_values = df[start_window:end_window, ]
window_type = "expanding"

for(v in variable){
  for(h in horizon_list){
    
    # The train set should be h steps behind the test set
    window_sup <- start_window - h - 1
    model_list = list()
    best = list() # Reset best params for new variable/horizon
    
    cat(sprintf("TEST: Estimating %s for horizon %d...\n", v, h))
    
    for (i in 1:nwindows){
      reoptimize_hyperparameters = FALSE
      
      if (window_type=="expanding"){
        Df = df[1:(window_sup+i),]
      } else { # rolling
        Df = df[i:(window_sup+i),]
      }
      
      # Re-optimize hyperparameters every 12 months
      if (i %% 12 == 1) {
        reoptimize_hyperparameters = TRUE
        cat(sprintf("  -> Re-optimizing at step %d\n", i))
      }
      
      # Run Model
      model = model_function(Df, horizon=h, variable=v, lag_orders=lag_orders, 
                             reoptimize_hyperparameters=reoptimize_hyperparameters, type="cv", best=best)
      
      model_list[[i]] <- list(pred=model$pred)
      best <- model$best # carry forward the hyperparams
    }
    
    # Extract forecasts
    forecast = head(unlist(lapply(model_list, function(x) x$pred)), nwindows)
    forecasts_list[[length(forecasts_list)+1]] <- list(variable=v, horizon=h, model=model_name, pred=forecast)
  }
}

# Save Results in the TESTS folder instead of output
saveRDS(forecasts_list, file = paste(paths$tests, paste0(model_name, ".rda"), sep = "/"))
saveRDS(actual_values, file = paste(paths$tests, "test_actual_values.rda", sep = "/"))

cat("Teste 02 completo! Results saved to 50_tests folder.\n")

# Carregando os arquivos pra revisao:
meus_testes_coulombe <- readRDS("~/tcc/Felipe_Dornelles_tcc/50_tests/test_coulombe_2srr.rda")
meus_valores_reais <- readRDS("~/tcc/Felipe_Dornelles_tcc/50_tests/test_actual_values.rda")
print(meus_valores_reais)
str(meus_testes_coulombe)

#verificar o arquivo do modelo do coloumbe se está dentro do desvio e de acordo com os dados reais. 
