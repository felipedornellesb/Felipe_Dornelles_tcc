# ============================================================
# 03_call_model_felipe.R
#
# Roda o modelo 2SRR (Coulombe 2019) no padrão do Medeiros.
# Salva os forecasts em forecasts/2SRR.rda (mesmo formato
# dos outros modelos do repositório) e os betas time-varying
# em forecasts/betas_2SRR.rda para análise posterior.
#
# Pré-requisito:
#   - (Caso queria atualizar a base de dados): 
#           Rodar 01_get_fred_data.R para ter data/data.rda
#   - Rodar 02_random_walk_oos_y.R para ter forecasts/rw.rda
#     e forecasts/yout.rda
#
# Ordem de execução:
#   01 → 02 → 03 (outros modelos) → 03_call_model_felipe.R → 04_eval_results.R
# ============================================================
library(here)
setwd(here("forecasting_inflation"))

install.packages("pracma")

library(pracma)
library(glmnet)
library(tidyverse)

# Funções do Medeiros (dataprep, accumulate_model, etc.)
source("functions/functions.R")

# Funções do 2SRR (run2srr, extract_betas_over_time)
source("functions/felipe_functions.R")

# Engine rolling window do Medeiros
source("functions/rolling_window.R")

# Funções do Coulombe — ajustar o caminho
source("coulombe/TVPRRcosso_v181120.R")
source("coulombe/zfun_v190304.R")
source("coulombe/fastZrot_v181125.R")
source("coulombe/CVGSBHK_v181127.R")


# ============================================================
# Parâmetros
# ============================================================

model_name <- "2SRR"
nwindows   <- 312     # mesmo valor do 03_call_model.R adaptado


# ============================================================
# Carrega data.rda
# ============================================================

load("data/data.rda")
dates <- data$date
data  <- data %>%
  select(-date) %>%
  as.matrix()
rownames(data) <- as.character(dates)

cat("Dados carregados:", nrow(data), "observações,",
    ncol(data), "variáveis\n")
cat("Período:", rownames(data)[1], "a",
    rownames(data)[nrow(data)], "\n\n")


# ============================================================
# Rolling window — horizontes 1 a 12
# ============================================================

model_list <- list()

for (i in 1:12) {
  cat(sprintf("Rodando horizonte h = %d ...\n", i))

  model_list[[i]] <- rolling_window(
    run2srr,
    data,
    nwindows + i - 1,
    i,
    "CPIAUCSL"
  )

  cat(sprintf("  Horizonte %d concluído.\n", i))
}


# ============================================================
# Consolida forecasts no formato padrão do Medeiros
# (mesmo formato de Ridge.rda, RF.rda, etc.)
# ============================================================

forecasts <- Reduce(
  cbind,
  lapply(model_list, function(x) head(x$forecast, nwindows))
)

forecasts <- accumulate_model(forecasts)

save(forecasts,
     file = paste0("forecasts/", model_name, ".rda"))

cat(sprintf("\nForecasts salvos em forecasts/%s.rda\n", model_name))


# ============================================================
# Salva betas time-varying do horizonte h = 1
# (usado em 05_compare_results.R para análise dos betas)
# ============================================================

df_betas <- extract_betas_over_time(
  model_list = model_list[[1]],
  df         = data,
  nwindows   = nwindows
)

save(df_betas,
     file = "forecasts/betas_2SRR.rda")

cat("Betas time-varying salvos em forecasts/betas_2SRR.rda\n")
cat(sprintf("  Dimensões: %d linhas (janelas x variáveis)\n",
            nrow(df_betas)))


# ============================================================
# Plot rápido de verificação (igual ao 03_call_model.R)
# ============================================================

plot(tail(data[, "CPIAUCSL"], nwindows),
     type = "l",
     main = "2SRR — Forecast h=1 vs Realizado",
     ylab = "CPIAUCSL",
     xlab = "Janelas OOS")
lines(forecasts[, 1], col = "blue", lwd = 1.5)
legend("topright",
       legend = c("Realizado", "2SRR h=1"),
       col    = c("black", "blue"),
       lty    = 1, lwd = c(1, 1.5),
       bty    = "n")