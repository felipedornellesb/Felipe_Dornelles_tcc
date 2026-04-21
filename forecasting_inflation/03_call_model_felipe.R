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

pkgs <- c("here", "glmnet", "tidyverse", "ks", "expm", "DistributionUtils", "rugarch")
new  <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(new)) install.packages(new)

library(here)
setwd(here("forecasting_inflation"))

library(rugarch)
library(glmnet)
library(tidyverse)

# Funções do Medeiros (dataprep, accumulate_model, etc.)
source("functions/functions.R")

# Implementação 2SRR (Coulombe 2024) — R moderno, sem dependências externas
source("functions/tvp_ridge_functions.R")

# Engine rolling window do Medeiros
source("functions/rolling_window.R")

# ============================================================
# Parâmetros
# ============================================================

model_name <- "2SRR"
# Medeiros: 180 windows
# Felipe/Hudson: 312 windows
nwindows <- 180   # alterar para 180 ou 312 conforme desejado

# ============================================================
# Carrega data.rda
# ============================================================

load("data/data.rda")
dates <- data$date
df    <- data %>%
  select(-date) %>%
  as.matrix()
rownames(df) <- as.character(dates)

cat("Dados carregados:", nrow(df), "observacoes,", ncol(df), "variaveis\n")
cat("Periodo:", rownames(df)[1], "a", rownames(df)[nrow(df)], "\n\n")

# ============================================================
# Diagnóstico de uma janela (descomentar para testar antes do loop)
# ============================================================

if (FALSE) {
  ind_teste <- 1:nwindows
  prep <- dataprep(ind_teste, df, "CPIAUCSL", horizon = 1, nofact = TRUE)
  cat("Xin:", dim(prep$Xin), "| yin:", length(prep$yin),
      "| NAs em Xout:", any(is.na(prep$Xout)), "\n")

  t0 <- proc.time()
  r1 <- run2srr(ind_teste, df, "CPIAUCSL", horizon = 1)
  cat(sprintf("Forecast h=1 teste: %.6f | K_pca: %d | Tempo: %.1fs\n",
              r1$forecast, r1$outputs$K_pca, (proc.time() - t0)[3]))
}

# ============================================================
# Rolling window — horizontes 1 a 12
# ============================================================

model_list <- list()
t_total    <- proc.time()

for (i in 1:12) {
  cat(sprintf("[%s] Rodando horizonte h = %d ...\n",
              format(Sys.time(), "%H:%M:%S"), i))
  t0 <- proc.time()

  model_list[[i]] <- rolling_window(
    run2srr, df, nwindows + i - 1, i, "CPIAUCSL"
  )

  cat(sprintf("[%s] h = %d concluido em %.1f min\n",
              format(Sys.time(), "%H:%M:%S"), i,
              (proc.time() - t0)[3] / 60))
}

cat(sprintf("\nTotal: %.1f min\n", (proc.time() - t_total)[3] / 60))

# ============================================================
# Consolida forecasts no formato padrão do Medeiros
# (mesmo formato de Ridge.rda, RF.rda, etc.)
# ============================================================

forecasts <- Reduce(
  cbind,
  lapply(model_list, function(x) head(x$forecast, nwindows))
)

forecasts <- accumulate_model(forecasts)

if (!dir.exists("forecasts")) dir.create("forecasts")
save(forecasts, file = paste0("forecasts/", model_name, ".rda"))
cat(sprintf("\nForecasts salvos em forecasts/%s.rda\n", model_name))

# ============================================================
# Salva betas time-varying do horizonte h = 1
# (usado em 04_eval_results_felipe.R para análise dos betas)
# ============================================================

df_betas <- extract_betas_over_time(
  model_list = model_list[[1]],
  df         = df,
  nwindows   = nwindows
)

save(df_betas, file = "forecasts/betas_2SRR.rda")
cat(sprintf("Betas time-varying salvos em forecasts/betas_2SRR.rda (%d linhas)\n",
            nrow(df_betas)))

# Exporta betas de todos os horizontes em CSVs separados
for (h in 1:12) {
  df_b_h <- extract_betas_over_time(
    model_list = model_list[[h]],
    df         = df,
    nwindows   = nwindows
  )
  write.csv(df_b_h,
            file      = sprintf("results/betas_2SRR_h%d.csv", h),
            row.names = FALSE)
}
cat("Betas de todos os horizontes exportados em results/\n")

# ============================================================
# Plot rápido de verificação (igual ao 03_call_model.R)
# ============================================================

plot(tail(df[, "CPIAUCSL"], nwindows),
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