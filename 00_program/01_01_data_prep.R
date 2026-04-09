# =============================================================================================================
# 01_DATA_PREP.R
# Descrição: Preparação dos mesmos dados analisados por Nathalia em seu TCC, para uso em análises subsequentes
# - Prepara os dados para análise, incluindo:
# - Importação e limpeza dos dados
# - Testes de estacionariedade e transformações necessárias para cada série
# - Exportação do dataset final para uso em análises subsequentes 
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
# PACKAGE MANAGEMENT
# =============================================================================================================

myPKGs <- c('dplyr', 'ipeadatar', 'readxl', 'lubridate', 'urca', 'tidyr', 'tseries')

InstalledPKGs      <- names(installed.packages()[, 'Package'])
InstallThesePKGs   <- myPKGs[!myPKGs %in% InstalledPKGs]
if (length(InstallThesePKGs) > 0)
  install.packages(InstallThesePKGs, repos = "http://cran.us.r-project.org")

invisible(lapply(myPKGs, library, character.only = TRUE))

source(paste(paths$functions, '00_Nathalia_functions.R', sep = '/'))

# =============================================================================================================
# 1. GET DATA
# =============================================================================================================

dataset_meta <- read_excel(
  paste(paths$data, "dataset.xlsx", sep = '/'),
  col_types = c("text", "text", "text", "date", "date", "numeric", "text")
)

UNRATE <- "SEADE12_TDTGSP12"
IPCA   <- "PRECOS12_IPCA12"
SPREAD <- "JPM366_EMBI366"

metadados        <- metadata(dataset_meta$codigo)
metadados_spread <- metadata("JPM366_EMBI366")

data        <- ipeadata(metadados$code)
data_spread <- ipeadata(metadados_spread$code)

# [BUG 5 CORRIGIDO] Força date para Date logo após pivot_wider
# pivot_wider pode criar list-column quando há múltiplos valores por data/código
df <- data %>%
  pivot_wider(names_from = "code") %>%
  select(-c(uname, tcode)) %>%
  mutate(date = as.Date(as.character(date)))

# [BUG 6 CORRIGIDO] Ordena por date já convertida, sem format desnecessário
df_sorted <- df[order(df$date), ]

# EMBI: diário -> mensal (média)
# [BUG 5 CORRIGIDO] Mesma correção aplicada ao data_spread
df_spread <- data_spread %>%
  pivot_wider(names_from = "code") %>%
  select(-c(uname, tcode)) %>%
  mutate(
    date       = as.Date(as.character(date)),
    year_month = format(date, "%Y-%m")
  ) %>%
  group_by(year_month) %>%
  summarise(JPM366_EMBI366 = mean(JPM366_EMBI366, na.rm = TRUE), .groups = "drop") %>%
  mutate(date = as.Date(paste0(year_month, "-01"))) %>%
  select(date, JPM366_EMBI366)

# Merge
df <- left_join(df_sorted, df_spread, by = "date")

# Diagnóstico: para se date ainda não for Date por algum motivo
stopifnot("Coluna date deve ser classe Date" = inherits(df$date, "Date"))

df_filtered <- df %>%
  filter(date >= as.Date("1996-01-01") & date < as.Date("2019-06-01")) %>%
  select(where(~ !any(is.na(.))))

cat(sprintf("df_filtered: %d obs x %d colunas (incluindo date)\n",
            nrow(df_filtered), ncol(df_filtered)))

# =============================================================================================================
# 2. STATIONARITY TESTS & TRANSFORMATIONS
# =============================================================================================================

non_stationary <- df_filtered

tipo <- c(4)
for (i in 2:ncol(non_stationary)) {
  col_vals <- non_stationary[[i]]
  tem_neg  <- sum(col_vals < 0,  na.rm = TRUE) > 0
  tem_zero <- sum(col_vals == 0, na.rm = TRUE) > 0
  n_total  <- sum(!is.na(col_vals))
  todo_neg <- sum(col_vals < 0,  na.rm = TRUE) == n_total
  
  if (tem_neg) {
    if (tem_zero)      tipo <- append(tipo, 3)
    else if (todo_neg) tipo <- append(tipo, 5)
    else               tipo <- append(tipo, 1)
  } else {
    if (tem_zero) tipo <- append(tipo, 2)
    else          tipo <- append(tipo, 0)
  }
}

test <- c(4)
for (i in 2:ncol(non_stationary)) {
  cat(sprintf("Testando variavel %d de %d\r", i, ncol(non_stationary)))
  X      <- na.exclude(as.matrix(non_stationary[[i]]))
  k      <- 0
  j      <- 0
  status <- "non-stationary"
  
  while (status == "non-stationary") {
    adf_test <- summary(ur.df(X, "none", lags = 12))
    
    estacionario_adf  <- adf_test@teststat[1] <= adf_test@cval[1, 2]
    estacionario_kpss <- kpss.test(X, null = "T")$p.value >= 0.05
    
    if (estacionario_adf && estacionario_kpss) {
      status <- "stationary"
      if (j == 1) k <- 0.5
    } else {
      if (j == 0) {
        X <- preparacao(X, i)
        j <- 1
      } else {
        k <- k + 1
        if (tipo[i] %in% c(0, 1, 2, 3)) X <- diff(X)
        if (tipo[i] == 5)               X <- cresc_discreto(X)
        j <- j + 1
      }
    }
    
    if (k >= 5) {
      warning(paste("Variavel", colnames(non_stationary)[i], "nao convergiu - k forcado a 2"))
      k <- 2
      break
    }
  }
  test <- append(test, k)
}
cat("\nTeste de estacionariedade concluido.\n")

transformation <- c(-1)
for (i in 2:ncol(non_stationary)) {
  t <- tipo[i]; k <- test[i]
  code <- NA
  if      (t == 0 && k == 0)   code <- 4
  else if (t == 0 && k == 0.5) code <- 2
  else if (t == 0 && k == 1)   code <- 5
  else if (t == 0 && k == 2)   code <- 6
  else if (t == 1 && k == 0)   code <- 1
  else if (t == 1 && k == 1)   code <- 2
  else if (t == 1 && k == 2)   code <- 3
  else if (t == 2 && k == 0)   code <- 1
  else if (t == 2 && k == 1)   code <- 2
  else if (t == 3 && k == 0)   code <- 1
  else if (t == 3 && k == 1)   code <- 2
  else if (t == 5 && k == 1)   code <- 7
  transformation <- append(transformation, code)
}

month     <- t(as.matrix(month(df_filtered$date)))
year      <- t(as.matrix(year(df_filtered$date)))
reference <- available_subjects()
metadados <- metadados %>% inner_join(reference, by = "code")
nome      <- names(df_filtered[2:ncol(df_filtered)])
metadados <- metadados[metadados$code %in% nome, ]

dataset <- list(
  data           = as.matrix(df_filtered[2:ncol(df_filtered)]),
  transformation = t(as.matrix(transformation[2:ncol(df_filtered)])),
  month          = month,
  year           = year,
  metadados      = metadados
)
dataset$names       <- names(as.data.frame(dataset$data))
names(dataset$data) <- NULL

# =============================================================================================================
# 3. APPLY TRANSFORMATIONS AND EXPORT
# =============================================================================================================

n_vars  <- ncol(dataset$data)
n_obs   <- nrow(dataset$data)
transformed_dataset <- data.frame(matrix(NA, nrow = n_obs, ncol = n_vars))

for (i in seq_len(n_vars)) {
  transformed_dataset[, i] <- transform_singlestep(dataset$data[, i], dataset$transformation[i])
}
names(transformed_dataset) <- dataset$names

max_nas <- max(sapply(seq_len(n_vars), function(i) {
  tid <- dataset$transformation[i]
  if (tid %in% c(3, 6))         2L
  else if (tid %in% c(2, 5, 7)) 1L
  else                           0L
}))

df <- transformed_dataset %>%
  slice(-seq_len(max_nas)) %>%
  select(where(~ sum(is.na(.)) == 0))

cat(sprintf("Dimensao final do df: %d obs x %d variaveis\n", nrow(df), ncol(df)))

save(df, file = paste(paths$data, "df.rda", sep = "/"))
cat("Preparacao concluida. Resultado salvo em 10_data/df.rda\n")