# =============================================================================
# 00_program/01_data_prep_v2.R
# Preparação de dados — adaptado de Nathalia Oreda (thesis_UFRGS)
# =============================================================================

rm(list = ls())

# ---- Caminhos ---------------------------------------------------------------
wd         <- "C:/Users/felip/OneDrive/Documentos/tcc/Felipe_Dornelles_tcc/"
PATH_DATA  <- "10_data/data_04_15_2026/"
PATH_FUNCS <- "20_tools/functions/"
PATH_TOOLS <- "20_tools/"

setwd(wd)
if (!dir.exists(PATH_DATA)) dir.create(PATH_DATA, recursive = TRUE)

# ---- Pacotes ----------------------------------------------------------------
library(dplyr)
library(ipeadatar)
library(readxl)
library(lubridate)
library(urca)
library(tidyr)
library(tseries)
library(zoo)

# ---- Funções auxiliares -----------------------------------------------------
func_files <- list.files(PATH_FUNCS, pattern = "\\.R$", full.names = TRUE)
if (length(func_files) == 0)
  func_files <- list.files(PATH_TOOLS, pattern = "\\.R$",
                           full.names = TRUE, recursive = FALSE)
invisible(lapply(func_files, source))
cat("Funções carregadas:", length(func_files), "arquivo(s)\n")

needed_fns  <- c("preparacao", "transform_singlestep")
missing_fns <- needed_fns[!sapply(needed_fns, exists)]
if (length(missing_fns) > 0)
  stop("Funções não encontradas: ", paste(missing_fns, collapse = ", "))

# ---- Log --------------------------------------------------------------------
LOG_FILE <- paste0(PATH_DATA, "prep_log.txt")
cat("", file = LOG_FILE)
log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ",
                paste(..., sep = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = LOG_FILE, append = TRUE)
}
log_msg("Iniciando 01_data_prep_v2.R")

# =============================================================================
# 1. CATÁLOGO
# =============================================================================
IPCA   <- "PRECOS12_IPCA12"
SPREAD <- "JPM366_EMBI366"

DATASET_XLSX <- paste0(PATH_DATA, "dataset.xlsx")
if (!file.exists(DATASET_XLSX))
  stop("dataset.xlsx não encontrado em: ", DATASET_XLSX)

dataset_meta <- read_excel(
  DATASET_XLSX,
  col_types = c("text","text","text","date","date","numeric","text")
)
log_msg("Catálogo: ", nrow(dataset_meta), " séries")

# =============================================================================
# 2. DOWNLOAD
# =============================================================================
log_msg("Download IPEA...")
metadados        <- metadata(dataset_meta$codigo)
metadados_spread <- metadata(SPREAD)
data_raw         <- ipeadata(metadados$code)
data_spread      <- ipeadata(metadados_spread$code)
log_msg("Download OK: ", length(unique(data_raw$code)), " séries")

# =============================================================================
# 3. WIDE FORMAT
# =============================================================================
df_wide <- data_raw %>%
  pivot_wider(names_from = "code") %>%
  select(-any_of(c("uname", "tcode")))

df_wide$date <- as.Date(df_wide$date)
df_wide      <- df_wide[order(df_wide$date), ]

# EMBI+ diário → mensal
df_spread_monthly <- data_spread %>%
  pivot_wider(names_from = "code") %>%
  select(-any_of(c("uname", "tcode"))) %>%
  mutate(date       = as.Date(date),
         year_month = format(date, "%Y-%m")) %>%
  group_by(year_month) %>%
  summarise(JPM366_EMBI366 = mean(JPM366_EMBI366, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(date = as.Date(paste0(year_month, "-01"))) %>%
  select(date, JPM366_EMBI366)

df <- left_join(df_wide, df_spread_monthly, by = "date")
stopifnot(inherits(df$date, "Date"))
log_msg("Wide: ", nrow(df), " obs x ", ncol(df),
        " colunas | classe date: ", class(df$date))

# =============================================================================
# 4. FILTRO TEMPORAL
# =============================================================================
df_filtered <- subset(df,
                      date >= as.Date("1996-01-01") &
                        date <  as.Date("2019-06-01"))
df_filtered <- df_filtered %>% select_if(~ !any(is.na(.)))
log_msg("Filtrado: ", nrow(df_filtered),
        " obs x ", ncol(df_filtered) - 1, " preditores")

# =============================================================================
# 5. CLASSIFICAÇÃO DE TIPO
# =============================================================================
non_stationary <- df

tipo <- c(4)
for (i in 2:ncol(non_stationary)) {
  col      <- non_stationary[, i]
  tem_neg  <- sum(col < 0,  na.rm = TRUE) > 0
  tem_zero <- sum(col == 0, na.rm = TRUE) > 0
  so_neg   <- tem_neg && (sum(col < 0, na.rm = TRUE) == sum(!is.na(col)))

  tipo <- append(tipo,
    if      (!tem_neg && !tem_zero)            0
    else if (!tem_neg &&  tem_zero)            2
    else if ( tem_neg && !tem_zero && !so_neg) 1
    else if ( tem_neg && !tem_zero &&  so_neg) 5
    else                                       3
  )
}
log_msg("Tipo OK: ", length(tipo) - 1, " colunas")

# =============================================================================
# 6. TESTES ADF + KPSS
# =============================================================================
test <- c(4)

for (i in 2:ncol(non_stationary)) {
  X      <- na.exclude(as.matrix(non_stationary[, i]))
  k      <- 0; j <- 0; status <- "non-stationary"; iter <- 0

  while (status == "non-stationary" && iter < 10) {
    iter     <- iter + 1
    adf_res  <- tryCatch(summary(ur.df(X, "none", lags = 12)), error = function(e) NULL)
    kpss_res <- tryCatch(kpss.test(X, null = "T"),              error = function(e) NULL)

    if (is.null(adf_res) || is.null(kpss_res)) {
      log_msg("AVISO: teste falhou col ", i, " — forçando k=1")
      k <- 1; status <- "stationary"; break
    }

    adf_ok  <- adf_res@teststat[1] <= adf_res@cval[1, 2]
    kpss_ok <- kpss_res$p.value >= 0.05

    if (adf_ok && kpss_ok) {
      status <- "stationary"
      if (j == 1) k <- 0.5
    } else {
      if (j == 0) { X <- preparacao(X, i); j <- 1
      } else      { k <- k + 1
                    X <- if (tipo[i] == 5) cresc_discreto(X) else diff(X)
                    j <- j + 1 }
    }
  }
  test <- append(test, k)
}
log_msg("ADF/KPSS OK")

# =============================================================================
# 7. MAPEAMENTO tipo + test → transformation
# =============================================================================
transformation <- c(-1)

for (i in 2:ncol(non_stationary)) {
  tr <-
    if      (tipo[i] == 0 && test[i] == 0)    4
    else if (tipo[i] == 0 && test[i] == 0.5)  2
    else if (tipo[i] == 0 && test[i] == 1)    5
    else if (tipo[i] == 0 && test[i] == 2)    6
    else if (tipo[i] == 1 && test[i] == 0)    1
    else if (tipo[i] == 1 && test[i] == 1)    2
    else if (tipo[i] == 1 && test[i] == 2)    3
    else if (tipo[i] == 2 && test[i] == 0)    1
    else if (tipo[i] == 2 && test[i] == 1)    2
    else if (tipo[i] == 3 && test[i] == 0)    1
    else if (tipo[i] == 3 && test[i] == 1)    2
    else if (tipo[i] == 5 && test[i] == 1)    7
    else                                       2

  transformation <- append(transformation, tr)
}
log_msg("Transformações mapeadas OK")

# =============================================================================
# 8. EMPACOTAR dataset
# =============================================================================
month_vec <- t(as.matrix(month(df_filtered$date)))
year_vec  <- t(as.matrix(year(df_filtered$date)))

dataset <- list(
  data           = as.matrix(df_filtered[2:ncol(df_filtered)]),
  transformation = t(as.matrix(transformation[2:ncol(df_filtered)])),
  month          = month_vec,
  year           = year_vec
)
dataset$names       <- names(as.data.frame(dataset$data))
names(dataset$data) <- NULL
dataset$transform   <- dataset$transformation

# =============================================================================
# 9. APLICAR TRANSFORMAÇÕES
# =============================================================================
transformed_dataset <- data.frame(
  matrix(NA, nrow = nrow(dataset$data), ncol = ncol(dataset$data))
)
for (i in 1:ncol(dataset$data)) {
  transformed_dataset[, i] <- tryCatch(
    transform_singlestep(dataset$data[, i], dataset$transformation[i]),
    error = function(e) {
      log_msg("ERRO transform col ", i, ": ", conditionMessage(e))
      rep(NA_real_, nrow(dataset$data))
    }
  )
}
names(transformed_dataset) <- dataset$names

# =============================================================================
# 10. REMOVER NAs INICIAIS — primeiras 2 linhas de diff
# =============================================================================
df_clean <- transformed_dataset[-c(1, 2), ] %>%
  select_if(~ !any(is.na(.)))

dates     <- df_filtered$date[-c(1, 2)]
month_out <- t(as.matrix(month(dates)))
year_out  <- t(as.matrix(year(dates)))

log_msg("Pós-transformação: ", nrow(df_clean), " obs x ",
        ncol(df_clean), " colunas")

# =============================================================================
# 11. LIMPEZA FINAL — opera em df_clean (pós-transformação)
# =============================================================================

# 11a. Séries degeneradas (zeros estruturais / var≈0)
DROP_FINAL <- c(
  "BM12_TJLP12",        # 82% zeros — TJLP extinta 2018
  "DIMAC_ECFLIQTOT12",  # var ≈ 4.5e-8 — série flat
  "SEADE12_TDOPSP12",   # 26% zeros estruturais
  "SEADE12_TDOTSP12"    # 19% zeros estruturais
)
n_before  <- ncol(df_clean)
df_clean  <- df_clean %>% select(-any_of(DROP_FINAL))
log_msg("Removidas ", n_before - ncol(df_clean), " colunas degeneradas")

# 11b. Colunas em nível absoluto — log-diferenciação segura
NIVEL_COLS    <- c("BM12_DEXGFN12", "BM12_DINGFN12",
                   "BPAG12_CF12",   "BPAG12_TC12",   "BPAG12_CK12")
nivel_present <- intersect(NIVEL_COLS, names(df_clean))
nivel_highvar <- nivel_present[
  sapply(nivel_present, function(col) var(df_clean[[col]], na.rm = TRUE) > 1000)
]

if (length(nivel_highvar) > 0) {
  log_msg("Log-diferenciando ", length(nivel_highvar), " colunas em nível: ",
          paste(nivel_highvar, collapse = ", "))
  df_clean <- df_clean %>%
    mutate(across(all_of(nivel_highvar),
                  ~ c(NA_real_, diff(log(abs(.) + 1)))))
  df_clean  <- df_clean[-1, ]   # remove linha NA gerada pelo diff
  dates     <- dates[-1]
  month_out <- t(as.matrix(month(dates)))
  year_out  <- t(as.matrix(year(dates)))
}

# 11c. Verificação final: var residual ≈ 0
pred_cols  <- setdiff(names(df_clean), IPCA)
var_resid  <- sapply(df_clean[pred_cols], var, na.rm = TRUE)
still_zero <- names(var_resid[var_resid < 1e-8])
if (length(still_zero) > 0) {
  log_msg("Var≈0 residual: ", paste(still_zero, collapse = ", "))
  df_clean <- df_clean %>% select(-any_of(still_zero))
}

log_msg("df_clean final: ", nrow(df_clean), " obs x ",
        ncol(df_clean), " colunas")

# =============================================================================
# 12. SALVAR df e df_model
# =============================================================================

# df sem date — interface idêntica ao script da Nathalia
df <- df_clean
save(df, file = paste0(PATH_DATA, "df.rda"))
log_msg("Salvo: df.rda")

# df_model com date — usado pelo 02_forecast.R
df_model <- cbind(data.frame(date = dates), df_clean)
save(df_model, dates, month_out, year_out,
     file = paste0(PATH_DATA, "df_model.rda"))
log_msg("Salvo: df_model.rda")

# CSV para inspeção (com date para facilitar conferência)
write.csv(df_model, paste0(PATH_DATA, "df.csv"), row.names = FALSE)
log_msg("Salvo: df.csv")

# =============================================================================
# 13. OBJETOS AUXILIARES DO FORECAST
#   02_forecast.R carrega 4 arquivos além de df.rda:
#     df_targets.rda   → mat_y com os alvos
#     df_panel_pca.rda → painel para fatores PCA
#     targets_br.rda   → nomes/índices dos alvos
#     all_options.rda  → grid V × H × M
# =============================================================================

# Alvos: IPCA (V1), SELIC overnight (V2), EMBI+ (V3)
# Ajuste os códigos se quiser outros alvos
ALVO_COLS <- c("PRECOS12_IPCA12", "BM12_TJOVER12", "JPM366_EMBI366")
ALVO_COLS <- intersect(ALVO_COLS, names(df_model))  # só os que existem

if (length(ALVO_COLS) == 0)
  stop("Nenhum alvo encontrado em df_model. Verifique os códigos em ALVO_COLS.")

# targets_br: lista nomeada V1, V2, ...
targets_br        <- as.list(ALVO_COLS)
names(targets_br) <- paste0("V", seq_along(ALVO_COLS))

# df_targets: date + colunas alvo
df_targets <- df_model %>% select(date, all_of(ALVO_COLS))

# df_panel_pca: date + preditores (sem os alvos)
df_panel_pca <- df_model %>% select(-all_of(ALVO_COLS))

# all_options: grid completo V × H × M
all_options <- expand.grid(
  V = seq_along(ALVO_COLS),   # 1 por alvo
  H = c(1L, 3L, 6L, 12L),    # horizontes
  M = 1L:4L                   # especificações
)

log_msg("targets_br: ", paste(ALVO_COLS, collapse = " | "))
log_msg("df_panel_pca: ", ncol(df_panel_pca) - 1, " preditores")
log_msg("all_options: ", nrow(all_options), " combinações")

save(df_targets,   file = paste0(PATH_DATA, "df_targets.rda"))
save(df_panel_pca, file = paste0(PATH_DATA, "df_panel_pca.rda"))
save(targets_br,   file = paste0(PATH_DATA, "targets_br.rda"))
save(all_options,  file = paste0(PATH_DATA, "all_options.rda"))

log_msg("Salvo: df_targets.rda | df_panel_pca.rda | targets_br.rda | all_options.rda")

# =============================================================================
# RESUMO FINAL
# =============================================================================
cat("\n========================================\n")
cat("  RESUMO — 01_data_prep_v2.R\n")
cat("========================================\n")
cat("  Observações :", nrow(df_model), "\n")
cat("  Preditores  :", ncol(df_panel_pca) - 1, "\n")       # -1 pela coluna date
cat("  Alvos       :", paste(ALVO_COLS, collapse = ", "), "\n")
cat("  Período     :", as.character(min(dates)),
    "a", as.character(max(dates)), "\n")
cat("  IPCA12      :", IPCA   %in% names(df_model), "\n")
cat("  EMBI+       :", SPREAD %in% names(df_model), "\n")
cat("  Combinações :", nrow(all_options), "\n")
cat("  CSV salvo   :", paste0(PATH_DATA, "df.csv"), "\n")
cat("  Log em      :", LOG_FILE, "\n")
cat("========================================\n")