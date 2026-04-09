# ============================================================
# 01_data_download.R
#
# Montagem de séries macroeconômicas brasileiras
# com mais dados, de modo a similarizar com o Coulombe, mas
# para usar no modelo TVP-2SRR, seguindo a lógica da Nathalia:
#   - targets (V1..V5): variáveis que quero prever
#   - painel auxiliar: séries que viram fatores PCA (os M's)
#
# Fontes:
#   - IPEADATA via {ipeadatar} (maioria das séries)
#   - BCB via {GetBCBData} (desemprego PME 1996–2011)
# ============================================================

rm(list = ls())

myPKGs <- c('dplyr', 'ipeadatar', 'readxl', 'lubridate', 'urca', 'tidyr', 'tseries')

InstalledPKGs      <- names(installed.packages()[, 'Package'])
InstallThesePKGs   <- myPKGs[!myPKGs %in% InstalledPKGs]
if (length(InstallThesePKGs) > 0)
  install.packages(InstallThesePKGs, repos = "http://cran.us.r-project.org")

invisible(lapply(myPKGs, library, character.only = TRUE))

library(ipeadatar)
library(GetBCBData)
library(dplyr)
library(tidyr)
library(purrr)
library(zoo)

# ============================================================
# SETUP — working directory, paths e pasta datada
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

# Cria pasta dentro de 10_data/ — ex.: 10_data/data_04_08_2026
run_folder <- file.path(paths$data, format(Sys.Date(), "data_%m_%d_%Y"))
if (!dir.exists(run_folder)) dir.create(run_folder)

cat(sprintf("Run folder: %s\n", run_folder))

start_date <- as.Date("1996-01-01")
end_date   <- as.Date("2025-12-01")

# ============================================================
# 1. CATÁLOGO DE SÉRIES — IPEADATA
# ------------------------------------------------------------
# Transformações:
#   "dl"  = diff(log(x))  → variação % mensal (séries em nível)
#   "d"   = diff(x)       → variação em p.p.  (taxas já em %)
#   "l"   = log(x)        → nível log
#   "niv" = x             → sem transformação (já estacionária)
# ============================================================

# ---- Grupo 1: Variáveis-alvo (os V's do grid M × V × H) ----
# O desemprego está separado — vem do BCB + IPEADATA (emenda).
targets_catalog <- tribble(
  ~code,               ~name,      ~transform,
  "BM12_PIB12",        "PIB",      "dl",
  "PRECOS12_IPCA12",   "IPCA",     "dl",
  "BM12_TJOVER12",     "SELIC",    "d",
  "BM12_ERC12",        "CAMBIO",   "dl"
)

# ---- Grupo 2: Painel auxiliar (alimenta o PCA → fatores M) ----
# Estas séries não são previstas — entram como regressores.
# O PCA extrai K=1..4 fatores delas para enriquecer o modelo.
panel_catalog <- tribble(
  ~code,                  ~name,           ~transform,

  # Atividade econômica
  "PIMPF12_QTIG12",       "PIM_geral",     "dl",
  "PIMPF12_QTBK12",       "PIM_bkap",      "dl",
  "PIMPF12_QTBCD12",      "PIM_bdur",      "dl",
  "PIMPF12_QTBI12",       "PIM_bint",      "dl",
  "BM12_CEEI12",          "EnergyConsump", "dl",
  "PMC12_VVTOT12",        "RetailSales",   "dl",

  # Mercado de trabalho
  "MTE12_SALDON12",       "CAGED_net",     "niv",
  "PMEN12_RRME12",        "RealWage",      "dl",

  # Preços e inflação
  "IGP12_IGPDIG12",       "IGPDI",         "dl",
  "IGP12_IPADIG12",       "IPA",           "dl",
  "IGP12_INCCD12",        "INCC",          "dl",
  "BM12_IPCAEXP1212",     "IPCA_exp",      "d",

  # Política monetária e crédito
  "BM12_M1MN12",          "M1",            "dl",
  "BM12_CRLIN12",         "Credit",        "dl",
  "BM12_SPREAD12",        "Spread",        "d",

  # Setor externo
  "FUNCEX12_XVTOT12",     "Exports",       "dl",
  "FUNCEX12_MVTOT12",     "Imports",       "dl",
  "BM12_RESERVAS12",      "FXReserves",    "dl",
  "JPM366_EMBI366",       "EMBI",          "d",

  # Mercado financeiro
  "GM366_IBVSP366",       "Ibovespa",      "dl",
  "BM12_ERREF12",         "RealFX",        "dl",

  # Setor fiscal
  "BM12_RNPSP12",         "PrimBalance",   "d",
  "BM12_DLSPN12",         "NetDebt_GDP",   "d"
)

# ============================================================
# 2. FUNÇÃO GENÉRICA DE DOWNLOAD — IPEADATA
# ------------------------------------------------------------
# Séries territoriais retornam múltiplas linhas por data
# (uma por UF/região). Filtra uname == "Brasil" ou vazio
# para manter apenas o agregado nacional antes de qualquer
# operação de data — esse era o motivo dos 17 FAILEDs.
# ============================================================

fetch_ipea <- function(code, name, start, end) {
  cat(sprintf("  [IPEA] %-30s %-16s ... ", code, name))
  tryCatch({
    df <- ipeadata(code)

    # Filtra agregado nacional em séries territoriais
    if ("uname" %in% names(df)) {
      nacional <- df |>
        filter(uname %in% c("", "Brasil") | is.na(uname))
      if (nrow(nacional) > 0) df <- nacional
    }

    df <- df |>
      filter(date >= start, date <= end) |>
      select(date, value) |>
      rename(!!name := value)

    cat(sprintf("%d obs.\n", nrow(df)))
    df
  }, error = function(e) {
    cat(sprintf("FAILED: %s\n", conditionMessage(e)))
    NULL
  })
}

# ============================================================
# 3. DESEMPREGO — EMENDA BCB (PME) + IPEADATA (PNAD-C)
# ------------------------------------------------------------
# Não existe série mensal de desemprego no IPEADATA que cubra
# 1996–2025 sem lacunas. Solução padrão na literatura macro BR:
#   - 1996–2011: PME/IBGE via BCB (código 24369)
#   - 2012–2025: PNADC12_TDESOCMD12 (mensal dessaz., IPEADATA)
# Aplico um offset de nível na PME para alinhar ao nível da
# PNAD-C, calculado pela diferença de médias no overlap
# jan/2012–dez/2015 — evita quebra estrutural artificial.
# use.memoise = FALSE evita erro 404 por cache corrompido.
# ============================================================

fetch_unemployment <- function(start, end) {
  cat("  Baixando desemprego (emenda PME/BCB + PNAD-C/IPEA)...\n")

  # PME via BCB — cobre até 2016, mas uso só até dez/2011
  cat("    [BCB  ] PME unemployment rate (code 24369) ... ")
  df_pme <- tryCatch({
    df <- gbcbd_get_series(
      id          = c(UNEMP_PME = 24369),
      first.date  = start,
      last.date   = as.Date("2011-12-01"),
      format.data = "wide",
      use.memoise = FALSE       # desativa cache corrompido
    ) |>
      rename(date = ref.date) |>
      select(date, DESEMPREGO = UNEMP_PME)
    cat(sprintf("%d obs. (PME up to Dec/2011)\n", nrow(df)))
    df
  }, error = function(e) {
    cat(sprintf("FAILED: %s\n", e$message)); NULL
  })

  # PNAD-C mensal dessaz. — cobre jan/2012 em diante
  cat("    [IPEA ] PNADC12_TDESOCMD12 (monthly seasonally adj.) ... ")
  df_pnadc <- tryCatch({
    raw <- ipeadata("PNADC12_TDESOCMD12")
    # Filtra nacional se série tiver dimensão territorial
    if ("uname" %in% names(raw)) {
      nac <- raw |> filter(uname %in% c("", "Brasil") | is.na(uname))
      if (nrow(nac) > 0) raw <- nac
    }
    df <- raw |>
      filter(date >= as.Date("2012-01-01"), date <= end) |>
      select(date, value) |>
      rename(DESEMPREGO = value)
    cat(sprintf("%d obs. (PNAD-C from Jan/2012)\n", nrow(df)))
    df
  }, error = function(e) {
    cat(sprintf("FAILED: %s\n", e$message)); NULL
  })

  # Calcula offset de nível no período de sobreposição (2012–2015)
  if (!is.null(df_pme) && !is.null(df_pnadc)) {
    cat("    Computing level offset for PME -> PNAD-C splice...\n")

    pme_overlap <- tryCatch(
      gbcbd_get_series(
        id          = c(PME_ext = 24369),
        first.date  = as.Date("2012-01-01"),
        last.date   = as.Date("2015-12-01"),
        format.data = "wide",
        use.memoise = FALSE     # mesmo ajuste
      ) |>
        rename(date = ref.date) |>
        pull(PME_ext),
      error = function(e) NULL
    )

    pnadc_overlap <- df_pnadc |>
      filter(date >= as.Date("2012-01-01"),
             date <= as.Date("2015-12-01")) |>
      pull(DESEMPREGO)

    if (!is.null(pme_overlap) && length(pme_overlap) > 0 &&
        length(pnadc_overlap) > 0) {
      offset <- mean(pnadc_overlap, na.rm = TRUE) -
                mean(pme_overlap,   na.rm = TRUE)
      cat(sprintf("    Level offset applied: %+.2f p.p.\n", offset))
      df_pme <- df_pme |> mutate(DESEMPREGO = DESEMPREGO + offset)
    }

    df_out <- bind_rows(df_pme, df_pnadc) |>
      arrange(date) |>
      distinct(date, .keep_all = TRUE)

    cat(sprintf("    -> DESEMPREGO total: %d obs. (%s to %s)\n",
                nrow(df_out),
                format(min(df_out$date), "%b/%Y"),
                format(max(df_out$date), "%b/%Y")))
    return(df_out)
  }

  # Fallback: retorna o que conseguiu baixar
  bind_rows(df_pme, df_pnadc) |> arrange(date)
}

# ============================================================
# 4. DOWNLOAD — TARGETS
# ============================================================

cat("\n=== Downloading targets (4 via IPEADATA + unemployment splice) ===\n")

targets_list <- map2(
  targets_catalog$code,
  targets_catalog$name,
  ~fetch_ipea(.x, .y, start_date, end_date)
)
names(targets_list) <- targets_catalog$name
targets_list        <- compact(targets_list)

targets_list$DESEMPREGO <- fetch_unemployment(start_date, end_date)

# ============================================================
# 5. DOWNLOAD — PAINEL AUXILIAR
# ============================================================

cat("\n=== Downloading auxiliary panel (~23 series for PCA) ===\n")

panel_list <- map2(
  panel_catalog$code,
  panel_catalog$name,
  ~fetch_ipea(.x, .y, start_date, end_date)
)
names(panel_list) <- panel_catalog$name
panel_list        <- compact(panel_list)

cat(sprintf("\n  -> %d/%d series downloaded successfully.\n",
            length(panel_list), nrow(panel_catalog)))

# ============================================================
# 6. MONTA O DATAFRAME WIDE
# ------------------------------------------------------------
# O full_join une todas as séries pela coluna date.
# Séries com frequência diária (EMBI, Ibovespa) geram muitas
# linhas extras — o filter final garante que só ficam as datas
# mensais dentro do período de estudo.
# ============================================================

cat("\n=== Building wide panel ===\n")

all_series <- c(targets_list, panel_list)

df_wide <- all_series |>
  reduce(full_join, by = "date") |>
  arrange(date) |>
  filter(date >= start_date, date <= end_date)

# Séries diárias (EMBI, Ibovespa) precisam ser colapsadas para
# frequência mensal antes do join — calcula média mensal
collapse_to_monthly <- function(df, col_name) {
  df |>
    mutate(date = as.Date(format(date, "%Y-%m-01"))) |>
    group_by(date) |>
    summarise(!!col_name := mean(.data[[col_name]], na.rm = TRUE),
              .groups = "drop")
}

# Recolapsa séries diárias se necessário (> 360 obs no período)
daily_series <- names(panel_list)[
  sapply(panel_list, function(x) nrow(x) > 400)
]

if (length(daily_series) > 0) {
  cat(sprintf("  Collapsing %d daily series to monthly frequency: %s\n",
              length(daily_series),
              paste(daily_series, collapse = ", ")))

  for (nm in daily_series) {
    panel_list[[nm]] <- collapse_to_monthly(panel_list[[nm]], nm)
  }

  # Reconstrói df_wide com séries já mensais
  all_series <- c(targets_list, panel_list)
  df_wide <- all_series |>
    reduce(full_join, by = "date") |>
    arrange(date) |>
    filter(date >= start_date, date <= end_date)
}

cat(sprintf("  Raw dimensions: %d rows x %d columns\n",
            nrow(df_wide), ncol(df_wide)))
cat("  NAs per column:\n")
print(colSums(is.na(df_wide)))

# ============================================================
# 7. APLICA TRANSFORMAÇÕES ESTACIONÁRIAS
# ------------------------------------------------------------
# Séries em nível (PIB, IPCA, câmbio) → dl = diff(log)
# Taxas em % (SELIC, spread, EMBI)    → d  = diff simples
# Séries já estacionárias (CAGED)     → niv = sem transform.
# ============================================================

cat("\n=== Applying stationarity transformations ===\n")

full_catalog <- bind_rows(
  targets_catalog,
  tibble(code = NA_character_, name = "DESEMPREGO", transform = "d"),
  panel_catalog
)

apply_transform <- function(x, trf) {
  switch(trf,
    "dl"  = c(NA, diff(log(x))),
    "d"   = c(NA, diff(x)),
    "l"   = log(x),
    "niv" = x,
    x
  )
}

df_transf <- df_wide
for (i in seq_len(nrow(full_catalog))) {
  nm  <- full_catalog$name[i]
  trf <- full_catalog$transform[i]
  if (!nm %in% names(df_transf)) next
  df_transf[[nm]] <- apply_transform(df_transf[[nm]], trf)
}

# Remove a primeira linha (NA gerada pela diferenciação)
df_transf <- df_transf |> slice(-1)

cat(sprintf("  Final dimensions: %d rows x %d columns\n",
            nrow(df_transf), ncol(df_transf)))
cat(sprintf("  Period: %s to %s\n",
            format(min(df_transf$date), "%b/%Y"),
            format(max(df_transf$date), "%b/%Y")))

# ============================================================
# 8. GRID DE COMBINAÇÕES M × V × H
# ------------------------------------------------------------
# M = número de fatores PCA a incluir como regressores
# V = índice do target (1=PIB, 2=IPCA, 3=SELIC, 4=CAMBIO, 5=DESEMP)
# H = horizonte de previsão em meses à frente
# Resultado: 4 × 5 × 3 = 60 combinações → 1 modelo por linha
# ============================================================

targets_br <- list(
  V1 = "PIB",
  V2 = "IPCA",
  V3 = "SELIC",
  V4 = "CAMBIO",
  V5 = "DESEMPREGO"
)

all_options <- expand.grid(
  M = 1:4,
  V = 1:5,
  H = c(1, 2, 4)
)

cat(sprintf("\n=== Estimation grid: %d combinations (M x V x H) ===\n",
            nrow(all_options)))

# ============================================================
# 9. SALVA — dentro de 10_data/data_MM_DD_YYYY/
# ============================================================

save(df_wide,     file = file.path(run_folder, "df_wide.rda"))
save(df_transf,   file = file.path(run_folder, "df_transf.rda"))
save(targets_br,  file = file.path(run_folder, "targets_br.rda"))
save(all_options, file = file.path(run_folder, "all_options.rda"))

cat(sprintf("\n=== Files saved to %s/ ===\n", run_folder))
cat(sprintf("  df_wide.rda     -> %d raw series\n",         ncol(df_wide) - 1))
cat(sprintf("  df_transf.rda   -> %d transformed series\n", ncol(df_transf) - 1))
cat("  targets_br.rda  -> 5 targets (V1..V5)\n")
cat("  all_options.rda -> 60 combinations M x V x H\n")
