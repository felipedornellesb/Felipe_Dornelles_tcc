# ============================================================
# 01_02_data_prep.R
#
# Montagem de séries macroeconômicas brasileiras
# para uso no modelo TVP-2SRR (Coulombe 2022).
#
#   targets (V1..V5) : variáveis previstas
#   painel auxiliar  : séries que viram fatores PCA (os M's)
#
# Fontes:
#   - IPEADATA via {ipeadatar}
#   - BCB via {GetBCBData}  (desemprego PME 1996–2011)
#
# Correções desta versão vs. 01_01_data_prep.R:
#   1) safe_transform() — protege log(x<=0) e limpa NaN/Inf
#   2) drop_bad_columns() — remove séries com cobertura < mínimo
#   3) impute_short_gaps() — na.approx só para lacunas curtas
#   4) collapse_to_monthly() — limpa NaN após mean(na.rm=TRUE)
#   5) df_panel_pca separado de df_targets (DESEMPREGO fora do PCA)
#   6) relatório de qualidade + CSVs auxiliares salvos
# ============================================================

rm(list = ls())

# ============================================================
# PACKAGES
# ============================================================

myPKGs <- c("dplyr", "ipeadatar", "readxl", "lubridate",
            "urca", "tidyr", "tseries", "purrr", "zoo", "GetBCBData")

InstalledPKGs    <- names(installed.packages()[, "Package"])
InstallThesePKGs <- myPKGs[!myPKGs %in% InstalledPKGs]
if (length(InstallThesePKGs) > 0)
  install.packages(InstallThesePKGs, repos = "http://cran.us.r-project.org")

invisible(lapply(myPKGs, library, character.only = TRUE))

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

run_folder <- file.path(paths$data, format(Sys.Date(), "data_%m_%d_%Y"))
if (!dir.exists(run_folder)) dir.create(run_folder, recursive = TRUE)

cat(sprintf("Run folder: %s\n", run_folder))

start_date <- as.Date("1996-01-01")
end_date   <- as.Date("2025-12-01")

# Parâmetros de limpeza
min_coverage_raw   <- 0.80   # cobertura mínima na base em nível
min_coverage_trans <- 0.85   # cobertura mínima após transformação
maxgap_interp      <- 2      # imputa apenas lacunas de até 2 meses
drop_unemployment_from_pca <- TRUE   # DESEMPREGO não entra no PCA

# ============================================================
# 1. CATÁLOGO DE SÉRIES
# ============================================================

# Transformações:
#   "dl"  = diff(log(x))  → variação % mensal
#   "d"   = diff(x)       → variação em p.p.
#   "l"   = log(x)        → nível log
#   "niv" = x             → sem transformação

targets_catalog <- tribble(
  ~code,               ~name,      ~transform,
  "BM12_PIB12",        "PIB",      "dl",
  "PRECOS12_IPCA12",   "IPCA",     "dl",
  "BM12_TJOVER12",     "SELIC",    "d",
  "BM12_ERC12",        "CAMBIO",   "dl"
)

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

full_catalog <- bind_rows(
  targets_catalog,
  tibble(code = NA_character_, name = "DESEMPREGO", transform = "d"),
  panel_catalog
)

# ============================================================
# 2. FUNÇÕES AUXILIARES
# ============================================================

# --- 2a. Download genérico IPEADATA ---
# Filtra agregado nacional antes de qualquer operação de data,
# evitando duplicatas de séries territoriais (um por UF).
fetch_ipea <- function(code, name, start, end) {
  cat(sprintf("  [IPEA] %-30s %-16s ... ", code, name))
  tryCatch({
    df <- ipeadata(code)

    if ("uname" %in% names(df)) {
      nacional <- df |>
        filter(uname %in% c("", "Brasil") | is.na(uname))
      if (nrow(nacional) > 0) df <- nacional
    }

    df <- df |>
      filter(date >= start, date <= end) |>
      select(date, value) |>
      mutate(value = as.numeric(value)) |>
      rename(!!name := value) |>
      arrange(date)

    cat(sprintf("%d obs.\n", nrow(df)))
    df
  }, error = function(e) {
    cat(sprintf("FAILED: %s\n", conditionMessage(e)))
    NULL
  })
}

# --- 2b. Colapsa série diária para mensal ---
# Limpa NaN que surge quando all(is.na) em mean(na.rm=TRUE).
collapse_to_monthly <- function(df, col_name) {
  df |>
    mutate(date = as.Date(format(date, "%Y-%m-01"))) |>
    group_by(date) |>
    summarise(
      !!col_name := mean(.data[[col_name]], na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(!!col_name := ifelse(is.nan(.data[[col_name]]), NA_real_,
                                .data[[col_name]]))
}

# --- 2c. Cobertura de uma coluna ---
calc_coverage <- function(x) mean(!is.na(x))

# --- 2d. Transformação estacionária com guarda contra valores inválidos ---
# Problemas corrigidos vs. apply_transform() do 01_01_data_prep:
#   - NaN e Inf são convertidos em NA antes de qualquer cálculo
#   - log() recebe apenas valores > 0; negativos/zeros viram NA
#   - NaN/Inf resultantes da transformação também são limpos
safe_transform <- function(x, trf) {
  x <- as.numeric(x)
  x[is.nan(x)]      <- NA_real_
  x[is.infinite(x)] <- NA_real_

  out <- switch(
    trf,
    "dl" = {
      x[x <= 0] <- NA_real_
      c(NA_real_, diff(log(x)))
    },
    "d"   = c(NA_real_, diff(x)),
    "l"   = { x[x <= 0] <- NA_real_; log(x) },
    "niv" = x,
    x
  )

  out[is.nan(out)]      <- NA_real_
  out[is.infinite(out)] <- NA_real_
  out
}

# --- 2e. Remove colunas com cobertura insuficiente ---
# keep_names: colunas obrigatórias que nunca são removidas
# (ex.: targets), mesmo com cobertura baixa — para que o
# forecast.R receba um aviso explícito, não um objeto vazio.
drop_bad_columns <- function(df, min_coverage, keep_names = character()) {
  cov_tbl  <- sapply(df |> select(-date), calc_coverage)
  keep_cov <- names(cov_tbl)[cov_tbl >= min_coverage]
  keep_all <- union(keep_names, keep_cov)
  keep_all <- intersect(c("date", keep_all), names(df))
  dropped  <- setdiff(names(df), keep_all)
  list(data    = df |> select(all_of(keep_all)),
       coverage = cov_tbl,
       dropped  = dropped)
}

# --- 2f. Imputa apenas lacunas curtas ---
# na.approx(maxgap=) preenche apenas blocos de NA com comprimento
# <= maxgap. Lacunas mais longas permanecem NA (não inventamos dados).
# na.locf preenche bordas residuais (início/fim da série).
impute_short_gaps <- function(x, maxgap = 2L) {
  x <- as.numeric(x)
  x[is.nan(x)]      <- NA_real_
  x[is.infinite(x)] <- NA_real_
  x <- zoo::na.approx(x, na.rm = FALSE, maxgap = maxgap)
  x <- zoo::na.locf(x, na.rm = FALSE)
  x <- zoo::na.locf(x, fromLast = TRUE, na.rm = FALSE)
  x
}

# ============================================================
# 3. DESEMPREGO — EMENDA BCB (PME) + IPEADATA (PNAD-C)
# ============================================================
# Não existe série mensal contínua 1996-2025 no IPEADATA.
# Solução padrão na literatura macro BR:
#   1996–2011 : PME/IBGE via BCB (código 24369)
#   2012–2025 : PNADC12_TDESOCMD12 (mensal dessaz., IPEADATA)
# Offset de nível calculado pela diferença de médias no overlap
# jan/2012–dez/2015 — evita quebra estrutural artificial.

fetch_unemployment <- function(start, end) {
  cat("  Baixando desemprego (emenda PME/BCB + PNAD-C/IPEA)...\n")

  cat("    [BCB  ] PME unemployment rate (code 24369) ... ")
  df_pme <- tryCatch({
    df <- gbcbd_get_series(
      id          = c(UNEMP_PME = 24369),
      first.date  = start,
      last.date   = as.Date("2011-12-01"),
      format.data = "wide",
      use.memoise = FALSE
    ) |>
      rename(date = ref.date) |>
      select(date, DESEMPREGO = UNEMP_PME) |>
      mutate(DESEMPREGO = as.numeric(DESEMPREGO)) |>
      arrange(date)
    cat(sprintf("%d obs. (PME até dez/2011)\n", nrow(df)))
    df
  }, error = function(e) {
    cat(sprintf("FAILED: %s\n", e$message)); NULL
  })

  cat("    [IPEA ] PNADC12_TDESOCMD12 (mensal dessaz.) ... ")
  df_pnadc <- tryCatch({
    raw <- ipeadata("PNADC12_TDESOCMD12")
    if ("uname" %in% names(raw)) {
      nac <- raw |> filter(uname %in% c("", "Brasil") | is.na(uname))
      if (nrow(nac) > 0) raw <- nac
    }
    df <- raw |>
      filter(date >= as.Date("2012-01-01"), date <= end) |>
      select(date, value) |>
      mutate(value = as.numeric(value)) |>
      rename(DESEMPREGO = value) |>
      arrange(date)
    cat(sprintf("%d obs. (PNAD-C desde jan/2012)\n", nrow(df)))
    df
  }, error = function(e) {
    cat(sprintf("FAILED: %s\n", e$message)); NULL
  })

  if (!is.null(df_pme) && !is.null(df_pnadc)) {
    cat("    Computing level offset for PME -> PNAD-C splice...\n")

    pme_overlap <- tryCatch(
      gbcbd_get_series(
        id          = c(PME_ext = 24369),
        first.date  = as.Date("2012-01-01"),
        last.date   = as.Date("2015-12-01"),
        format.data = "wide",
        use.memoise = FALSE
      ) |>
        rename(date = ref.date) |>
        pull(PME_ext),
      error = function(e) NULL
    )

    pnadc_overlap <- df_pnadc |>
      filter(date >= as.Date("2012-01-01"),
             date <= as.Date("2015-12-01")) |>
      pull(DESEMPREGO)

    if (!is.null(pme_overlap) &&
        length(pme_overlap) > 0 &&
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

  bind_rows(df_pme, df_pnadc) |>
    arrange(date) |>
    distinct(date, .keep_all = TRUE)
}

# ============================================================
# 4. DOWNLOAD — TARGETS
# ============================================================

cat("\n=== Downloading targets (4 via IPEADATA + unemployment splice) ===\n")

targets_list <- map2(
  targets_catalog$code,
  targets_catalog$name,
  ~ fetch_ipea(.x, .y, start_date, end_date)
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
  ~ fetch_ipea(.x, .y, start_date, end_date)
)
names(panel_list) <- panel_catalog$name
panel_list        <- compact(panel_list)

cat(sprintf("\n  -> %d/%d series downloaded successfully.\n",
            length(panel_list), nrow(panel_catalog)))

# ============================================================
# 6. COLAPSA SÉRIES DIÁRIAS PARA MENSAL
# ============================================================

daily_series <- names(panel_list)[
  sapply(panel_list, function(x) nrow(x) > 400)
]

if (length(daily_series) > 0) {
  cat(sprintf("  Collapsing %d daily series to monthly: %s\n",
              length(daily_series), paste(daily_series, collapse = ", ")))
  for (nm in daily_series)
    panel_list[[nm]] <- collapse_to_monthly(panel_list[[nm]], nm)
}

# ============================================================
# 7. MONTA BASE WIDE
# ============================================================

cat("\n=== Building wide panel ===\n")

all_series <- c(targets_list, panel_list)

df_wide <- all_series |>
  reduce(full_join, by = "date") |>
  arrange(date) |>
  filter(date >= start_date, date <= end_date) |>
  distinct(date, .keep_all = TRUE)

cat(sprintf("  Raw dimensions: %d rows x %d columns\n",
            nrow(df_wide), ncol(df_wide)))
cat("  NAs per column (raw):\n")
print(colSums(is.na(df_wide)))

# ============================================================
# 8. REMOVE SÉRIES COM BAIXA COBERTURA BRUTA
# ============================================================

keep_mandatory <- c("PIB", "IPCA", "SELIC", "CAMBIO", "DESEMPREGO")

raw_filter       <- drop_bad_columns(df_wide, min_coverage_raw, keep_mandatory)
df_wide_filtered <- raw_filter$data
raw_coverage_tbl <- raw_filter$coverage
dropped_raw      <- setdiff(raw_filter$dropped, "date")

cat("\n=== Coverage check (raw data) ===\n")
print(sort(raw_coverage_tbl))

if (length(dropped_raw) > 0) {
  cat("\nDropped in raw stage (coverage < ", min_coverage_raw, "):\n", sep = "")
  print(dropped_raw)
}

# ============================================================
# 9. APLICA TRANSFORMAÇÕES ESTACIONÁRIAS
# ============================================================

cat("\n=== Applying stationarity transformations ===\n")

df_transf <- df_wide_filtered

for (i in seq_len(nrow(full_catalog))) {
  nm  <- full_catalog$name[i]
  trf <- full_catalog$transform[i]
  if (!nm %in% names(df_transf)) next
  df_transf[[nm]] <- safe_transform(df_transf[[nm]], trf)
}

# Remove 1ª linha (NA gerado pela diferenciação)
df_transf <- df_transf |> slice(-1)

cat(sprintf("  After transform: %d rows x %d columns\n",
            nrow(df_transf), ncol(df_transf)))
cat(sprintf("  Period: %s to %s\n",
            format(min(df_transf$date), "%b/%Y"),
            format(max(df_transf$date), "%b/%Y")))

# ============================================================
# 10. REMOVE SÉRIES COM BAIXA COBERTURA APÓS TRANSFORMAÇÃO
# ============================================================

trans_filter         <- drop_bad_columns(df_transf, min_coverage_trans,
                                          keep_mandatory)
df_transf_filtered   <- trans_filter$data
trans_coverage_tbl   <- trans_filter$coverage
dropped_trans        <- setdiff(trans_filter$dropped, "date")

cat("\n=== Coverage check (transformed data) ===\n")
print(sort(trans_coverage_tbl))

if (length(dropped_trans) > 0) {
  cat("\nDropped after transform (coverage < ", min_coverage_trans, "):\n", sep = "")
  print(dropped_trans)
}

# ============================================================
# 11. IMPUTA APENAS LACUNAS CURTAS
# ============================================================
# Preenche buracos de até maxgap_interp meses consecutivos.
# Lacunas mais longas permanecem NA — o forecast.R já trata
# isso filtrando colunas via is.finite() em dataprep_generic().

cat("\n=== Imputing short gaps only (maxgap =", maxgap_interp, "months) ===\n")

df_model <- df_transf_filtered |>
  arrange(date) |>
  mutate(across(-date, ~ impute_short_gaps(.x, maxgap = maxgap_interp)))

# Remove linhas ainda quase totalmente vazias
row_na_count <- rowSums(is.na(df_model |> select(-date)))
df_model     <- df_model[row_na_count < (ncol(df_model) - 2L), ]

cat(sprintf("  Model base: %d rows x %d columns\n",
            nrow(df_model), ncol(df_model)))
cat("  Remaining NAs per column after short-gap imputation:\n")
print(colSums(is.na(df_model)))

# ============================================================
# 12. SEPARA TARGETS E PAINEL PCA
# ============================================================

targets_br <- list(
  V1 = "PIB",
  V2 = "IPCA",
  V3 = "SELIC",
  V4 = "CAMBIO",
  V5 = "DESEMPREGO"
)

target_names <- unname(unlist(targets_br))
panel_vars   <- setdiff(names(df_model), c("date", target_names))

# Opcionalmente exclui DESEMPREGO do PCA — metodologicamente
# preferível porque a série emendada tem nível deslocado.
if (drop_unemployment_from_pca)
  panel_vars <- setdiff(panel_vars, "DESEMPREGO")

df_targets   <- df_model |> select(date, any_of(target_names))
df_panel_pca <- df_model |> select(date, all_of(panel_vars))

cat("\n=== Final modeling objects ===\n")
cat(sprintf("  df_targets   : %d variables\n", ncol(df_targets) - 1L))
cat(sprintf("  df_panel_pca : %d variables\n", ncol(df_panel_pca) - 1L))

# ============================================================
# 13. GRID M × V × H
# ============================================================

all_options <- expand.grid(
  M = 1:4,
  V = 1:5,
  H = c(1, 2, 4)
)

cat(sprintf("\n=== Estimation grid: %d combinations (M x V x H) ===\n",
            nrow(all_options)))

# ============================================================
# 14. RELATÓRIO DE QUALIDADE
# ============================================================

series_status <- tibble(
  variable             = setdiff(names(df_model), "date"),
  in_model             = TRUE,
  coverage_raw         = raw_coverage_tbl[variable],
  coverage_transformed = trans_coverage_tbl[variable],
  remaining_na         = colSums(is.na(df_model[variable]))
)

dropped_all <- unique(c(dropped_raw, dropped_trans))
dropped_status <- if (length(dropped_all) > 0) {
  tibble(variable = setdiff(dropped_all, "date"), in_model = FALSE)
} else tibble(variable = character(), in_model = logical())

# ============================================================
# 15. SALVA — 10_data/data_MM_DD_YYYY/
# ============================================================

save(df_wide,            file = file.path(run_folder, "df_wide.rda"))
save(df_wide_filtered,   file = file.path(run_folder, "df_wide_filtered.rda"))
save(df_transf,          file = file.path(run_folder, "df_transf.rda"))
save(df_transf_filtered, file = file.path(run_folder, "df_transf_filtered.rda"))
save(df_model,           file = file.path(run_folder, "df_model.rda"))
save(df_targets,         file = file.path(run_folder, "df_targets.rda"))
save(df_panel_pca,       file = file.path(run_folder, "df_panel_pca.rda"))
save(targets_br,         file = file.path(run_folder, "targets_br.rda"))
save(all_options,        file = file.path(run_folder, "all_options.rda"))
save(series_status,      file = file.path(run_folder, "series_status.rda"))
save(dropped_status,     file = file.path(run_folder, "dropped_status.rda"))

write.csv(df_model,       file.path(run_folder, "df_model.csv"),      row.names = FALSE)
write.csv(df_targets,     file.path(run_folder, "df_targets.csv"),    row.names = FALSE)
write.csv(df_panel_pca,   file.path(run_folder, "df_panel_pca.csv"),  row.names = FALSE)
write.csv(series_status,  file.path(run_folder, "series_status.csv"), row.names = FALSE)
write.csv(dropped_status, file.path(run_folder, "dropped_status.csv"),row.names = FALSE)

cat(sprintf("\n=== Files saved to %s/ ===\n", run_folder))
cat(sprintf("  df_wide.rda            -> %d raw series\n",               ncol(df_wide) - 1L))
cat(sprintf("  df_wide_filtered.rda   -> %d filtered raw series\n",      ncol(df_wide_filtered) - 1L))
cat(sprintf("  df_transf.rda          -> %d transformed series\n",       ncol(df_transf) - 1L))
cat(sprintf("  df_transf_filtered.rda -> %d transformed+filtered series\n", ncol(df_transf_filtered) - 1L))
cat(sprintf("  df_model.rda           -> %d final model series\n",       ncol(df_model) - 1L))
cat(sprintf("  df_targets.rda         -> %d target series (V1..V5)\n",   ncol(df_targets) - 1L))
cat(sprintf("  df_panel_pca.rda       -> %d PCA panel series\n",         ncol(df_panel_pca) - 1L))
cat("  targets_br.rda         -> 5 targets (V1..V5)\n")
cat("  all_options.rda        -> 60 combinations M x V x H\n")
cat("  series_status.csv      -> quality report\n")
cat("  dropped_status.csv     -> dropped series log\n")

