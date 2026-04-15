# ============================================================
# 01_data_prep_v1.R
#
# Montagem do painel ampliado (~50 séries) para o modelo
# TVP-2SRR / MSRRd — versão 1 do TCC.
#
# Alterações em relação ao data_prep.R original:
#   - IBC-Br substitui o PIB trimestral como alvo principal
#   - ~20 novas séries adicionadas ao painel auxiliar
#   - ICp2 (Bai & Ng 2002) determina número ótimo de fatores
#   - Janela inicial (tau) = 120 meses (2006-01 a 2015-12)
#   - Limiar de cobertura: 0.80 na base, 0.85 após transf.
#   - Imputação apenas de lacunas curtas (maxgap = 2)
# ============================================================

rm(list = ls())

# ============================================================
# 0. PACOTES
# ============================================================

myPKGs <- c("ipeadatar", "GetBCBData", "dplyr", "tidyr",
            "purrr", "zoo")
need   <- myPKGs[!myPKGs %in% names(installed.packages()[, "Package"])]
if (length(need) > 0)
  install.packages(need, repos = "http://cran.us.r-project.org")
invisible(lapply(myPKGs, library, character.only = TRUE))

# ============================================================
# 1. PATHS
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

# ============================================================
# 2. JANELA TEMPORAL
# ============================================================

# IBC-Br só existe desde jan/2003; com 4 lags + h=4 a série
# efetiva só começa em mai/2003. Adotamos jan/2003 como início
# para ter ao menos 120 observações na primeira janela OOS
# (jan/2003 – dez/2012 = 120 meses), tornando tau = 120.
# Justificativa detalhada em 02_forecast_v1.R.

start_date <- as.Date("2003-01-01")
end_date   <- as.Date("2025-12-01")

# Parâmetros de limpeza
min_coverage_raw   <- 0.80
min_coverage_trans <- 0.85
maxgap_interp      <- 2
drop_unemployment_from_pca <- TRUE

# ============================================================
# 3. CATÁLOGO DE SÉRIES
# ============================================================

# ----------------------------------------------------------
# ALVOS (targets) — 5 séries
# ----------------------------------------------------------
# IBC-Br (SGS 24363) substitui o PIB trimestral como proxy
# mensal de atividade (disponibilidade mensal completa desde
# jan/2003; metodologia BCB baseada em PIM + PMC + PMS).

targets_catalog <- tibble::tribble(
  ~source, ~code,              ~name,         ~transform,
  "bcb",   "24363",           "IBC_BR",       "dl",   # IBC-Br
  "ipea",  "PRECOS12_IPCA12", "IPCA",         "dl",
  "ipea",  "BM12_TJOVER12",   "SELIC",        "d",
  "ipea",  "BM12_ERC12",      "CAMBIO",       "dl",
  "ipea",  "MTE12_DESOCD12",  "DESEMPREGO",   "d"
)

# ----------------------------------------------------------
# PAINEL AUXILIAR
# ----------------------------------------------------------
panel_catalog <- tibble::tribble(
  ~source, ~code,                     ~name,             ~transform,

  # --- Atividade ---
  "ipea",  "PIMPF12_QTIG12",          "PIM_geral",        "dl",
  "ipea",  "PIMPF12_QTBK12",          "PIM_bkap",         "dl",
  "ipea",  "PIMPF12_QTBCD12",         "PIM_bdur",         "dl",
  "ipea",  "PIMPF12_QTBI12",          "PIM_bint",         "dl",
  "ipea",  "BM12_CEEI12",             "EnergyConsump",    "dl",
  "ipea",  "PMC12_VVTOT12",           "RetailSales",      "dl",

  # --- Mercado de trabalho ---
  "ipea",  "PMEN12_RRME12",           "RealWage",         "dl",
  "ipea",  "MTE12_SALDON12",          "CAGED_net",        "niv",

  # --- Preços ---
  "ipea",  "IGP12_IGPDIG12",          "IGPDI",            "dl",
  "ipea",  "IGP12_IPADIG12",          "IPA",              "dl",
  "ipea",  "IGP12_INCCD12",           "INCC",             "dl",
  "ipea",  "BM12_INPAM12",            "IGPM",             "dl",
  "ipea",  "PRECOS12_IPCA15M12",      "IPCA15",           "dl",
  "ipea",  "BM12_IPCAEXP1212",        "IPCA_exp",         "d",

  # --- Condições financeiras ---
  "ipea",  "BM12_M1MN12",             "M1",               "dl",
  "ipea",  "BM12_CRLIN12",            "Credit",           "dl",
  "ipea",  "BM12_SPREAD12",           "Spread",           "d",

  # --- Setor externo ---
  "ipea",  "FUNCEX12_XVTOT12",        "Exports",          "dl",
  "ipea",  "FUNCEX12_MVTOT12",        "Imports",          "dl",
  "ipea",  "BM12_RESERVAS12",         "FXReserves",       "dl",
  "ipea",  "JPM366_EMBI366",          "EMBI",             "d",
  "ipea",  "BM12_BALANP12",           "BOP",              "d",

  # --- Mercados financeiros ---
  "ipea",  "GM366_IBVSP366",          "Ibovespa",         "dl",
  "ipea",  "BM12_ERREF12",            "RealFX",           "dl",

  # --- Fiscal ---
  "ipea",  "BM12_RNPSP12",            "PrimBalance",      "d",
  "ipea",  "BM12_DLSPN12",            "NetDebt_GDP",      "d",

  # --- Expectativas / Confiça ---
  "ipea",  "SPIINDCF",                "Conf_Industria",   "d",
  "ipea",  "BM12_EXPIB12",            "PIB_exp_Focus",    "d"

  # NOTA: Sondagem de serviços (FGV) não tem código IPEA público
  # confirmado; adicionar via GetFGVData ou leitura manual de CSV.
)

full_catalog <- bind_rows(
  targets_catalog,
  panel_catalog
)

# ============================================================
# 4. FUNÇÕES AUXILIARES
# ============================================================

fetch_ipea <- function(code, name, start, end) {
  cat(sprintf("  [IPEA] %-32s %-18s ... ", code, name))
  tryCatch({
    df <- ipeadata(code)
    if ("uname" %in% names(df)) {
      nac <- df |> filter(uname %in% c("", "Brasil") | is.na(uname))
      if (nrow(nac) > 0) df <- nac
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

fetch_bcb <- function(code, name, start, end) {
  cat(sprintf("  [BCB ] %-32s %-18s ... ", code, name))
  tryCatch({
    df <- gbcbd_get_series(
      id          = setNames(as.integer(code), name),
      first.date  = start,
      last.date   = end,
      format.data = "wide",
      use.memoise = FALSE
    ) |>
      rename(date = ref.date) |>
      select(date, all_of(name)) |>
      mutate(across(-date, as.numeric)) |>
      arrange(date)
    cat(sprintf("%d obs.\n", nrow(df)))
    df
  }, error = function(e) {
    cat(sprintf("FAILED: %s\n", conditionMessage(e)))
    NULL
  })
}

collapse_to_monthly <- function(df, col_name) {
  df |>
    mutate(date = as.Date(format(date, "%Y-%m-01"))) |>
    group_by(date) |>
    summarise(
      !!col_name := mean(.data[[col_name]], na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(!!col_name := ifelse(is.nan(.data[[col_name]]), NA_real_, .data[[col_name]]))
}

calc_coverage <- function(x) mean(!is.na(x))

safe_transform <- function(x, trf) {
  x <- as.numeric(x)
  x[is.nan(x)]      <- NA
  x[is.infinite(x)] <- NA
  out <- switch(
    trf,
    "dl"  = { x[x <= 0] <- NA; c(NA, diff(log(x))) },
    "d"   = c(NA, diff(x)),
    "l"   = { x[x <= 0] <- NA; log(x) },
    "niv" = x,
    x
  )
  out[is.nan(out)]      <- NA
  out[is.infinite(out)] <- NA
  out
}

impute_short_gaps <- function(x, maxgap = 2) {
  x <- as.numeric(x)
  x[is.nan(x)]      <- NA
  x[is.infinite(x)] <- NA
  x <- zoo::na.approx(x, na.rm = FALSE, maxgap = maxgap)
  x <- zoo::na.locf(x, na.rm = FALSE)
  x <- zoo::na.locf(x, fromLast = TRUE, na.rm = FALSE)
  x
}

drop_bad_columns <- function(df, min_coverage, keep_names = character()) {
  cov_tbl   <- sapply(df |> select(-date), calc_coverage)
  keep_cov  <- names(cov_tbl)[cov_tbl >= min_coverage]
  keep_all  <- union(keep_names, keep_cov)
  keep_all  <- intersect(c("date", keep_all), names(df))
  dropped   <- setdiff(names(df), keep_all)
  list(data = df |> select(all_of(keep_all)),
       coverage = cov_tbl,
       dropped = dropped)
}

# ============================================================
# 5. DOWNLOAD
# ============================================================

cat("\n=== Baixando alvos (targets) ===\n")

targets_list <- purrr::map2(
  targets_catalog$code,
  targets_catalog$name,
  function(code, name) {
    src <- targets_catalog$source[targets_catalog$name == name]
    if (src == "bcb")
      fetch_bcb(code, name, start_date, end_date)
    else
      fetch_ipea(code, name, start_date, end_date)
  }
)
names(targets_list) <- targets_catalog$name
targets_list        <- purrr::compact(targets_list)

cat("\n=== Baixando painel auxiliar ===\n")

panel_list <- purrr::map2(
  panel_catalog$code,
  panel_catalog$name,
  function(code, name) {
    src <- panel_catalog$source[panel_catalog$name == name]
    if (src == "bcb")
      fetch_bcb(code, name, start_date, end_date)
    else
      fetch_ipea(code, name, start_date, end_date)
  }
)
names(panel_list) <- panel_catalog$name
panel_list        <- purrr::compact(panel_list)

cat(sprintf("\n  -> %d/%d séries do painel baixadas com sucesso.\n",
            length(panel_list), nrow(panel_catalog)))

# ============================================================
# 6. COLAPSA DIÁRIAS PARA MENSAL
# ============================================================

daily_series <- names(panel_list)[
  sapply(panel_list, function(x) !is.null(x) && nrow(x) > 400)
]
if (length(daily_series) > 0) {
  cat(sprintf("  Colapsando séries diárias: %s\n",
              paste(daily_series, collapse = ", ")))
  for (nm in daily_series)
    panel_list[[nm]] <- collapse_to_monthly(panel_list[[nm]], nm)
}

# ============================================================
# 7. MONTA BASE WIDE
# ============================================================

cat("\n=== Construindo painel wide ===\n")

all_series <- c(targets_list, panel_list)

df_wide <- all_series |>
  purrr::reduce(full_join, by = "date") |>
  arrange(date) |>
  filter(date >= start_date, date <= end_date) |>
  distinct(date, .keep_all = TRUE)

cat(sprintf("  Dimensões brutas: %d x %d\n", nrow(df_wide), ncol(df_wide)))

# ============================================================
# 8. FILTRA BAIXA COBERTURA (BASE BRUTA)
# ============================================================

keep_mandatory <- targets_catalog$name

raw_filter       <- drop_bad_columns(df_wide, min_coverage_raw, keep_mandatory)
df_wide_filtered <- raw_filter$data
raw_coverage_tbl <- raw_filter$coverage
dropped_raw      <- setdiff(raw_filter$dropped, "date")

cat("\n=== Cobertura (base bruta) ===\n")
print(round(sort(raw_coverage_tbl), 3))
if (length(dropped_raw) > 0) { cat("\nDropped (bruto):\n"); print(dropped_raw) }

# ============================================================
# 9. TRANSFORMAÇÕES
# ============================================================

cat("\n=== Aplicando transformações ===\n")

df_transf <- df_wide_filtered

for (i in seq_len(nrow(full_catalog))) {
  nm  <- full_catalog$name[i]
  trf <- full_catalog$transform[i]
  if (!nm %in% names(df_transf)) next
  df_transf[[nm]] <- safe_transform(df_transf[[nm]], trf)
}

df_transf <- df_transf |> slice(-1)
cat(sprintf("  Após transformação: %d x %d\n", nrow(df_transf), ncol(df_transf)))

# ============================================================
# 10. FILTRA BAIXA COBERTURA (APÓS TRANSFORMAÇÃO)
# ============================================================

trans_filter         <- drop_bad_columns(df_transf, min_coverage_trans, keep_mandatory)
df_transf_filtered   <- trans_filter$data
trans_coverage_tbl   <- trans_filter$coverage
dropped_trans        <- setdiff(trans_filter$dropped, "date")

cat("\n=== Cobertura (após transf.) ===\n")
print(round(sort(trans_coverage_tbl), 3))
if (length(dropped_trans) > 0) { cat("\nDropped (transf.):\n"); print(dropped_trans) }

# ============================================================
# 11. IMPUTA LACUNAS CURTAS
# ============================================================

cat("\n=== Imputando lacunas curtas ===\n")

df_model <- df_transf_filtered |>
  arrange(date) |>
  mutate(across(-date, ~ impute_short_gaps(.x, maxgap = maxgap_interp)))

row_na <- rowSums(is.na(df_model |> select(-date)))
df_model <- df_model[row_na < (ncol(df_model) - 2), ]

cat(sprintf("  Dimensões do modelo: %d x %d\n", nrow(df_model), ncol(df_model)))
cat("  NAs restantes por coluna:\n")
print(colSums(is.na(df_model)))

# ============================================================
# 12. OBJETOS FINAIS
# ============================================================

targets_br   <- as.list(setNames(targets_catalog$name, paste0("V", seq_along(targets_catalog$name))))
target_names <- unname(unlist(targets_br))
panel_vars   <- setdiff(names(df_model), c("date", target_names))

if (drop_unemployment_from_pca && "DESEMPREGO" %in% panel_vars)
  panel_vars <- setdiff(panel_vars, "DESEMPREGO")

df_targets   <- df_model |> select(date, any_of(target_names))
df_panel_pca <- df_model |> select(date, all_of(panel_vars))

cat(sprintf("\n=== Objetos finais ===\n"))
cat(sprintf("  Alvos (targets)    : %d séries\n", ncol(df_targets) - 1))
cat(sprintf("  Painel PCA         : %d séries\n", ncol(df_panel_pca) - 1))

# ============================================================
# 13. GRID M x V x H
# ============================================================

all_options <- expand.grid(
  M = 1:4,
  V = seq_along(target_names),
  H = c(1, 2, 4)
)
cat(sprintf("  Grid de estimativa : %d combinações (M x V x H)\n", nrow(all_options)))

# ============================================================
# 14. RELATÓRIO DE QUALIDADE
# ============================================================

series_status <- tibble::tibble(
  variable             = setdiff(names(df_model), "date"),
  in_model             = TRUE,
  coverage_raw         = raw_coverage_tbl[variable],
  coverage_transformed = trans_coverage_tbl[variable],
  remaining_na         = colSums(is.na(df_model[variable]))
)

dropped_status <- tibble::tibble(
  variable = setdiff(unique(c(dropped_raw, dropped_trans)), "date"),
  in_model = FALSE
)

# ============================================================
# 15. SALVA
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

cat(sprintf("\n=== Arquivos salvos em %s ===\n", run_folder))
cat(sprintf("  df_model         : %d x %d\n", nrow(df_model), ncol(df_model)))
cat(sprintf("  df_targets       : %d séries\n",  ncol(df_targets) - 1))
cat(sprintf("  df_panel_pca     : %d séries\n",  ncol(df_panel_pca) - 1))
cat(sprintf("  all_options      : %d combinações\n", nrow(all_options)))
