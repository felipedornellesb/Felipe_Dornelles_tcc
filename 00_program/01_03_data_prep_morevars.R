# ============================================================
# 01_data_download.R
#
# Montagem de séries macroeconômicas brasileiras
# para usar no modelo TVP-2SRR.
# ============================================================

rm(list = ls())

myPKGs <- c('dplyr', 'ipeadatar', 'GetBCBData', 'readxl', 'lubridate', 'urca', 'tidyr', 'tseries', 'purrr', 'zoo')

InstalledPKGs      <- names(installed.packages()[, 'Package'])
InstallThesePKGs   <- myPKGs[!myPKGs %in% InstalledPKGs]
if (length(InstallThesePKGs) > 0)
  install.packages(InstallThesePKGs, repos = "http://cran.us.r-project.org")

invisible(lapply(myPKGs, library, character.only = TRUE))

# ============================================================
# SETUP — working directory, paths e pasta datada
# ============================================================

wd <- "C:/Users/00341048/Downloads/Felipe_Dornelles_tcc-main/Felipe_Dornelles_tcc-main/"
setwd(wd)

paths <- list(
  program   = "00_program",
  data      = "10_data",
  tools     = "20_tools",
  functions = "20_tools/functions",
  output    = "30_output",
  results   = "40_results"
)

# REVISÃO: recursive = TRUE garante a criação das pastas-mãe se não existirem
run_folder <- file.path(paths$data, format(Sys.Date(), "data_%m_%d_%Y"))
if (!dir.exists(run_folder)) dir.create(run_folder, recursive = TRUE)

cat(sprintf("Run folder: %s\n", run_folder))

start_date <- as.Date("1996-01-01")
end_date   <- as.Date("2025-12-01")

# ============================================================
# 1. CATÁLOGO DE SÉRIES — IPEADATA
# ============================================================

targets_catalog <- tribble(
  ~code,               ~name,      ~transform,
  "BM12_PIB12",        "PIB",      "dl",
  "PRECOS12_IPCA12",   "IPCA",     "dl",
  "BM12_TJOVER12",     "SELIC",    "d",
  "BM12_ERC12",        "CAMBIO",   "dl"
)

panel_catalog <- tribble(
  ~code,                  ~name,           ~transform,
  "PIMPF12_QTIG12",       "PIM_geral",     "dl",
  "PIMPF12_QTBK12",       "PIM_bkap",      "dl",
  "PIMPF12_QTBCD12",      "PIM_bdur",      "dl",
  "PIMPF12_QTBI12",       "PIM_bint",      "dl",
  "BM12_CEEI12",          "EnergyConsump", "dl",
  "PMC12_VVTOT12",        "RetailSales",   "dl",
  "MTE12_SALDON12",       "CAGED_net",     "niv",
  "PMEN12_RRME12",        "RealWage",      "dl",
  "IGP12_IGPDIG12",       "IGPDI",         "dl",
  "IGP12_IPADIG12",       "IPA",           "dl",
  "IGP12_INCCD12",        "INCC",          "dl",
  "BM12_IPCAEXP1212",     "IPCA_exp",      "d",
  "BM12_M1MN12",          "M1",            "dl",
  "BM12_CRLIN12",         "Credit",        "dl",
  "BM12_SPREAD12",        "Spread",        "d",
  "FUNCEX12_XVTOT12",     "Exports",       "dl",
  "FUNCEX12_MVTOT12",     "Imports",       "dl",
  "BM12_RESERVAS12",      "FXReserves",    "dl",
  "JPM366_EMBI366",       "EMBI",          "d",
  "GM366_IBVSP366",       "Ibovespa",      "dl",
  "BM12_ERREF12",         "RealFX",        "dl",
  "BM12_RNPSP12",         "PrimBalance",   "d",
  "BM12_DLSPN12",         "NetDebt_GDP",   "d"
)

# ============================================================
# 2. FUNÇÃO GENÉRICA DE DOWNLOAD — IPEADATA
# ============================================================

fetch_ipea <- function(code, name, start, end) {
  cat(sprintf("  [IPEA] %-30s %-16s ... ", code, name))
  tryCatch({
    df <- ipeadata(code)
    
    # REVISÃO: Verifica explicitamente se a série existe e tem a coluna 'date'
    if (nrow(df) == 0 || !"date" %in% names(df)) {
      cat("FAILED: Série descontinuada ou inexistente.\n")
      return(NULL)
    }
    
    if ("uname" %in% names(df)) {
      nacional <- df |>
        filter(uname %in% c("", "Brasil") | is.na(uname))
      if (nrow(nacional) > 0) df <- nacional
    }
    
    # REVISÃO: Uso do .data$date blinda contra confusões de escopo com base::date()
    df <- df |>
      filter(.data$date >= start, .data$date <= end) |>
      select(.data$date, value) |>
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
# ============================================================

fetch_unemployment <- function(start, end) {
  cat("  Baixando desemprego (emenda PME/BCB + PNAD-C/IPEA)...\n")
  
  # REVISÃO: O código 24369 é PNAD-C (só tem pós 2012). O código da PME antiga é 21774.
  cat("    [BCB  ] PME unemployment rate (code 21774) ... ")
  df_pme <- tryCatch({
    df <- gbcbd_get_series(
      id          = c(UNEMP_PME = 21774),
      first.date  = start,
      last.date   = as.Date("2011-12-01"),
      format.data = "wide",
      use.memoise = FALSE 
    ) |>
      rename(date = ref.date) |>
      select(date, DESEMPREGO = UNEMP_PME)
    cat(sprintf("%d obs. (PME up to Dec/2011)\n", nrow(df)))
    df
  }, error = function(e) {
    cat(sprintf("FAILED: %s\n", e$message)); NULL
  })
  
  cat("    [IPEA ] PNADC12_TDESOCMD12 (monthly seasonally adj.) ... ")
  df_pnadc <- tryCatch({
    raw <- ipeadata("PNADC12_TDESOCMD12")
    if ("uname" %in% names(raw)) {
      nac <- raw |> filter(uname %in% c("", "Brasil") | is.na(uname))
      if (nrow(nac) > 0) raw <- nac
    }
    df <- raw |>
      filter(.data$date >= as.Date("2012-01-01"), .data$date <= end) |>
      select(.data$date, value) |>
      rename(DESEMPREGO = value)
    cat(sprintf("%d obs. (PNAD-C from Jan/2012)\n", nrow(df)))
    df
  }, error = function(e) {
    cat(sprintf("FAILED: %s\n", e$message)); NULL
  })
  
  if (!is.null(df_pme) && !is.null(df_pnadc)) {
    cat("    Computing level offset for PME -> PNAD-C splice...\n")
    pme_overlap <- tryCatch(
      gbcbd_get_series(
        id          = c(PME_ext = 21774),
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
      filter(.data$date >= as.Date("2012-01-01"),
             .data$date <= as.Date("2015-12-01")) |>
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
  
  bind_rows(df_pme, df_pnadc) |> arrange(date)
}

# ============================================================
# 4. DOWNLOAD — TARGETS & PANEL
# ============================================================

cat("\n=== Downloading targets ===\n")
targets_list <- map2(
  targets_catalog$code,
  targets_catalog$name,
  ~fetch_ipea(.x, .y, start_date, end_date)
)
names(targets_list) <- targets_catalog$name
targets_list        <- compact(targets_list)
targets_list$DESEMPREGO <- fetch_unemployment(start_date, end_date)

cat("\n=== Downloading auxiliary panel ===\n")
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
# ============================================================

cat("\n=== Building wide panel ===\n")

all_series <- c(targets_list, panel_list)

df_wide <- all_series |>
  reduce(full_join, by = "date") |>
  arrange(date) |>
  filter(date >= start_date, date <= end_date)

collapse_to_monthly <- function(df, col_name) {
  df |>
    mutate(date = as.Date(format(date, "%Y-%m-01"))) |>
    group_by(date) |>
    summarise(!!col_name := mean(.data[[col_name]], na.rm = TRUE),
              .groups = "drop")
}

daily_series <- names(panel_list)[sapply(panel_list, function(x) nrow(x) > 400)]

if (length(daily_series) > 0) {
  cat(sprintf("  Collapsing %d daily series to monthly frequency: %s\n",
              length(daily_series),
              paste(daily_series, collapse = ", ")))
  
  for (nm in daily_series) {
    panel_list[[nm]] <- collapse_to_monthly(panel_list[[nm]], nm)
  }
  
  all_series <- c(targets_list, panel_list)
  df_wide <- all_series |>
    reduce(full_join, by = "date") |>
    arrange(date) |>
    filter(date >= start_date, date <= end_date)
}

cat(sprintf("  Raw dimensions: %d rows x %d columns\n", nrow(df_wide), ncol(df_wide)))

# ============================================================
# 7. APLICA TRANSFORMAÇÕES ESTACIONÁRIAS
# ============================================================

cat("\n=== Applying stationarity transformations ===\n")

full_catalog <- bind_rows(
  targets_catalog,
  tibble(code = NA_character_, name = "DESEMPREGO", transform = "d"),
  panel_catalog
)

# REVISÃO: Transforma valores negativos ou zero em NA antes do logaritmo
apply_transform <- function(x, trf) {
  switch(trf,
         "dl"  = {
           x[x <= 0] <- NA
           c(NA, diff(log(x)))
         },
         "d"   = c(NA, diff(x)),
         "l"   = {
           x[x <= 0] <- NA
           log(x)
         },
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

df_transf <- df_transf |> slice(-1)

# ============================================================
# 8 e 9. GRID DE COMBINAÇÕES E EXPORT
# ============================================================

targets_br <- list(
  V1 = "PIB", V2 = "IPCA", V3 = "SELIC", V4 = "CAMBIO", V5 = "DESEMPREGO"
)

all_options <- expand.grid(M = 1:4, V = 1:5, H = c(1, 2, 4))

save(df_wide,     file = file.path(run_folder, "df_wide.rda"))
save(df_transf,   file = file.path(run_folder, "df_transf.rda"))
save(targets_br,  file = file.path(run_folder, "targets_br.rda"))
save(all_options, file = file.path(run_folder, "all_options.rda"))

cat(sprintf("\n=== Files saved to %s/ ===\n", run_folder))
