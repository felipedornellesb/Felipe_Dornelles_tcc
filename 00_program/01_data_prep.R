# ============================================================
# 01_data_prep.R
#
# Baixa dados via ipeadatar, aplica transformações FRED-MD,
# exporta dataset.xlsx e os .rda para 10_data/data_MM_DD_YYYY/
#
# Outputs em 10_data/data_MM_DD_YYYY/:
#   dataset.xlsx        <- planilha bruta (formato FRED-MD)
#   df_model.rda        <- série estacionária + coluna date
#   df_targets.rda      <- somente as colunas-alvo (estacionárias)
#   df_panel_pca.rda    <- painel de covariadas (estacionárias)
#   targets_br.rda      <- lista nomeada V1..V4 com nomes das variáveis
#   all_options.rda     <- grade M x V x H (expand.grid)
#   fullcatalog.rda     <- catálogo completo de séries
# ============================================================

rm(list = ls())

# ============================================================
# 0. PATHS
# ============================================================

# Detecta a raiz do projeto (assume que o script está em 00_program/)
wd <- normalizePath(file.path(dirname(sys.frame(1)$ofile), ".."),
                    mustWork = FALSE)
if (!dir.exists(wd)) wd <- getwd()   # fallback interativo
setwd(wd)

paths <- list(
  data    = "10_data",
  tools   = "20_tools",
  output  = "30_output",
  results = "40_results"
)

run_date   <- format(Sys.Date(), "%m_%d_%Y")
data_folder <- file.path(paths$data, paste0("data_", run_date))
for (p in c(data_folder, paths$output, paths$results))
  if (!dir.exists(p)) dir.create(p, recursive = TRUE)

cat(sprintf("Pasta de dados: %s\n", data_folder))

# ============================================================
# 1. PACOTES
# ============================================================

myPKGs <- c("dplyr", "ipeadatar", "openxlsx", "lubridate",
            "tidyr", "purrr", "zoo")
need   <- myPKGs[!myPKGs %in% names(installed.packages()[, "Package"])]
if (length(need) > 0)
  install.packages(need, repos = "http://cran.us.r-project.org")
invisible(lapply(myPKGs, library, character.only = TRUE))

# ============================================================
# 2. CATÁLOGOS
# ============================================================

start_date <- as.Date("1996-01-01")
end_date   <- as.Date("2025-12-01")

# Tcodes (FRED-MD):
#  1 = nível, 2 = Δ, 3 = ΔΔ, 4 = log, 5 = Δlog, 6 = ΔΔlog, 7 = Δ(log-1)
targets_catalog <- tribble(
  ~code,              ~name,    ~tcode,
  "BM12_PIB12",       "PIB",    5,
  "PRECOS12_IPCA12",  "IPCA",   5,
  "BM12_TJOVER12",    "SELIC",  2,
  "BM12_ERC12",       "CAMBIO", 5
)

panel_catalog <- tribble(
  ~code,                   ~name,           ~tcode,
  "PIMPF12_QTIG12",        "PIMgeral",       5,
  "PIMPF12_QTBK12",        "PIMbkap",        5,
  "PIMPF12_QTBCD12",       "PIMbdur",        5,
  "PIMPF12_QTBI12",        "PIMbint",        5,
  "BM12_CEEI12",           "EnergyConsump",  5,
  "PMC12_VVTOT12",         "RetailSales",    5,
  "MTE12_SALDON12",        "CAGEDnet",       1,
  "PMEN12_RRME12",         "RealWage",       5,
  "IGP12_IGPDI_G12",       "IGPDI",          5,
  "IGP12_IPADI_G12",       "IPA",            5,
  "IGP12_INCC_D12",        "INCC",           5,
  "BM12_IPCAEXP1212",      "IPCAexp",        2,
  "BM12_M1MN12",           "M1",             5,
  "BM12_CRLIN12",          "Credit",         5,
  "BM12_SPREAD12",         "Spread",         2,
  "FUNCEX12_XVTOT12",      "Exports",        5,
  "FUNCEX12_MVTOT12",      "Imports",        5,
  "BM12_RESERVAS12",       "FXReserves",     5,
  "JPM366_EMBI366",        "EMBI",           2,
  "GM366_IBVSP366",        "Ibovespa",       5,
  "BM12_ERREF12",          "RealFX",         5,
  "BM12_RNPSP12",          "PrimBalance",    2,
  "BM12_DLSPN12",          "NetDebtGDP",     2
)

full_catalog <- bind_rows(targets_catalog, panel_catalog)

# ============================================================
# 3. DOWNLOAD VIA IPEADATAR
# ============================================================

fetch_ipea <- function(code, name, start, end) {
  cat(sprintf("  %-30s %-16s ", code, name))
  tryCatch({
    raw <- ipeadatar::ipeadata(code)
    if (nrow(raw) == 0 || !"date" %in% names(raw)) {
      cat("VAZIA\n"); return(NULL)
    }
    # Filtra apenas Brasil quando coluna uname existe
    if ("uname" %in% names(raw)) {
      bra <- dplyr::filter(raw, uname %in% c("", "Brasil") | is.na(uname))
      if (nrow(bra) > 0) raw <- bra
    }
    out <- raw |>
      dplyr::filter(date >= start, date <= end) |>
      dplyr::select(date, value) |>
      dplyr::rename(!!name := value)
    cat(sprintf("%d obs\n", nrow(out)))
    out
  }, error = function(e) {
    cat(sprintf("ERRO: %s\n", conditionMessage(e)))
    NULL
  })
}

cat("\n--- Baixando TARGETS ---\n")
tgt_list <- Map(fetch_ipea,
                targets_catalog$code,
                targets_catalog$name,
                MoreArgs = list(start = start_date, end = end_date))
tgt_list <- setNames(Filter(Negate(is.null), tgt_list),
                     targets_catalog$name[!sapply(tgt_list, is.null)])

cat("\n--- Baixando PAINEL ---\n")
pan_list <- Map(fetch_ipea,
                panel_catalog$code,
                panel_catalog$name,
                MoreArgs = list(start = start_date, end = end_date))
pan_list <- setNames(Filter(Negate(is.null), pan_list),
                     panel_catalog$name[!sapply(pan_list, is.null)])

# ============================================================
# 4. COLAPSO PARA FREQUÊNCIA MENSAL (séries diárias)
# ============================================================

collapse_monthly <- function(df, nm) {
  df |>
    dplyr::mutate(date = as.Date(format(date, "%Y-%m-01"))) |>
    dplyr::group_by(date) |>
    dplyr::summarise(!!nm := mean(.data[[nm]], na.rm = TRUE), .groups = "drop")
}

is_daily <- function(df) nrow(df) > 400  # proxy: >400 obs → diária

all_series <- c(tgt_list, pan_list)
for (nm in names(all_series)) {
  if (is_daily(all_series[[nm]])) {
    cat(sprintf("  Colapsando mensal: %s\n", nm))
    all_series[[nm]] <- collapse_monthly(all_series[[nm]], nm)
  }
}

# ============================================================
# 5. MONTAR PAINEL WIDE + EXPORTAR dataset.xlsx
# ============================================================

df_wide <- all_series |>
  purrr::reduce(dplyr::full_join, by = "date") |>
  dplyr::arrange(date) |>
  dplyr::filter(date >= start_date, date <= end_date)

cat(sprintf("\nDimensão bruta: %d obs × %d colunas\n",
            nrow(df_wide), ncol(df_wide)))

# Linha de tcodes (formato FRED-MD)
tcode_row <- df_wide[1, ]
tcode_row[1, ] <- NA
for (nm in names(df_wide)[-1]) {
  tc <- full_catalog$tcode[full_catalog$name == nm]
  tcode_row[1, nm] <- if (length(tc) > 0) tc else NA
}
tcode_row$date <- NA

df_xlsx             <- dplyr::bind_rows(tcode_row, df_wide)
df_xlsx$date        <- as.character(df_xlsx$date)
df_xlsx$date[1]     <- "Transform"
openxlsx::write.xlsx(df_xlsx,
                     file = file.path(data_folder, "dataset.xlsx"))
cat("dataset.xlsx exportado\n")

# ============================================================
# 6. TRANSFORMAÇÕES FRED-MD → ESTACIONARIEDADE
# ============================================================

apply_tcode <- function(x, tcode) {
  x <- as.numeric(x)
  switch(as.character(tcode),
    "1" = x,
    "2" = c(NA, diff(x)),
    "3" = c(NA, NA, diff(diff(x))),
    "4" = { x[x <= 0] <- NA; log(x) },
    "5" = { x[x <= 0] <- NA; c(NA, diff(log(x))) },
    "6" = { x[x <= 0] <- NA; c(NA, NA, diff(diff(log(x)))) },
    "7" = c(NA, diff(x / dplyr::lag(x) - 1)),
    x  # fallback: nível
  )
}

df_stat <- df_wide
for (nm in names(df_stat)[-1]) {
  tc <- full_catalog$tcode[full_catalog$name == nm]
  if (length(tc) == 1) df_stat[[nm]] <- apply_tcode(df_stat[[nm]], tc)
}

# Remove linhas com NA em TODOS os alvos (burn-in das diferenças)
target_names <- targets_catalog$name[
  targets_catalog$name %in% names(df_stat)
]
df_stat <- df_stat |>
  dplyr::filter(dplyr::if_any(dplyr::all_of(target_names),
                              ~ !is.na(.x)))

# Alinha todos ao mesmo período
df_stat <- df_stat |> dplyr::filter(date >= start_date, date <= end_date)

cat(sprintf("Dimensão estacionária: %d obs × %d colunas\n",
            nrow(df_stat), ncol(df_stat)))
cat(sprintf("Período: %s  a  %s\n",
            min(df_stat$date), max(df_stat$date)))

# ============================================================
# 7. SEPARAR df_model / df_targets / df_panel_pca
# ============================================================

panel_names <- panel_catalog$name[
  panel_catalog$name %in% names(df_stat)
]

df_model      <- df_stat[, c("date", target_names, panel_names)]
df_targets    <- df_stat[, c("date", target_names)]
df_panel_pca  <- df_stat[, c("date", panel_names)]

# ============================================================
# 8. METADADOS
# ============================================================

targets_br <- as.list(target_names)
names(targets_br) <- paste0("V", seq_along(targets_br))

all_options <- expand.grid(
  M = 1:4,
  V = seq_along(targets_br),
  H = c(1L, 3L, 6L, 12L),
  stringsAsFactors = FALSE
)

# ============================================================
# 9. SALVAR .rda
# ============================================================

save(df_model,     file = file.path(data_folder, "df_model.rda"))
save(df_targets,   file = file.path(data_folder, "df_targets.rda"))
save(df_panel_pca, file = file.path(data_folder, "df_panel_pca.rda"))
save(targets_br,   file = file.path(data_folder, "targets_br.rda"))
save(all_options,  file = file.path(data_folder, "all_options.rda"))
save(full_catalog, file = file.path(data_folder, "fullcatalog.rda"))

cat("\n=== data_prep.R CONCLUÍDO ===\n")
cat(sprintf("Arquivos salvos em: %s\n", data_folder))
cat("  df_model.rda\n")
cat("  df_targets.rda\n")
cat("  df_panel_pca.rda\n")
cat("  targets_br.rda\n")
cat("  all_options.rda\n")
cat("  fullcatalog.rda\n")
cat("  dataset.xlsx\n")
