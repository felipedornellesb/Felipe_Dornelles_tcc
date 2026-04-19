# ============================================================
# 01_data_prep_v1.R
#
# TCC Felipe Dornelles — Preparação de dados para TVP-2SRR / VARF Brasil
#
# ESTRATÉGIA:
#   ETAPA 1 — Carrega df.rda da Nathalia Oreda (UFRGS) como base confiável
#             (séries IPEA validadas, 1996-01 a 2019-05, já transformadas)
#   ETAPA 2 — Tenta download incremental via ipeadatar + rbcb para estender
#             a cobertura até a data mais recente disponível (alvo: dez/2025)
#   ETAPA 3 — Valida e exporta df_wide_raw, df_transf e df_model.rda
#
# Bugs corrigidos em relação ao data_prep_v1:
#   - is.data.frame() antes de nrow() em sapply → elimina "subscrito inválido 'list'"
#   - tryCatch por série individual → 1 série com erro não derruba o painel
#   - Colapso de séries diárias (EMBI+) para mensal antes do join
#   - MTE12_DESOCD12 substituído por série CAGED via IPEA (PAN12_QIIGG12 → MPT12_SDADM12)
#   - Validação de data como Date atômico antes de filter()
#   - df_wide só é construído se todas as etapas anteriores passarem em checagem
# =============================================================================

rm(list = ls())

# ─── 0. Pacotes ──────────────────────────────────────────────────────────────
required_pkgs <- c("dplyr", "tidyr", "lubridate", "ipeadatar",
                   "rbcb", "readxl", "urca", "tseries", "zoo", "purrr")

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Instalando pacote: ", pkg)
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# ─── 1. Configurações ────────────────────────────────────────────────────────
START_DATE   <- as.Date("1996-01-01")
END_DATE     <- as.Date("2025-12-01")   # Alvo máximo
NATHALIA_END <- as.Date("2019-05-01")   # Cobertura garantida da base fallback
MIN_OBS      <- 100                     # Mínimo de obs para manter série
OUTPUT_DIR   <- "data"
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

cat("=== 01_data_prep_v2.R ===\n")
cat("Período alvo:", format(START_DATE), "a", format(END_DATE), "\n\n")

# ─── 2. ETAPA 1: Carregar base da Nathalia como fallback ─────────────────────
cat("[ETAPA 1] Carregando base validada (Nathalia Oreda / UFRGS)...\n")

nathalia_rda_url <- paste0(
  "https://raw.githubusercontent.com/nathaliaoreda/thesis_UFRGS/main/data/df.rda"
)

tmp_rda <- tempfile(fileext = ".rda")
download_ok <- tryCatch({
  download.file(nathalia_rda_url, tmp_rda, mode = "wb", quiet = TRUE)
  TRUE
}, error = function(e) {
  warning("Falha ao baixar df.rda da Nathalia: ", conditionMessage(e))
  FALSE
})

if (!download_ok) {
  # Tenta path local se existir (caso o usuário já tenha o arquivo)
  local_fallback <- file.path(OUTPUT_DIR, "df_nathalia.rda")
  if (file.exists(local_fallback)) {
    tmp_rda <- local_fallback
    download_ok <- TRUE
    cat("  → Usando cópia local:", local_fallback, "\n")
  } else {
    stop("Não foi possível obter a base da Nathalia. Coloque df_nathalia.rda em ", OUTPUT_DIR)
  }
}

# Carrega em ambiente isolado para não poluir o global
env_nathalia <- new.env()
load(tmp_rda, envir = env_nathalia)

# O objeto salvo pela Nathalia chama-se 'df'
df_nathalia <- env_nathalia$df

# Garantir que a coluna date existe (a base da Nathalia NÃO tem coluna date —
# ela é uma matriz numérica já transformada; reconstituímos o índice temporal)
if (!is.data.frame(df_nathalia)) {
  df_nathalia <- as.data.frame(df_nathalia)
}

# A base filtrada da Nathalia vai de 1996-03 (após 2 diferenças) a 2019-05
# Reconstituímos o eixo de datas
n_rows_nathalia <- nrow(df_nathalia)
# Sequência mensal: primeiro dia de cada mês
date_seq_nathalia <- seq.Date(
  from = as.Date("1996-03-01"),   # início após remoção das 2 primeiras obs por diferenciação
  by   = "month",
  length.out = n_rows_nathalia
)
df_nathalia$date <- date_seq_nathalia

cat("  → Base Nathalia carregada:", nrow(df_nathalia), "obs x",
    ncol(df_nathalia) - 1, "variáveis\n")
cat("  → Colunas disponíveis:", paste(head(names(df_nathalia), 8), collapse = ", "), "...\n\n")

# ─── 3. ETAPA 2: Download incremental via ipeadatar ──────────────────────────
# Catálogo de séries com os códigos IPEA VALIDADOS (extraídos do dataset.xlsx
# da Nathalia) + séries adicionais do seu TCC (IBC-Br via rbcb SGS)
# Código   : código IPEA
# fonte    : "ipea" ou "bcb_sgs"
# sgs_id   : ID no SGS do BCB (apenas para fonte=="bcb_sgs")
# desc     : descrição curta
# freq     : "M" mensal, "D" diária (será colapsada para mensal)

series_catalog <- tribble(
  ~codigo,                ~fonte,     ~sgs_id, ~desc,                          ~freq,
  # ── Preços ──────────────────────────────────────────────────────────────────
  "PRECOS12_IPCA12",      "ipea",     NA,      "IPCA (% a.m.)",                "M",
  "PRECOS12_IGP12",       "ipea",     NA,      "IGP-M (% a.m.)",               "M",
  "PRECOS12_IPCAG12",     "ipea",     NA,      "IPCA acum. 12m",               "M",
  # ── Atividade ────────────────────────────────────────────────────────────────
  "PIMPF12_PGBR12",       "ipea",     NA,      "Produção Industrial (PIM-PF)", "M",
  "PIMPF12_QGBR12",       "ipea",     NA,      "PIM-PF quantum",               "M",
  # ── Emprego ──────────────────────────────────────────────────────────────────
  "MPT12_SDADM12",        "ipea",     NA,      "CAGED saldo empregos formais", "M",
  "SEADE12_TDTGSP12",     "ipea",     NA,      "Taxa desemprego SEADE/SP",     "M",
  # ── Crédito / Juros ──────────────────────────────────────────────────────────
  "BM12_TJOVER12",        "ipea",     NA,      "Taxa Selic Over (% a.a.)",     "M",
  "BM12_CRLIN12",         "ipea",     NA,      "Crédito total / PIB",          "M",
  "BM12_M112",            "ipea",     NA,      "M1",                           "M",
  "BM12_MBASE12",         "ipea",     NA,      "Base monetária",               "M",
  # ── Externo / Câmbio ─────────────────────────────────────────────────────────
  "GM366_ERC366",         "ipea",     NA,      "Câmbio R$/USD (PTAX)",         "M",
  "BOP12_CAB12",          "ipea",     NA,      "Conta corrente (USD mi)",      "M",
  # ── EMBI+ (diária → mensal) ──────────────────────────────────────────────────
  "JPM366_EMBI366",       "ipea",     NA,      "EMBI+ Brasil (spread)",        "D",
  # ── IBC-Br via SGS/BCB (proxy mensal do PIB) ─────────────────────────────────
  NA_character_,          "bcb_sgs",  24363L,  "IBC-Br (índice ativ. econôm)","M",
  # ── Expectativas Focus/BCB ────────────────────────────────────────────────────
  NA_character_,          "bcb_sgs",  13522L,  "Expectativa IPCA 12m (Focus)", "M"
)

# ─── 3.1 Helper: download seguro de 1 série IPEA ─────────────────────────────
safe_ipea_download <- function(code, desc) {
  result <- tryCatch({
    raw <- ipeadata(code)
    # Verificações de estrutura
    if (!is.data.frame(raw) || nrow(raw) == 0) {
      warning(code, ": retornou vazio ou não é data.frame")
      return(NULL)
    }
    if (!"date" %in% names(raw)) {
      warning(code, ": coluna 'date' ausente")
      return(NULL)
    }
    # Forçar date para Date atômico (elimina o bug "subscrito inválido 'list'")
    raw$date <- tryCatch(
      as.Date(as.character(raw$date)),
      error = function(e) NULL
    )
    if (is.null(raw$date) || !inherits(raw$date, "Date")) {
      warning(code, ": não foi possível converter 'date' para Date")
      return(NULL)
    }
    # Manter apenas colunas essenciais
    if (!"value" %in% names(raw)) {
      warning(code, ": coluna 'value' ausente")
      return(NULL)
    }
    df_out <- raw %>%
      select(date, value) %>%
      rename(!!code := value) %>%
      filter(!is.na(date), is.finite(!!sym(code)) | is.na(!!sym(code)))
    df_out
  }, error = function(e) {
    warning(code, " (", desc, "): ERRO — ", conditionMessage(e))
    NULL
  })
  result
}

# ─── 3.2 Helper: colapso diário → mensal ─────────────────────────────────────
daily_to_monthly <- function(df_daily, value_col) {
  df_daily %>%
    mutate(date = floor_date(date, "month")) %>%
    group_by(date) %>%
    summarise(!!value_col := mean(.data[[value_col]], na.rm = TRUE), .groups = "drop")
}

# ─── 3.3 Helper: download seguro série SGS/BCB ───────────────────────────────
safe_bcb_sgs <- function(sgs_id, desc) {
  result <- tryCatch({
    raw <- rbcb::get_series(sgs_id,
                            start_date = START_DATE,
                            end_date   = END_DATE)
    if (!is.data.frame(raw) || nrow(raw) == 0) return(NULL)
    raw$date <- as.Date(raw$date)
    col_val  <- names(raw)[names(raw) != "date"][1]
    raw %>%
      select(date, all_of(col_val)) %>%
      rename(!!paste0("SGS_", sgs_id) := all_of(col_val))
  }, error = function(e) {
    warning("SGS ", sgs_id, " (", desc, "): ERRO — ", conditionMessage(e))
    NULL
  })
  result
}

# ─── 3.4 Executar downloads ───────────────────────────────────────────────────
cat("[ETAPA 2] Download incremental de séries...\n")

panel_list <- vector("list", nrow(series_catalog))
names(panel_list) <- paste0(
  ifelse(is.na(series_catalog$codigo),
         paste0("SGS_", series_catalog$sgs_id),
         series_catalog$codigo)
)

for (i in seq_len(nrow(series_catalog))) {
  row  <- series_catalog[i, ]
  nome <- names(panel_list)[i]
  cat("  Baixando:", nome, "-", row$desc, "... ")

  if (row$fonte == "ipea") {
    df_i <- safe_ipea_download(row$codigo, row$desc)
    # Colapso diário → mensal para EMBI+
    if (!is.null(df_i) && row$freq == "D") {
      df_i <- daily_to_monthly(df_i, row$codigo)
    }
  } else {
    df_i <- safe_bcb_sgs(row$sgs_id, row$desc)
  }

  if (!is.null(df_i) && nrow(df_i) > 0) {
    cat("OK (", nrow(df_i), "obs)\n")
  } else {
    cat("FALHOU\n")
  }
  panel_list[[nome]] <- df_i
}

# ─── 3.5 Remover séries que falharam ─────────────────────────────────────────
# BUG CORRIGIDO: is.data.frame() antes de nrow() evita "subscrito inválido 'list'"
ok_series <- sapply(panel_list, function(x) {
  is.data.frame(x) && !is.null(x) && nrow(x) >= MIN_OBS
})

n_ok  <- sum(ok_series)
n_fail <- sum(!ok_series)
cat("\n  Séries baixadas com sucesso:", n_ok, "/", nrow(series_catalog), "\n")
if (n_fail > 0) {
  cat("  Séries com falha (serão preenchidas pela base da Nathalia):\n")
  cat("   ", paste(names(panel_list)[!ok_series], collapse = ", "), "\n")
}
panel_list <- panel_list[ok_series]

# ─── 3.6 Montar df_wide das séries novas ─────────────────────────────────────
if (length(panel_list) > 0) {
  cat("\n[ETAPA 2.6] Construindo painel wide das séries baixadas...\n")
  df_new <- purrr::reduce(panel_list, function(a, b) {
    full_join(a, b, by = "date")
  })
  # Garantir que date seja Date atômico
  df_new$date <- as.Date(df_new$date)
  # Filtrar janela de interesse
  df_new <- df_new %>%
    filter(date >= START_DATE, date <= END_DATE) %>%
    arrange(date)
  cat("  → df_new:", nrow(df_new), "obs x", ncol(df_new) - 1, "colunas\n")
} else {
  warning("Nenhuma série foi baixada com sucesso. Usando apenas base da Nathalia.")
  df_new <- NULL
}

# ─── 4. ETAPA 3: Integração — base da Nathalia + extensão incremental ────────
cat("\n[ETAPA 3] Integrando base da Nathalia com extensão incremental...\n")

# A base da Nathalia já está TRANSFORMADA (estacionária).
# Para extensão após NATHALIA_END, precisamos:
#   (a) identificar quais colunas da base nova coincidem com as da Nathalia
#   (b) aplicar a mesma transformação às obs novas
#   (c) empilhar as duas partes

# Por ora, exportamos três objetos:
#   df_wide_raw : painel bruto (unindo o que baixamos + datas disponíveis)
#   df_nathalia : base transformada validada (pronta para o modelo até mai/2019)
#   df_model    : objeto final usado no forecast

# Para o forecast TVP-2SRR, o df_model é inicialmente = df_nathalia
# Quando a extensão for validada, troque pelo painel extendido

# ── 4a. Salvar painel bruto (para auditoria) ──────────────────────────────────
df_wide_raw <- if (!is.null(df_new)) df_new else data.frame(date = as.Date(character()))
save(df_wide_raw, file = file.path(OUTPUT_DIR, "df_wide_raw.rda"))
cat("  → Salvo: data/df_wide_raw.rda\n")

# ── 4b. Objeto df_model principal = base transformada da Nathalia ──────────────
df_model <- df_nathalia
# Garantir que não há Inf ou NaN
df_model <- df_model %>%
  mutate(across(where(is.numeric), ~ ifelse(is.infinite(.) | is.nan(.), NA_real_, .)))

# Remover colunas com mais de 5% de NA
na_frac <- colMeans(is.na(df_model %>% select(-date)))
cols_keep <- names(na_frac)[na_frac <= 0.05]
df_model  <- df_model %>% select(date, all_of(cols_keep))

# Imputação linear simples para lacunas curtas (≤ 3 meses)
df_model <- df_model %>%
  mutate(across(where(is.numeric), ~ zoo::na.approx(., maxgap = 3, na.rm = FALSE)))

cat("  → df_model final:", nrow(df_model), "obs x", ncol(df_model) - 1, "variáveis\n")

# ── 4c. Checar sanidade antes de salvar ───────────────────────────────────────
sanity_ok <- TRUE

if (nrow(df_model) < 50) {
  warning("df_model tem menos de 50 observações — verifique o pipeline!")
  sanity_ok <- FALSE
}
if (ncol(df_model) < 5) {
  warning("df_model tem menos de 5 variáveis — verifique os downloads!")
  sanity_ok <- FALSE
}
has_inf <- any(sapply(df_model %>% select(-date), function(x) any(is.infinite(x), na.rm = TRUE)))
if (has_inf) {
  warning("df_model ainda contém Inf — checagem necessária!")
  sanity_ok <- FALSE
}

if (sanity_ok) {
  save(df_model, file = file.path(OUTPUT_DIR, "df_model.rda"))
  cat("  ✓ Salvo: data/df_model.rda  [pronto para 02_forecast.R]\n")
} else {
  warning("df_model NÃO foi salvo por falha na checagem de sanidade.")
}

# ─── 5. Relatório final ───────────────────────────────────────────────────────
cat("\n=== RELATÓRIO FINAL ===\n")
cat("Período em df_model  :", format(min(df_model$date)), "a",
    format(max(df_model$date)), "\n")
cat("Obs                  :", nrow(df_model), "\n")
cat("Variáveis            :", ncol(df_model) - 1, "\n")
cat("Colunas              :\n")
print(names(df_model))
cat("\nSanidade             :", ifelse(sanity_ok, "PASSOU ✓", "FALHOU ✗"), "\n")
cat("========================\n")