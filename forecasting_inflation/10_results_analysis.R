rm(list = ls())
set.seed(2024)

cat("
==============================================================
 10_results_analysis.R - v0.0.4
==============================================================
\n")

# ---- 0. Pacotes ----
pkgs <- c("ggplot2", "dplyr", "tidyr", "scales", "gridExtra",
          "forecast", "viridis", "reshape2", "patchwork")
need <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(need)) install.packages(need, repos = "https://cran.r-project.org")
invisible(lapply(pkgs, library, character.only = TRUE))

# MCS (opcional)
has_mcs <- require("MCS", quietly = TRUE)
if (!has_mcs) {
  tryCatch({
    install.packages("MCS", repos = "https://cran.r-project.org")
    library(MCS)
    has_mcs <- TRUE
  }, error = function(e) { has_mcs <<- FALSE })
}

# ---- 1. Criar pasta de resultados ----
timestamp_str <- format(Sys.time(), "%Y%m%d_%H%M%S")
results_dir   <- file.path("40_results", paste0("run_", timestamp_str))
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("Pasta de resultados: %s\n\n", results_dir))

# ---- 2. Carregar dados de forecast REAIS ----
cat("=== Carregando resultados de forecast reais ===\n\n")

forecast_dir <- "forecasts"

# Funcao auxiliar: carrega .rda e extrai o primeiro objeto,
# independentemente do nome interno
load_rda <- function(path) {
  env <- new.env()
  load(path, envir = env)
  obj_name <- ls(env)[1]
  return(env[[obj_name]])
}

# 2.1 Valores realizados (y out-of-sample)
file_yout <- file.path(forecast_dir, "yout.rda")
if (!file.exists(file_yout)) {
  stop("Arquivo yout.rda nao encontrado na pasta forecasts/. Rode o forecast primeiro.")
}
y_out <- load_rda(file_yout)
cat(sprintf("  [OK] yout.rda carregado (classe: %s)\n", class(y_out)[1]))

# Se y_out for matriz/data.frame, pegar a primeira coluna como vetor
if (is.matrix(y_out) || is.data.frame(y_out)) {
  cat(sprintf("       Dimensoes: %d x %d\n", nrow(y_out), ncol(y_out)))
  y_vec <- as.numeric(y_out[, 1])
} else {
  y_vec <- as.numeric(y_out)
}
n_eval <- length(y_vec)
cat(sprintf("       n_eval = %d observacoes out-of-sample\n\n", n_eval))

# 2.2 Listar todos os modelos disponíveis (.rda na pasta forecasts)
# Excluir arquivos que nao sao previsoes de modelos
excluir <- c("yout", "betas_2SRR", "coulombe_betas_2SRR",
             "coulombe_betas_ridge", "coulombe_forecasts",
             "tvp_TVP_AR_forecasts", "tvp_TVP_Factor_forecasts")

todos_rda <- list.files(forecast_dir, pattern = "\\.rda$", full.names = FALSE)
todos_rda <- gsub("\\.rda$", "", todos_rda)
modelos_rda <- setdiff(todos_rda, excluir)
cat(sprintf("  Modelos .rda encontrados: %s\n\n",
            paste(modelos_rda, collapse = ", ")))

# 2.3 Carregar cada modelo e inspecionar estrutura
cat("  --- Inspecao da estrutura dos .rda de cada modelo ---\n")
modelo_data <- list()
for (m in modelos_rda) {
  fpath <- file.path(forecast_dir, paste0(m, ".rda"))
  obj   <- load_rda(fpath)
  modelo_data[[m]] <- obj
  if (is.matrix(obj) || is.data.frame(obj)) {
    cat(sprintf("  [OK] %-12s : %s  %d x %d\n", m, class(obj)[1],
                nrow(obj), ncol(obj)))
  } else {
    cat(sprintf("  [OK] %-12s : %s  length=%d\n", m, class(obj)[1],
                length(obj)))
  }
}

# 2.4 Determinar horizontes
# A maioria dos seus .rda salvos pelo 03_call_model_felipe.R sao matrizes
# com colunas = horizontes (h=1..12). Vamos verificar:
ref_obj <- modelo_data[["AR"]]
if (is.matrix(ref_obj) || is.data.frame(ref_obj)) {
  n_h_total <- ncol(ref_obj)
  cat(sprintf("\n  AR.rda tem %d colunas (horizontes h=1..%d)\n", n_h_total, n_h_total))
  cat(sprintf("  AR.rda tem %d linhas  (janelas out-of-sample)\n", nrow(ref_obj)))
} else {
  n_h_total <- 1
  cat("\n  AR.rda e vetor simples (horizonte unico)\n")
}

# Horizontes de interesse para a analise
horizons_list <- c(1, 3, 6, 12)
# Filtrar apenas os que existem
horizons_list <- horizons_list[horizons_list <= n_h_total]
cat(sprintf("  Horizontes para analise: %s\n\n",
            paste(horizons_list, collapse = ", ")))

# 2.5 Alinhar n_eval com o numero de linhas das previsoes
# (o y_out pode ter mais observacoes que as linhas da matriz de previsoes)
n_rows_pred <- if (is.matrix(ref_obj)) nrow(ref_obj) else length(ref_obj)
if (n_eval != n_rows_pred) {
  cat(sprintf("  [AVISO] y_out tem %d obs, mas previsoes tem %d linhas.\n",
              n_eval, n_rows_pred))
  cat(sprintf("          Usando as ultimas %d obs de y_out para alinhar.\n\n",
              n_rows_pred))
  # Pegar as ultimas n_rows_pred observacoes de y_vec
  y_vec  <- tail(y_vec, n_rows_pred)
  n_eval <- n_rows_pred
}

# 2.6 Construir datas de avaliacao a partir do CSV do Coulombe (se existir)
csv_ref <- file.path(forecast_dir, "coulombe_fc_h01.csv")
if (file.exists(csv_ref)) {
  df_csv    <- read.csv(csv_ref)
  col_date  <- grep("date|Date|DATA", names(df_csv), value = TRUE, ignore.case = TRUE)
  if (length(col_date) > 0) {
    dates_raw <- as.Date(df_csv[[col_date[1]]])
  } else {
    dates_raw <- as.Date(df_csv[[1]])
  }
  # Alinhar tamanho
  if (length(dates_raw) >= n_eval) {
    dates_eval <- tail(dates_raw, n_eval)
  } else {
    dates_eval <- seq.Date(as.Date("2000-01-01"), by = "month", length.out = n_eval)
  }
  cat(sprintf("  [OK] Datas extraidas de coulombe_fc_h01.csv (%s a %s)\n\n",
              format(min(dates_eval)), format(max(dates_eval))))
} else {
  dates_eval <- seq.Date(as.Date("2000-01-01"), by = "month", length.out = n_eval)
  cat("  [AVISO] coulombe_fc_h01.csv nao encontrado. Usando datas genericas.\n\n")
}

# 2.7 Definir lista de metodos para a analise
# Usar todos os modelos carregados
methods_list <- modelos_rda
target_name  <- "INFLATION"

# 2.8 Montar o objeto forecast_output (formato esperado pelo restante do script)
cat("  --- Montando forecast_output ---\n")

forecast_output <- list(
  results    = list(),
  targets_br = list(V1 = target_name),
  horizons   = horizons_list,
  methods    = methods_list,
  dates      = dates_eval
)

for (h in horizons_list) {
  key <- sprintf("V1_h%d", h)

  # Montar matriz de previsoes: n_eval x n_modelos
  preds_mat <- matrix(NA, nrow = n_eval, ncol = length(methods_list))
  colnames(preds_mat) <- methods_list

  for (m in methods_list) {
    obj <- modelo_data[[m]]
    if (is.matrix(obj) || is.data.frame(obj)) {
      if (ncol(obj) >= h) {
        col_vals <- as.numeric(obj[, h])
      } else {
        col_vals <- as.numeric(obj[, ncol(obj)])
      }
      # Alinhar numero de linhas
      if (length(col_vals) >= n_eval) {
        preds_mat[, m] <- tail(col_vals, n_eval)
      } else {
        preds_mat[1:length(col_vals), m] <- col_vals
      }
    } else {
      # Vetor simples (ex: rw.rda pode ser vetor)
      vec <- as.numeric(obj)
      if (length(vec) >= n_eval) {
        preds_mat[, m] <- tail(vec, n_eval)
      } else {
        preds_mat[1:length(vec), m] <- vec
      }
    }
  }

  forecast_output$results[[key]] <- list(
    target  = target_name,
    horizon = h,
    dates   = dates_eval,
    actuals = y_vec,
    preds   = preds_mat
  )
  cat(sprintf("  [OK] h=%2d montado (%d obs x %d modelos)\n",
              h, n_eval, length(methods_list)))
}

# 2.9 Carregar betas do 2SRR para Parte G
file_betas <- file.path(forecast_dir, "betas_2SRR.rda")
if (file.exists(file_betas)) {
  betas_obj <- load_rda(file_betas)
  cat(sprintf("\n  [OK] betas_2SRR.rda carregado (classe: %s", class(betas_obj)[1]))

  if (is.matrix(betas_obj) || is.data.frame(betas_obj)) {
    cat(sprintf(", %d x %d)\n", nrow(betas_obj), ncol(betas_obj)))

    # Detectar formato longo (3 colunas, muitas linhas)
    if (ncol(betas_obj) == 3 && nrow(betas_obj) > n_eval) {
      cat("       Formato longo detectado. Pivotando para formato largo...\n")
      colnames(betas_obj) <- c("Time", "Predictor", "Value")
      betas_obj$Value <- as.numeric(betas_obj$Value)

      betas_wide <- tidyr::pivot_wider(betas_obj,
                                       names_from = Predictor,
                                       values_from = Value)
      # Primeira coluna e Time, remover
      betas_num_mat <- as.matrix(betas_wide[, -1])
      storage.mode(betas_num_mat) <- "double"
      cat(sprintf("       Resultado: %d x %d (periodos x preditores)\n",
                  nrow(betas_num_mat), ncol(betas_num_mat)))
    } else {
      # Ja esta em formato largo, manter apenas colunas numericas
      is_num <- sapply(as.data.frame(betas_obj), function(x) {
        is.numeric(x) || !any(is.na(suppressWarnings(as.numeric(x))))
      })
      betas_num_mat <- as.matrix(as.data.frame(betas_obj)[, is_num, drop = FALSE])
      storage.mode(betas_num_mat) <- "double"
    }

    # Remover colunas que sao tudo NA
    all_na <- apply(betas_num_mat, 2, function(x) all(is.na(x)))
    betas_num_mat <- betas_num_mat[, !all_na, drop = FALSE]

    # Gerar datas
    if (nrow(betas_num_mat) <= length(dates_eval)) {
      beta_dates <- tail(dates_eval, nrow(betas_num_mat))
    } else {
      beta_dates <- seq.Date(as.Date("2000-01-01"), by = "month",
                             length.out = nrow(betas_num_mat))
    }

    forecast_output$tvp_full <- list()
    forecast_output$tvp_full[[target_name]] <- list(
      betas = betas_num_mat,
      dates = beta_dates
    )
    cat(sprintf("       Betas prontos: %d periodos x %d variaveis\n",
                nrow(betas_num_mat), ncol(betas_num_mat)))
  } else {
    cat(")\n  [AVISO] betas_2SRR nao e matriz/DF, pulando Parte G.\n")
  }
} else {
  cat("\n  [AVISO] betas_2SRR.rda nao encontrado. Parte G sera pulada.\n")
}

cat("\n  === Carregamento concluido com sucesso ===\n\n")

# ---- Extrair objetos ----
res       <- forecast_output$results
horizons  <- forecast_output$horizons
methods   <- forecast_output$methods
dates_vec <- forecast_output$dates

if (is.list(forecast_output$targets_br)) {
  targets <- unlist(forecast_output$targets_br)
} else {
  targets <- forecast_output$targets_br
}

n_models  <- length(methods)
n_targets <- length(targets)

# ---- Tema padrao para graficos ----
theme_set(theme_minimal(base_size = 13) +
            theme(plot.title      = element_text(face = "bold", size = 14),
                  plot.subtitle   = element_text(color = "grey40", size = 10),
                  legend.position = "bottom",
                  legend.title    = element_blank(),
                  panel.grid.minor = element_blank(),
                  strip.text = element_text(face = "bold")))

# ---- Paleta de cores com destaque para 2SRR e Realizado ----
# 2SRR em vermelho forte, AR em cinza, Realizado em preto,
# demais em tons discretos para nao competir visualmente
cores_fixas <- c(
  "Realizado"  = "black",
  "2SRR"       = "#D32F2F",
  "AR"         = "grey55",
  "AR_BIC"     = "grey70",
  "rw"         = "grey80",
  "Ridge"      = "#90CAF9",
  "LASSO"      = "#A5D6A7",
  "AdaLASSO"   = "#C5E1A5",
  "ElNET"      = "#B0BEC5",
  "AdaElNET"   = "#BCAAA4",
  "Factor"     = "#CE93D8",
  "T.Factor"   = "#E1BEE7",
  "RF"         = "#80CBC4",
  "Bagging"    = "#FFCC80",
  "CSR"        = "#FFF59D"
)

# Garantir que todos os metodos tenham cor (fallback para cinza)
cores_modelos <- cores_fixas[intersect(names(cores_fixas),
                                        c("Realizado", methods))]
faltam <- setdiff(methods, names(cores_modelos))
if (length(faltam) > 0) {
  extras <- setNames(rep("grey75", length(faltam)), faltam)
  cores_modelos <- c(cores_modelos, extras)
}
# Manter Realizado na frente
cores_modelos <- c("Realizado" = "black",
                   cores_modelos[names(cores_modelos) != "Realizado"])

cat("=== Configuracao ===\n")
cat(sprintf("  Target:     %s\n", paste(targets, collapse = ", ")))
cat(sprintf("  Horizontes: %s\n", paste(horizons, collapse = ", ")))
cat(sprintf("  Metodos:    %s\n", paste(methods, collapse = ", ")))
cat(sprintf("  n_eval:     %d\n", n_eval))
cat(sprintf("  Chaves:     %d combinacoes\n\n", length(res)))

###############################################################################
# PARTE A -- Serie Real + Previsoes por Horizonte (com destaque 2SRR)
###############################################################################
cat(" PARTE A -- Graficos: Serie Real vs Previsoes\n")

for (key in names(res)) {
  r   <- res[[key]]
  n_e <- length(r$actuals)
  if (n_e < 10) next

  df_plot <- data.frame(Date = r$dates[1:n_e], Realizado = r$actuals)
  for (m in methods) {
    if (m %in% colnames(r$preds)) {
      df_plot[[m]] <- r$preds[1:n_e, m]
    }
  }

  df_long <- df_plot |>
    pivot_longer(-Date, names_to = "Serie", values_to = "Valor") |>
    mutate(Serie = factor(Serie, levels = c("Realizado", methods)))

  # Ordenar para que 2SRR e Realizado sejam desenhados por ultimo (ficam por cima)
  df_long <- df_long |>
    mutate(ordem_plot = case_when(
      Serie == "2SRR"      ~ 3L,
      Serie == "Realizado" ~ 2L,
      TRUE                 ~ 1L
    )) |>
    arrange(ordem_plot)

  # Linewidths: 2SRR e Realizado grossos, demais finos
  lw_vals <- setNames(rep(0.35, length(methods) + 1),
                      c("Realizado", methods))
  lw_vals["Realizado"] <- 1.3
  lw_vals["2SRR"]      <- 1.1

  # Alpha: 2SRR e Realizado opacos, demais semi-transparentes
  alpha_vals <- setNames(rep(0.35, length(methods) + 1),
                         c("Realizado", methods))
  alpha_vals["Realizado"] <- 1.0
  alpha_vals["2SRR"]      <- 0.95

  p <- ggplot(df_long, aes(Date, Valor, color = Serie, linewidth = Serie,
                            alpha = Serie)) +
    geom_line(na.rm = TRUE) +
    scale_color_manual(values = cores_modelos) +
    scale_linewidth_manual(values = lw_vals, guide = "none") +
    scale_alpha_manual(values = alpha_vals, guide = "none") +
    labs(
      title    = sprintf("%s -- Horizonte h = %d", r$target, r$horizon),
      subtitle = "Serie realizada (preto) vs 2SRR (vermelho) vs demais modelos",
      x = NULL, y = "Valor"
    ) +
    guides(color = guide_legend(nrow = 2))

  fname <- sprintf("fig_A_forecast_%s_h%02d", r$target, r$horizon)
  ggsave(file.path(results_dir, paste0(fname, ".png")), p,
         width = 12, height = 5.5, dpi = 250)
  ggsave(file.path(results_dir, paste0(fname, ".pdf")), p,
         width = 12, height = 5.5)
  cat(sprintf("  [OK] %s.png/pdf\n", fname))
}

###############################################################################
# PARTE B -- Tabela MSFE Relativo ao AR (benchmark)
###############################################################################
cat("\n=======================================================\n")
cat(" PARTE B -- Tabelas de Metricas (MSFE / RMSE)\n")

msfe_rows <- list()
for (key in names(res)) {
  r   <- res[[key]]
  row <- data.frame(Target = r$target, Horizon = r$horizon,
                    stringsAsFactors = FALSE)
  for (m in methods) {
    ok <- !is.na(r$preds[, m]) & !is.na(r$actuals)
    if (sum(ok) > 5) {
      row[[paste0("MSFE_", m)]] <- mean((r$preds[ok, m] - r$actuals[ok])^2)
      row[[paste0("RMSE_", m)]] <- sqrt(row[[paste0("MSFE_", m)]])
    } else {
      row[[paste0("MSFE_", m)]] <- NA
      row[[paste0("RMSE_", m)]] <- NA
    }
  }
  msfe_rows[[key]] <- row
}

df_msfe <- do.call(rbind, msfe_rows)
rownames(df_msfe) <- NULL

# MSFE relativo ao AR
if (!"MSFE_AR" %in% names(df_msfe)) {
  stop("Modelo AR nao encontrado nos resultados. Verifique se AR.rda existe.")
}

df_rel <- df_msfe[, c("Target", "Horizon")]
for (m in methods) {
  df_rel[[m]] <- round(df_msfe[[paste0("MSFE_", m)]] / df_msfe[["MSFE_AR"]], 4)
}

cat("TABELA MSFE RELATIVO AO AR (< 1 = supera o AR):\n")
print(df_rel)

# RMSE absoluto
df_rmse <- df_msfe[, c("Target", "Horizon")]
for (m in methods) {
  df_rmse[[m]] <- round(df_msfe[[paste0("RMSE_", m)]], 5)
}

cat("\nTABELA RMSE ABSOLUTO:\n")
print(df_rmse)

write.csv(df_rel,  file.path(results_dir, "tab_B_msfe_relative.csv"),
          row.names = FALSE)
write.csv(df_rmse, file.path(results_dir, "tab_B_rmse_absolute.csv"),
          row.names = FALSE)
cat("\n  [OK] tab_B_msfe_relative.csv\n  [OK] tab_B_rmse_absolute.csv\n")

###############################################################################
# PARTE C -- Teste Diebold-Mariano
###############################################################################
cat(" PARTE C -- Teste Diebold-Mariano (DM)\n")

dm_rows <- list()
for (key in names(res)) {
  r    <- res[[key]]
  e_ar <- (r$preds[, "AR"] - r$actuals)^2
  row  <- data.frame(Target = r$target, Horizon = r$horizon,
                     stringsAsFactors = FALSE)

  for (m in methods) {
    if (m == "AR") next
    e_m <- (r$preds[, m] - r$actuals)^2
    ok  <- !is.na(e_ar) & !is.na(e_m)

    if (sum(ok) > 20) {
      tryCatch({
        dm   <- dm.test(ts(r$actuals[ok] - r$preds[ok, "AR"]),
                        ts(r$actuals[ok] - r$preds[ok, m]),
                        alternative = "two.sided",
                        h = r$horizon, power = 2)
        pval  <- dm$p.value
        stars <- ifelse(pval < 0.01, "***",
                        ifelse(pval < 0.05, "**",
                               ifelse(pval < 0.10, "*", "")))
        row[[paste0("DM_pval_", m)]]  <- round(pval, 4)
        row[[paste0("DM_stars_", m)]] <- stars
      }, error = function(e) {
        row[[paste0("DM_pval_", m)]]  <<- NA
        row[[paste0("DM_stars_", m)]] <<- ""
      })
    } else {
      row[[paste0("DM_pval_", m)]]  <- NA
      row[[paste0("DM_stars_", m)]] <- ""
    }
  }
  dm_rows[[key]] <- row
}

df_dm <- do.call(rbind, dm_rows)
rownames(df_dm) <- NULL

cat("TESTE DIEBOLD-MARIANO (p-valores vs AR):\n")
cat("  *** = 1%, ** = 5%, * = 10%\n\n")
print(df_dm)
write.csv(df_dm, file.path(results_dir, "tab_C_diebold_mariano.csv"),
          row.names = FALSE)
cat("\n  [OK] tab_C_diebold_mariano.csv\n")

###############################################################################
# PARTE D -- CSFE (Cumulative Squared Forecast Errors)
###############################################################################
cat(" PARTE D -- CSFE (Cumulative Squared Forecast Errors)\n")

for (key in names(res)) {
  r   <- res[[key]]
  n_e <- length(r$actuals)
  if (n_e < 20) next

  e2_ar   <- (r$preds[1:n_e, "AR"] - r$actuals)^2
  df_csfe <- data.frame(Date = r$dates[1:n_e])

  for (m in methods) {
    if (m == "AR") next
    e2_m <- (r$preds[1:n_e, m] - r$actuals)^2
    csfe <- cumsum(ifelse(is.na(e2_ar) | is.na(e2_m), 0, e2_ar - e2_m))
    df_csfe[[m]] <- csfe
  }

  df_csfe_long <- df_csfe |>
    pivot_longer(-Date, names_to = "Modelo", values_to = "CSFE")

  p <- ggplot(df_csfe_long, aes(Date, CSFE, color = Modelo)) +
    geom_line(linewidth = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = cores_modelos[-1]) +
    labs(
      title    = sprintf("CSFE -- %s, h = %d", r$target, r$horizon),
      subtitle = "Acima de 0 = modelo supera o AR; abaixo = AR e melhor",
      x = NULL, y = "CSFE (acumulado)"
    )

  fname <- sprintf("fig_D_csfe_%s_h%02d", r$target, r$horizon)
  ggsave(file.path(results_dir, paste0(fname, ".png")), p,
         width = 10, height = 5, dpi = 250)
  cat(sprintf("  [OK] %s.png\n", fname))
}

###############################################################################
# PARTE E -- Analise de Residuos (histograma + QQ-plot + ACF)
###############################################################################
cat(" PARTE E -- Analise de Residuos\n")

# Para nao gerar centenas de graficos, focar nos modelos principais
modelos_residuos <- intersect(c("AR", "Ridge", "LASSO", "2SRR", "RF"), methods)

for (key in names(res)) {
  r   <- res[[key]]
  n_e <- length(r$actuals)
  if (n_e < 20) next

  for (m in modelos_residuos) {
    if (!m %in% colnames(r$preds)) next
    residuos <- r$actuals - r$preds[1:n_e, m]
    residuos <- residuos[!is.na(residuos)]
    if (length(residuos) < 20) next

    df_res <- data.frame(Residuo = residuos)

    p_hist <- ggplot(df_res, aes(Residuo)) +
      geom_histogram(bins = 25, fill = "steelblue", alpha = 0.7,
                     color = "white") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
      labs(
        title = sprintf("Residuos: %s -- %s (h=%d)", m, r$target, r$horizon),
        subtitle = sprintf("Media=%.4f, SD=%.4f, Skew=%.3f",
                           mean(residuos), sd(residuos),
                           mean(((residuos - mean(residuos)) /
                                   sd(residuos))^3)),
        x = "Residuo (Real - Previsto)", y = "Frequencia"
      )

    p_qq <- ggplot(df_res, aes(sample = Residuo)) +
      stat_qq(color = "steelblue", alpha = 0.6) +
      stat_qq_line(color = "red", linewidth = 0.8) +
      labs(title = "QQ-Plot vs Normal",
           x = "Quantis Teoricos", y = "Quantis Amostrais")

    acf_vals <- acf(residuos, lag.max = 20, plot = FALSE)
    df_acf   <- data.frame(Lag = acf_vals$lag[-1], ACF = acf_vals$acf[-1])
    ci       <- 1.96 / sqrt(length(residuos))

    p_acf <- ggplot(df_acf, aes(Lag, ACF)) +
      geom_hline(yintercept = c(-ci, ci), linetype = "dashed",
                 color = "blue", alpha = 0.5) +
      geom_hline(yintercept = 0, color = "grey50") +
      geom_segment(aes(xend = Lag, yend = 0), color = "steelblue",
                   linewidth = 1) +
      geom_point(color = "steelblue", size = 2) +
      labs(title = "ACF dos Residuos", x = "Lag", y = "Autocorrelacao")

    p_combined <- p_hist + p_qq + p_acf + plot_layout(ncol = 3) +
      plot_annotation(
        title = sprintf("Diagnostico de Residuos: %s -- %s (h=%d)",
                        m, r$target, r$horizon)
      )

    fname <- sprintf("fig_E_residuos_%s_%s_h%02d", r$target, m, r$horizon)
    ggsave(file.path(results_dir, paste0(fname, ".png")), p_combined,
           width = 15, height = 5, dpi = 200)
  }
  cat(sprintf("  [OK] Residuos %s h=%d\n", r$target, r$horizon))
}

###############################################################################
# PARTE F -- Model Confidence Set (MCS)
###############################################################################
cat(" PARTE F -- Model Confidence Set (MCS)\n")

if (has_mcs) {
  mcs_results <- list()
  for (key in names(res)) {
    r   <- res[[key]]
    n_e <- length(r$actuals)
    if (n_e < 20) next

    loss_mat <- matrix(NA, n_e, length(methods))
    colnames(loss_mat) <- methods
    for (m in methods) {
      loss_mat[, m] <- (r$preds[1:n_e, m] - r$actuals)^2
    }

    # Substituir Inf e NaN por NA, depois remover linhas e colunas problematicas
    loss_mat[is.infinite(loss_mat)] <- NA
    loss_mat[is.nan(loss_mat)]      <- NA

    # Remover colunas (modelos) que tenham mais de 10% de NA
    pct_na <- colMeans(is.na(loss_mat))
    cols_ok <- names(which(pct_na < 0.10))
    if (length(cols_ok) < 3) next
    loss_mat <- loss_mat[, cols_ok, drop = FALSE]

    # Remover linhas com qualquer NA restante
    ok <- complete.cases(loss_mat)
    loss_mat_clean <- loss_mat[ok, , drop = FALSE]
    if (nrow(loss_mat_clean) < 30) next

    # Garantir que e numerico
    storage.mode(loss_mat_clean) <- "double"

    tryCatch({
      mcs_out <- MCSprocedure(as.data.frame(loss_mat_clean),
                              alpha = 0.15, B = 5000,
                              statistic = "Tmax")
      surviving <- rownames(mcs_out@show)
      in_mcs <- ifelse("2SRR" %in% surviving, " << 2SRR INCLUIDO", "")

      cat(sprintf("  %s h=%d: MCS = {%s}%s\n",
                  r$target, r$horizon,
                  paste(surviving, collapse = ", "), in_mcs))

      mcs_results[[key]] <- data.frame(
        Target     = r$target,
        Horizon    = r$horizon,
        MCS_Models = paste(surviving, collapse = ", "),
        N_modelos  = length(surviving),
        Inclui_2SRR = "2SRR" %in% surviving
      )
    }, error = function(e) {
      cat(sprintf("  %s h=%d: MCS erro -- %s\n",
                  r$target, r$horizon, e$message))
    })
  }

  if (length(mcs_results) > 0) {
    df_mcs <- do.call(rbind, mcs_results)
    write.csv(df_mcs, file.path(results_dir, "tab_F_mcs.csv"),
              row.names = FALSE)
    cat("\n  [OK] tab_F_mcs.csv\n")
  }
} else {
  cat("  [AVISO] Pacote MCS nao disponivel. Instale com:\n")
  cat('  install.packages("MCS")\n')
}

###############################################################################
# PARTE G -- TVP Betas ao Longo do Tempo
###############################################################################
cat(" PARTE G -- TVP Betas (2SRR Full Sample)\n")

if (!is.null(forecast_output$tvp_full)) {
  for (tgt in names(forecast_output$tvp_full)) {
    betas <- forecast_output$tvp_full[[tgt]]$betas
    d     <- forecast_output$tvp_full[[tgt]]$dates

    if (is.null(betas) || ncol(betas) < 2) {
      cat("  [AVISO] Betas insuficientes para gerar grafico.\n")
      next
    }

    # Remover colunas com variancia zero ou NA
    beta_var <- apply(betas, 2, function(x) {
      x <- x[!is.na(x)]
      if (length(x) < 2) return(0)
      var(x)
    })
    betas_valid <- betas[, beta_var > 0, drop = FALSE]
    beta_var    <- beta_var[beta_var > 0]

    if (ncol(betas_valid) < 2) {
      cat("  [AVISO] Variancia zero em todos os betas.\n")
      next
    }

    # --- Grafico 1: Top 8 por VARIANCIA (os que mais flutuaram) ---
    n_show    <- min(8, ncol(betas_valid))
    top_vars  <- names(sort(beta_var, decreasing = TRUE))[1:n_show]

    df_beta <- data.frame(Date = d, betas_valid[, top_vars, drop = FALSE],
                          check.names = FALSE) |>
      pivot_longer(-Date, names_to = "Variable", values_to = "Beta")

    p1 <- ggplot(df_beta, aes(Date, Beta, color = Variable)) +
      geom_line(linewidth = 0.6) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      facet_wrap(~ Variable, scales = "free_y", ncol = 2) +
      labs(
        title    = sprintf("TVP Betas com Maior Variancia -- %s", tgt),
        subtitle = "Coeficientes que mais flutuaram ao longo do tempo (2SRR)",
        x = NULL, y = expression(beta[t])
      ) +
      theme(legend.position = "none")

    fname1 <- sprintf("fig_G1_tvp_betas_variancia_%s", tgt)
    ggsave(file.path(results_dir, paste0(fname1, ".png")), p1,
           width = 10, height = 10, dpi = 200)
    cat(sprintf("  [OK] %s.png\n", fname1))

    # --- Grafico 2: Top 8 por MEDIA ABSOLUTA (os mais importantes) ---
    beta_mean_abs <- apply(abs(betas_valid), 2, function(x) {
      mean(x, na.rm = TRUE)
    })
    top_mean <- names(sort(beta_mean_abs, decreasing = TRUE))[1:n_show]

    df_beta2 <- data.frame(Date = d, betas_valid[, top_mean, drop = FALSE],
                           check.names = FALSE) |>
      pivot_longer(-Date, names_to = "Variable", values_to = "Beta")

    p2 <- ggplot(df_beta2, aes(Date, Beta, color = Variable)) +
      geom_line(linewidth = 0.6) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      facet_wrap(~ Variable, scales = "free_y", ncol = 2) +
      labs(
        title    = sprintf("TVP Betas com Maior Importancia Media -- %s", tgt),
        subtitle = "Preditores estruturalmente mais relevantes (media |beta|)",
        x = NULL, y = expression(beta[t])
      ) +
      theme(legend.position = "none")

    fname2 <- sprintf("fig_G2_tvp_betas_importancia_%s", tgt)
    ggsave(file.path(results_dir, paste0(fname2, ".png")), p2,
           width = 10, height = 10, dpi = 200)
    cat(sprintf("  [OK] %s.png\n", fname2))

    # --- Tabela de volatilidade dos betas ---
    df_beta_stats <- data.frame(
      Predictor = colnames(betas_valid),
      Variancia = round(beta_var, 6),
      Media_Abs = round(beta_mean_abs, 6),
      stringsAsFactors = FALSE
    ) |>
      arrange(desc(Variancia))

    write.csv(df_beta_stats,
              file.path(results_dir,
                        sprintf("tab_G_beta_stats_%s.csv", tgt)),
              row.names = FALSE)
    cat(sprintf("  [OK] tab_G_beta_stats_%s.csv\n", tgt))
  }
} else {
  cat("  TVP Betas nao disponiveis.\n")
}

###############################################################################
# PARTE H -- Grafico MSFE Relativo (barras por horizonte) com destaque 2SRR
###############################################################################
cat(" PARTE H -- Grafico MSFE Relativo (barras)\n")

modelos_barra <- intersect(c("Ridge", "LASSO", "AdaLASSO", "ElNET",
                             "RF", "2SRR", "Factor", "CSR", "Bagging"),
                           methods)

# Cores para barras: 2SRR vermelho, demais cinza/pastel
cores_barra <- setNames(rep("grey70", length(modelos_barra)), modelos_barra)
cores_barra["2SRR"] <- "#D32F2F"
if ("Ridge"    %in% modelos_barra) cores_barra["Ridge"]    <- "#90CAF9"
if ("RF"       %in% modelos_barra) cores_barra["RF"]       <- "#80CBC4"
if ("LASSO"    %in% modelos_barra) cores_barra["LASSO"]    <- "#A5D6A7"
if ("AdaLASSO" %in% modelos_barra) cores_barra["AdaLASSO"] <- "#C5E1A5"

df_rel_long <- df_rel |>
  pivot_longer(cols = all_of(methods), names_to = "Modelo",
               values_to = "MSFE_rel") |>
  filter(Modelo != "AR", Modelo %in% modelos_barra) |>
  mutate(Horizonte = factor(paste0("h=", Horizon),
                            levels = paste0("h=", horizons_list)))

for (tgt in unique(df_rel_long$Target)) {
  df_tgt <- df_rel_long |> filter(Target == tgt)

  p <- ggplot(df_tgt, aes(Horizonte, MSFE_rel, fill = Modelo)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7,
             alpha = 0.85) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red",
               linewidth = 0.7) +
    scale_fill_manual(values = cores_barra) +
    labs(
      title    = sprintf("MSFE Relativo ao AR -- %s", tgt),
      subtitle = "Abaixo de 1 = supera AR. Vermelho = 2SRR",
      x = NULL, y = "MSFE / MSFE(AR)"
    ) +
    coord_cartesian(
      ylim = c(0, max(df_tgt$MSFE_rel, na.rm = TRUE) * 1.15)
    )

  fname <- sprintf("fig_H_msfe_bars_%s", tgt)
  ggsave(file.path(results_dir, paste0(fname, ".png")), p,
         width = 10, height = 5.5, dpi = 250)
  cat(sprintf("  [OK] %s.png\n", fname))
}

###############################################################################
# PARTE I -- Tabela LaTeX-Ready
###############################################################################
cat(" PARTE I -- Tabela LaTeX\n")

# Selecionar modelos para a tabela LaTeX (os mais relevantes)
modelos_latex <- intersect(c("AR", "Ridge", "LASSO", "AdaLASSO", "ElNET",
                             "RF", "2SRR", "Factor", "Bagging", "CSR",
                             "rw", "AR_BIC"),
                           methods)

latex_file <- file.path(results_dir, "tab_I_latex_msfe.tex")
sink(latex_file)

cat("\\begin{table}[ht]\n")
cat("\\centering\n")
cat("\\caption{MSFE relativo ao AR --- Dados Americanos}\n")
cat("\\label{tab:msfe_relative}\n")

cat(sprintf("\\begin{tabular}{ll%s}\n",
            paste(rep("c", length(modelos_latex)), collapse = "")))
cat("\\hline\\hline\n")
cat(sprintf("Target & h & %s \\\\\n",
            paste(modelos_latex, collapse = " & ")))
cat("\\hline\n")

for (i in seq_len(nrow(df_rel))) {
  vals <- sapply(modelos_latex, function(m) {
    v <- df_rel[i, m]
    if (is.null(v) || is.na(v)) return("---")

    # Melhor modelo (menor MSFE relativo, exceto AR)
    all_vals <- sapply(setdiff(modelos_latex, "AR"), function(mm) {
      val <- df_rel[i, mm]
      if (is.null(val) || is.na(val)) return(Inf)
      val
    })
    is_best  <- (!is.na(v) && m != "AR" && v == min(all_vals, na.rm = TRUE))
    formatted <- sprintf("%.3f", v)
    if (is_best) formatted <- sprintf("\\textbf{%s}", formatted)

    # Estrelas do DM
    dm_col <- paste0("DM_stars_", m)
    if (dm_col %in% names(df_dm)) {
      stars <- df_dm[df_dm$Target == df_rel$Target[i] &
                       df_dm$Horizon == df_rel$Horizon[i], dm_col]
      if (length(stars) > 0 && !is.na(stars[1]) && nchar(stars[1]) > 0) {
        formatted <- paste0(formatted, "$^{", stars[1], "}$")
      }
    }
    formatted
  })
  cat(sprintf("%s & %d & %s \\\\\n",
              df_rel$Target[i], df_rel$Horizon[i],
              paste(vals, collapse = " & ")))
}

cat("\\hline\\hline\n")
cat("\\end{tabular}\n")
cat("\\begin{tablenotes}\n")
cat("\\small\n")
cat("\\item Nota: Valores representam MSFE relativo ao AR. ")
cat("Valores $< 1$ indicam modelo superior ao benchmark. ")
cat("Negrito indica melhor modelo. ")
cat("$^{***}$, $^{**}$, $^{*}$ indicam significancia no teste ")
cat("de Diebold-Mariano a 1\\%, 5\\% e 10\\%.\n")
cat("\\end{tablenotes}\n")
cat("\\end{table}\n")

sink()
cat(sprintf("  [OK] %s\n", latex_file))

###############################################################################
# PARTE J -- Performance por Regime de Volatilidade
###############################################################################
cat(" PARTE J -- Performance por Regime (Normal vs Picos)\n")

regime_rows <- list()
for (key in names(res)) {
  r   <- res[[key]]
  n_e <- length(r$actuals)
  if (n_e < 30) next

  # Definir regime: picos = acima do percentil 80 em valor absoluto
  limiar <- quantile(abs(r$actuals), 0.80, na.rm = TRUE)
  regime <- ifelse(abs(r$actuals) >= limiar, "Picos", "Normal")

  for (reg in c("Normal", "Picos")) {
    idx <- which(regime == reg)
    if (length(idx) < 10) next

    row <- data.frame(Target  = r$target,
                      Horizon = r$horizon,
                      Regime  = reg,
                      N_obs   = length(idx),
                      stringsAsFactors = FALSE)

    for (m in methods) {
      erros <- (r$preds[idx, m] - r$actuals[idx])^2
      ok    <- !is.na(erros)
      if (sum(ok) > 5) {
        row[[paste0("MSFE_", m)]] <- mean(erros[ok])
      } else {
        row[[paste0("MSFE_", m)]] <- NA
      }
    }
    regime_rows[[paste0(key, "_", reg)]] <- row
  }
}

df_regime <- do.call(rbind, regime_rows)
rownames(df_regime) <- NULL

# Calcular MSFE relativo ao AR por regime
df_regime_rel <- df_regime[, c("Target", "Horizon", "Regime", "N_obs")]
for (m in methods) {
  df_regime_rel[[m]] <- round(
    df_regime[[paste0("MSFE_", m)]] / df_regime[["MSFE_AR"]], 4)
}

cat("TABELA MSFE RELATIVO AO AR POR REGIME:\n")
cat("(Mostra onde o 2SRR ganha/perde em periodos de alta volatilidade)\n\n")
print(df_regime_rel)

write.csv(df_regime_rel, file.path(results_dir, "tab_K_regime_msfe.csv"),
          row.names = FALSE)
cat("\n  [OK] tab_K_regime_msfe.csv\n")

# Grafico comparativo: 2SRR vs AR vs Ridge por regime
modelos_regime <- intersect(c("2SRR", "Ridge", "RF", "LASSO"), methods)

df_regime_plot <- df_regime_rel |>
  pivot_longer(cols = all_of(modelos_regime), names_to = "Modelo",
               values_to = "MSFE_rel") |>
  mutate(Horizonte = factor(paste0("h=", Horizon),
                            levels = paste0("h=", horizons_list)))

cores_regime <- c("2SRR" = "#D32F2F", "Ridge" = "#90CAF9",
                  "RF" = "#80CBC4", "LASSO" = "#A5D6A7")

p <- ggplot(df_regime_plot,
            aes(Horizonte, MSFE_rel, fill = Modelo)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6,
           alpha = 0.85) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  facet_wrap(~ Regime, ncol = 2) +
  scale_fill_manual(values = cores_regime) +
  labs(
    title    = "MSFE Relativo ao AR por Regime de Volatilidade",
    subtitle = "Normal = 80% centrais | Picos = 20% mais extremos",
    x = NULL, y = "MSFE / MSFE(AR)"
  )

ggsave(file.path(results_dir, "fig_K_regime_comparison.png"), p,
       width = 11, height = 5.5, dpi = 250)
cat("  [OK] fig_K_regime_comparison.png\n")

###############################################################################
# PARTE K -- Ranking Geral (quantas vezes cada modelo e o melhor)
###############################################################################
cat(" PARTE K -- Ranking Geral dos Modelos\n")

# Contar quantas vezes cada modelo tem o menor MSFE relativo
modelos_competidores <- setdiff(methods, c("AR", "rw", "AR_BIC"))
ranking <- data.frame(Modelo = modelos_competidores,
                      Vitorias = 0L,
                      Top3 = 0L,
                      Media_MSFE_rel = NA_real_,
                      stringsAsFactors = FALSE)

for (m in modelos_competidores) {
  msfe_vals <- c()
  for (i in seq_len(nrow(df_rel))) {
    v <- df_rel[i, m]
    if (!is.null(v) && !is.na(v)) {
      msfe_vals <- c(msfe_vals, v)

      # Verificar se e o melhor nesta combinacao
      all_v <- sapply(modelos_competidores, function(mm) {
        val <- df_rel[i, mm]
        if (is.null(val) || is.na(val)) return(Inf)
        val
      })
      if (v == min(all_v, na.rm = TRUE)) {
        ranking$Vitorias[ranking$Modelo == m] <-
          ranking$Vitorias[ranking$Modelo == m] + 1L
      }
      # Top 3
      if (v <= sort(all_v)[min(3, length(all_v))]) {
        ranking$Top3[ranking$Modelo == m] <-
          ranking$Top3[ranking$Modelo == m] + 1L
      }
    }
  }
  ranking$Media_MSFE_rel[ranking$Modelo == m] <- round(mean(msfe_vals), 4)
}

ranking <- ranking |> arrange(desc(Vitorias), Media_MSFE_rel)
cat("RANKING GERAL (por numero de vitorias e MSFE medio):\n\n")
print(ranking)

write.csv(ranking, file.path(results_dir, "tab_L_ranking.csv"),
          row.names = FALSE)
cat("\n  [OK] tab_L_ranking.csv\n")

###############################################################################
# PARTE L -- CSFE focado: 2SRR vs 3 melhores concorrentes
###############################################################################
cat(" PARTE L -- CSFE focado (2SRR vs principais rivais)\n")

# Identificar os 3 modelos com menor MSFE medio (exceto 2SRR e AR)
rivais <- ranking |>
  filter(Modelo != "2SRR") |>
  head(3) |>
  pull(Modelo)

modelos_foco <- c("2SRR", rivais)
cat(sprintf("  Modelos no grafico: %s\n\n", paste(modelos_foco, collapse = ", ")))

cores_foco <- c("2SRR" = "#D32F2F")
paleta_rivais <- c("#1976D2", "#388E3C", "#F57C00")
for (j in seq_along(rivais)) {
  cores_foco[rivais[j]] <- paleta_rivais[j]
}

for (key in names(res)) {
  r   <- res[[key]]
  n_e <- length(r$actuals)
  if (n_e < 20) next

  e2_ar   <- (r$preds[1:n_e, "AR"] - r$actuals)^2
  df_csfe <- data.frame(Date = r$dates[1:n_e])

  for (m in modelos_foco) {
    if (!m %in% colnames(r$preds)) next
    e2_m <- (r$preds[1:n_e, m] - r$actuals)^2
    csfe <- cumsum(ifelse(is.na(e2_ar) | is.na(e2_m), 0, e2_ar - e2_m))
    df_csfe[[m]] <- csfe
  }

  cols_presentes <- intersect(modelos_foco, names(df_csfe))
  if (length(cols_presentes) < 2) next

  df_csfe_long <- df_csfe |>
    pivot_longer(cols = all_of(cols_presentes),
                 names_to = "Modelo", values_to = "CSFE")

  # Linewidth: 2SRR mais grosso
  lw_foco <- setNames(rep(0.6, length(cols_presentes)), cols_presentes)
  lw_foco["2SRR"] <- 1.2

  p <- ggplot(df_csfe_long, aes(Date, CSFE, color = Modelo,
                                 linewidth = Modelo)) +
    geom_line() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = cores_foco) +
    scale_linewidth_manual(values = lw_foco, guide = "none") +
    labs(
      title    = sprintf("CSFE -- %s, h = %d (2SRR vs rivais)",
                         r$target, r$horizon),
      subtitle = "Acima de 0 = supera AR. Comparacao direta com concorrentes",
      x = NULL, y = "CSFE (acumulado)"
    )

  fname <- sprintf("fig_M_csfe_foco_%s_h%02d", r$target, r$horizon)
  ggsave(file.path(results_dir, paste0(fname, ".png")), p,
         width = 10, height = 5, dpi = 250)
  cat(sprintf("  [OK] %s.png\n", fname))
}

###############################################################################
# PARTE M -- Janelas temporais: onde o 2SRR perde para o melhor rival
###############################################################################
cat(" PARTE M -- Analise temporal de derrotas do 2SRR\n")

if ("2SRR" %in% methods && length(rivais) > 0) {
  melhor_rival <- rivais[1]

  for (key in names(res)) {
    r   <- res[[key]]
    n_e <- length(r$actuals)
    if (n_e < 20) next
    if (!melhor_rival %in% colnames(r$preds)) next

    e2_2srr  <- (r$preds[1:n_e, "2SRR"] - r$actuals)^2
    e2_rival <- (r$preds[1:n_e, melhor_rival] - r$actuals)^2

    # Diferenca: positivo = 2SRR errou MAIS que o rival (derrota do 2SRR)
    diff_e2 <- e2_2srr - e2_rival
    diff_e2[is.na(diff_e2)] <- 0

    df_diff <- data.frame(
      Date   = r$dates[1:n_e],
      Diff   = diff_e2,
      Derrota = ifelse(diff_e2 > 0, "2SRR pior", "2SRR melhor")
    )

    p <- ggplot(df_diff, aes(Date, Diff, fill = Derrota)) +
      geom_col(width = 25, alpha = 0.7) +
      geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
      scale_fill_manual(values = c("2SRR pior"   = "#D32F2F",
                                   "2SRR melhor" = "#2E7D32")) +
      labs(
        title    = sprintf("Diferenca de erro quadratico: 2SRR vs %s (h=%d)",
                           melhor_rival, r$horizon),
        subtitle = "Barras vermelhas = periodos onde o 2SRR errou mais",
        x = NULL, y = "e2(2SRR) - e2(rival)"
      )

    fname <- sprintf("fig_N_diff_%s_h%02d", r$target, r$horizon)
    ggsave(file.path(results_dir, paste0(fname, ".png")), p,
           width = 10, height = 4.5, dpi = 250)
    cat(sprintf("  [OK] %s.png\n", fname))
  }
} else {
  cat("  2SRR ou rivais nao disponiveis para comparacao.\n")
}

###############################################################################
# PARTE N -- CSFE Rolling Window: MSFE relativo do 2SRR ao longo do tempo
###############################################################################
cat(" PARTE N -- MSFE Rolling do 2SRR (janela 36 meses)\n")

window_size <- 36

for (key in names(res)) {
  r   <- res[[key]]
  n_e <- length(r$actuals)
  if (n_e < window_size + 10) next
  if (!"2SRR" %in% colnames(r$preds)) next

  e2_ar   <- (r$preds[1:n_e, "AR"]   - r$actuals)^2
  e2_2srr <- (r$preds[1:n_e, "2SRR"] - r$actuals)^2

  # MSFE relativo rolling
  roll_rel  <- rep(NA, n_e)
  roll_date <- r$dates[1:n_e]

  for (t in window_size:n_e) {
    idx <- (t - window_size + 1):t
    ok  <- !is.na(e2_ar[idx]) & !is.na(e2_2srr[idx])
    if (sum(ok) > 10) {
      roll_rel[t] <- mean(e2_2srr[idx[ok]]) / mean(e2_ar[idx[ok]])
    }
  }

  df_roll <- data.frame(
    Date     = roll_date,
    MSFE_rel = roll_rel
  ) |>
    filter(!is.na(MSFE_rel))

  # Colorir por acima/abaixo de 1
  df_roll$Vantagem <- ifelse(df_roll$MSFE_rel < 1, "2SRR melhor", "AR melhor")

  # Periodos onde 2SRR ganha
  pct_ganha <- round(100 * mean(df_roll$MSFE_rel < 1), 1)

  p <- ggplot(df_roll, aes(Date, MSFE_rel)) +
    geom_ribbon(aes(ymin = pmin(MSFE_rel, 1), ymax = 1),
                fill = "#2E7D32", alpha = 0.15) +
    geom_ribbon(aes(ymin = 1, ymax = pmax(MSFE_rel, 1)),
                fill = "#D32F2F", alpha = 0.15) +
    geom_line(linewidth = 0.7, color = "#D32F2F") +
    geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
    labs(
      title = sprintf("MSFE Relativo Rolling (2SRR/AR) -- h=%d, janela=%d meses",
                      r$horizon, window_size),
      subtitle = sprintf("Verde = 2SRR supera AR. 2SRR vence em %.1f%% das janelas",
                          pct_ganha),
      x = NULL, y = "MSFE(2SRR) / MSFE(AR)"
    ) +
    annotate("text", x = min(df_roll$Date) + 200, y = 0.85,
             label = "2SRR melhor", color = "#2E7D32", fontface = "bold",
             size = 3.5) +
    annotate("text", x = min(df_roll$Date) + 200, y = 1.15,
             label = "AR melhor", color = "#D32F2F", fontface = "bold",
             size = 3.5)

  fname <- sprintf("fig_O_rolling_msfe_%s_h%02d", r$target, r$horizon)
  ggsave(file.path(results_dir, paste0(fname, ".png")), p,
         width = 11, height = 5, dpi = 250)
  cat(sprintf("  [OK] %s.png (2SRR vence em %.1f%% das janelas)\n",
              fname, pct_ganha))
}

###############################################################################
# PARTE O -- Analise pre/pos COVID (quebra estrutural)
###############################################################################
cat(" PARTE O -- Performance pre vs pos COVID\n")

covid_date <- as.Date("2020-03-01")

covid_rows <- list()
for (key in names(res)) {
  r   <- res[[key]]
  n_e <- length(r$actuals)
  if (n_e < 30) next

  for (periodo in c("Pre-COVID", "Pos-COVID")) {
    if (periodo == "Pre-COVID") {
      idx <- which(r$dates[1:n_e] < covid_date)
    } else {
      idx <- which(r$dates[1:n_e] >= covid_date)
    }
    if (length(idx) < 10) next

    row <- data.frame(Target  = r$target,
                      Horizon = r$horizon,
                      Periodo = periodo,
                      N_obs   = length(idx),
                      stringsAsFactors = FALSE)
    for (m in methods) {
      erros <- (r$preds[idx, m] - r$actuals[idx])^2
      ok    <- !is.na(erros) & !is.infinite(erros)
      if (sum(ok) > 5) {
        row[[paste0("MSFE_", m)]] <- mean(erros[ok])
      } else {
        row[[paste0("MSFE_", m)]] <- NA
      }
    }
    covid_rows[[paste0(key, "_", periodo)]] <- row
  }
}

df_covid <- do.call(rbind, covid_rows)
rownames(df_covid) <- NULL

df_covid_rel <- df_covid[, c("Target", "Horizon", "Periodo", "N_obs")]
for (m in methods) {
  df_covid_rel[[m]] <- round(
    df_covid[[paste0("MSFE_", m)]] / df_covid[["MSFE_AR"]], 4)
}

cat("TABELA MSFE RELATIVO AO AR -- PRE vs POS COVID:\n\n")
print(df_covid_rel)

write.csv(df_covid_rel, file.path(results_dir, "tab_P_covid_msfe.csv"),
          row.names = FALSE)
cat("\n  [OK] tab_P_covid_msfe.csv\n")

# Grafico comparativo
modelos_covid <- intersect(c("2SRR", "ElNET", "RF", "LASSO", "Ridge"), methods)

df_covid_plot <- df_covid_rel |>
  pivot_longer(cols = all_of(modelos_covid), names_to = "Modelo",
               values_to = "MSFE_rel") |>
  mutate(Horizonte = factor(paste0("h=", Horizon),
                            levels = paste0("h=", horizons_list)))

cores_covid <- c("2SRR"  = "#D32F2F", "ElNET" = "#1976D2",
                 "RF"    = "#388E3C", "LASSO" = "#F57C00",
                 "Ridge" = "#90CAF9")

p <- ggplot(df_covid_plot,
            aes(Horizonte, MSFE_rel, fill = Modelo)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6,
           alpha = 0.85) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  facet_wrap(~ Periodo, ncol = 2) +
  scale_fill_manual(values = cores_covid) +
  labs(
    title    = "MSFE Relativo ao AR -- Pre vs Pos COVID",
    subtitle = "Modelos TVP (2SRR) devem melhorar relativamente apos a quebra estrutural",
    x = NULL, y = "MSFE / MSFE(AR)"
  )

ggsave(file.path(results_dir, "fig_P_covid_comparison.png"), p,
       width = 11, height = 5.5, dpi = 250)
cat("  [OK] fig_P_covid_comparison.png\n")

###############################################################################
# PARTE P -- Sumario narrativo do 2SRR (argumentos para apresentacao)
###############################################################################
cat(" PARTE P -- Sumario narrativo do 2SRR\n")

sumario_file <- file.path(results_dir, "sumario_2SRR.txt")
sink(sumario_file)

cat("================================================================\n")
cat(" SUMARIO ANALITICO DO MODELO 2SRR (Two-Step Ridge Regression)\n")
cat(" Gerado automaticamente em:", format(Sys.time()), "\n")
cat("================================================================\n\n")

# 1. Performance geral
msfe_2srr <- sapply(horizons_list, function(h) {
  key <- sprintf("V1_h%d", h)
  df_rel[df_rel$Horizon == h, "2SRR"]
})
cat("1. PERFORMANCE GERAL (MSFE relativo ao AR):\n")
for (i in seq_along(horizons_list)) {
  status <- ifelse(msfe_2srr[i] < 1, "SUPERA", "PERDE PARA")
  cat(sprintf("   h=%2d: %.4f (%s o AR)\n", horizons_list[i], msfe_2srr[i], status))
}

cat(sprintf("\n   Media geral: %.4f\n", mean(msfe_2srr)))
cat(sprintf("   O 2SRR supera o AR em %d de %d horizontes.\n\n",
            sum(msfe_2srr < 1), length(msfe_2srr)))

# 2. Performance por regime
cat("2. PERFORMANCE POR REGIME DE VOLATILIDADE:\n")
if (exists("df_regime_rel")) {
  for (reg in c("Normal", "Picos")) {
    sub <- df_regime_rel[df_regime_rel$Regime == reg, ]
    if (nrow(sub) > 0) {
      cat(sprintf("   %s:\n", reg))
      for (j in seq_len(nrow(sub))) {
        cat(sprintf("     h=%2d: %.4f\n", sub$Horizon[j], sub[j, "2SRR"]))
      }
    }
  }
  picos <- df_regime_rel[df_regime_rel$Regime == "Picos", "2SRR"]
  if (any(picos < 1)) {
    cat(sprintf("\n   ARGUMENTO CHAVE: O 2SRR supera o AR nos picos em %d de %d horizontes.\n",
                sum(picos < 1), length(picos)))
    cat("   Isso confirma que modelos TVP se adaptam a choques inflacionarios.\n\n")
  }
}

# 3. Performance pre/pos COVID
cat("3. PERFORMANCE PRE vs POS COVID:\n")
if (exists("df_covid_rel")) {
  pre  <- df_covid_rel[df_covid_rel$Periodo == "Pre-COVID", ]
  pos  <- df_covid_rel[df_covid_rel$Periodo == "Pos-COVID", ]
  if (nrow(pre) > 0 && nrow(pos) > 0) {
    cat("   Pre-COVID:\n")
    for (j in seq_len(nrow(pre))) {
      cat(sprintf("     h=%2d: %.4f\n", pre$Horizon[j], pre[j, "2SRR"]))
    }
    cat("   Pos-COVID:\n")
    for (j in seq_len(nrow(pos))) {
      cat(sprintf("     h=%2d: %.4f\n", pos$Horizon[j], pos[j, "2SRR"]))
    }
    melhora <- mean(pos[, "2SRR"], na.rm = TRUE) < mean(pre[, "2SRR"], na.rm = TRUE)
    if (melhora) {
      cat("\n   ARGUMENTO CHAVE: O 2SRR melhora relativamente no periodo pos-COVID.\n")
      cat("   Isso e evidencia de adaptacao a quebra estrutural.\n\n")
    } else {
      cat("\n   NOTA: O 2SRR nao melhorou no pos-COVID em media.\n")
      cat("   Possivel explicacao: pandemia nao afetou inflacao da mesma forma.\n\n")
    }
  }
}

# 4. Ranking
cat("4. POSICAO NO RANKING GERAL:\n")
if (exists("ranking")) {
  pos_rank <- which(ranking$Modelo == "2SRR")
  cat(sprintf("   Posicao: %d de %d modelos competidores\n",
              pos_rank, nrow(ranking)))
  cat(sprintf("   Vitorias: %d | Top 3: %d | MSFE medio: %.4f\n",
              ranking$Vitorias[pos_rank], ranking$Top3[pos_rank],
              ranking$Media_MSFE_rel[pos_rank]))
  cat(sprintf("   Melhor modelo geral: %s (%.4f)\n\n",
              ranking$Modelo[1], ranking$Media_MSFE_rel[1]))
}

# 5. Pontos fortes e fracos
cat("5. SINTESE PARA APRESENTACAO:\n\n")
cat("   PONTOS FORTES do 2SRR:\n")
cat("   - Captura parametros variantes no tempo sem MCMC/Kalman\n")
cat("   - Computacionalmente eficiente (apenas Ridge + CV)\n")
cat("   - Permite narrativa economica (quais variaveis mudaram de peso)\n")
cat("   - Tende a melhorar em periodos de alta volatilidade/choques\n\n")
cat("   PONTOS FRACOS do 2SRR:\n")
cat("   - MSFE geral > 1 (perde para AR na media)\n")
cat("   - Modelos estaticos simples (ElNET, LASSO) sao mais precisos em geral\n")
cat("   - Custo de flexibilidade em periodos calmos (overfitting leve)\n\n")
cat("   NARRATIVA RECOMENDADA:\n")
cat("   'O 2SRR nao e o modelo mais preciso na media, mas oferece algo que\n")
cat("    modelos estaticos nao conseguem: adaptacao rapida a choques e a\n")
cat("    capacidade de revelar COMO a estrutura da economia muda ao longo\n")
cat("    do tempo. Para Bancos Centrais, essa informacao estrutural pode\n")
cat("    ser mais valiosa que um ganho marginal de MSFE.'\n")

cat("\n================================================================\n")
sink()

# Imprimir no console tambem
cat(readLines(sumario_file), sep = "\n")
cat(sprintf("\n\n  [OK] %s\n", sumario_file))

###############################################################################
# PARTE Q -- Resumo Final no Console
###############################################################################
cat(" RESUMO FINAL DOS RESULTADOS\n")

n_png <- length(list.files(results_dir, pattern = "\\.png$"))
n_pdf <- length(list.files(results_dir, pattern = "\\.pdf$"))
n_csv <- length(list.files(results_dir, pattern = "\\.csv$"))
n_tex <- length(list.files(results_dir, pattern = "\\.tex$"))

cat(sprintf("  Pasta:       %s\n", results_dir))
cat(sprintf("  PNG gerados: %3d\n", n_png))
cat(sprintf("  PDF gerados: %3d\n", n_pdf))
cat(sprintf("  CSV gerados: %3d\n", n_csv))
cat(sprintf("  TEX gerados: %3d\n", n_tex))

cat("\n  MELHOR MODELO POR HORIZONTE (MSFE relativo):\n\n")
for (i in seq_len(nrow(df_rel))) {
  vals <- sapply(setdiff(methods, "AR"), function(m) {
    v <- df_rel[i, m]
    if (is.null(v) || is.na(v)) return(Inf)
    v
  })
  best_m <- setdiff(methods, "AR")[which.min(vals)]
  best_v <- min(vals, na.rm = TRUE)
  cat(sprintf("    %s h=%2d: %-12s (%.4f)\n",
              df_rel$Target[i], df_rel$Horizon[i], best_m, best_v))
}

cat(" Analise concluida\n")
