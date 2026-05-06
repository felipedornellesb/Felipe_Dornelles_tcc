###############################################################################
rm(list = ls())
set.seed(2024)

cat("
==============================================================
 11_results_analysis.R - v0.0.2
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
    # Gerar datas para os betas (mesmo tamanho que nrow)
    if (nrow(betas_obj) <= length(dates_eval)) {
      beta_dates <- tail(dates_eval, nrow(betas_obj))
    } else {
      beta_dates <- seq.Date(as.Date("2000-01-01"), by = "month",
                             length.out = nrow(betas_obj))
    }
    forecast_output$tvp_full <- list()
    forecast_output$tvp_full[[target_name]] <- list(
      betas = as.matrix(betas_obj),
      dates = beta_dates
    )
  } else {
    cat(")\n  [AVISO] betas_2SRR nao e matriz, pulando Parte G.\n")
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

# Paleta de cores: gerar automaticamente para N modelos
set.seed(42)
if (n_models <= 12) {
  paleta_base <- c("black", "grey55", "steelblue", "#D32F2F", "#2E7D32",
                   "#FF9800", "#9C27B0", "#00BCD4", "#795548", "#E91E63",
                   "#3F51B5", "#CDDC39")
} else {
  paleta_base <- c("black", "grey55",
                   scales::hue_pal()(n_models - 1))
}
cores_modelos <- setNames(
  c("black", paleta_base[2:(n_models + 1)]),
  c("Realizado", methods)
)

cat("=== Configuracao ===\n")
cat(sprintf("  Target:     %s\n", paste(targets, collapse = ", ")))
cat(sprintf("  Horizontes: %s\n", paste(horizons, collapse = ", ")))
cat(sprintf("  Metodos:    %s\n", paste(methods, collapse = ", ")))
cat(sprintf("  n_eval:     %d\n", n_eval))
cat(sprintf("  Chaves:     %d combinacoes\n\n", length(res)))

###############################################################################
# PARTE A -- Serie Real + Previsoes por Horizonte
###############################################################################
cat("=======================================================\n")
cat(" PARTE A -- Graficos: Serie Real vs Previsoes\n")

for (key in names(res)) {
  r   <- res[[key]]
  n_e <- length(r$actuals)
  if (n_e < 10) next

  # Montar dataframe
  df_plot <- data.frame(Date = r$dates[1:n_e], Realizado = r$actuals)
  for (m in methods) {
    if (m %in% colnames(r$preds)) {
      df_plot[[m]] <- r$preds[1:n_e, m]
    }
  }

  df_long <- df_plot |>
    pivot_longer(-Date, names_to = "Serie", values_to = "Valor") |>
    mutate(Serie = factor(Serie, levels = c("Realizado", methods)))

  # Linewidths
  lw_vals <- c(Realizado = 1.3)
  for (m in methods) lw_vals[m] <- 0.5

  p <- ggplot(df_long, aes(Date, Valor, color = Serie, linewidth = Serie)) +
    geom_line(na.rm = TRUE, alpha = 0.85) +
    scale_color_manual(values = cores_modelos) +
    scale_linewidth_manual(values = lw_vals, guide = "none") +
    labs(
      title    = sprintf("%s -- Horizonte h = %d", r$target, r$horizon),
      subtitle = "Serie realizada vs previsoes pseudo-out-of-sample",
      x = NULL, y = "Valor"
    )

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
cat("\n=======================================================\n")
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
cat("\n=======================================================\n")
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
cat("\n=======================================================\n")
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
cat("\n=======================================================\n")
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

    ok             <- complete.cases(loss_mat)
    loss_mat_clean <- loss_mat[ok, , drop = FALSE]
    if (nrow(loss_mat_clean) < 20) next

    tryCatch({
      # Remover cl=NULL que causa erro em algumas versoes do MCS
      mcs_out <- MCSprocedure(as.data.frame(loss_mat_clean),
                              alpha = 0.15, B = 5000,
                              statistic = "Tmax")

      # Extrair modelos sobreviventes
      show_df <- mcs_out@show
      if ("Rank_M" %in% colnames(show_df)) {
        surviving <- rownames(show_df)
      } else {
        surviving <- rownames(show_df)
      }

      cat(sprintf("  %s h=%d: MCS = {%s}\n",
                  r$target, r$horizon, paste(surviving, collapse = ", ")))

      mcs_results[[key]] <- data.frame(
        Target     = r$target,
        Horizon    = r$horizon,
        MCS_Models = paste(surviving, collapse = ", ")
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
cat("\n=======================================================\n")
cat(" PARTE G -- TVP Betas (2SRR Full Sample)\n")

if (!is.null(forecast_output$tvp_full)) {
  for (tgt in names(forecast_output$tvp_full)) {
    betas <- forecast_output$tvp_full[[tgt]]$betas
    d     <- forecast_output$tvp_full[[tgt]]$dates

    if (is.null(betas) || ncol(betas) < 2) {
      cat("  [AVISO] Betas insuficientes para gerar grafico.\n")
      next
    }

    beta_var <- apply(betas, 2, var)
    n_show   <- min(8, ncol(betas))
    top_vars <- names(sort(beta_var, decreasing = TRUE))[1:n_show]

    df_beta <- data.frame(Date = d, betas[, top_vars, drop = FALSE]) |>
      pivot_longer(-Date, names_to = "Variable", values_to = "Beta")

    p <- ggplot(df_beta, aes(Date, Beta, color = Variable)) +
      geom_line(linewidth = 0.6) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      facet_wrap(~ Variable, scales = "free_y", ncol = 2) +
      labs(
        title    = sprintf("TVP Betas -- Target: %s (Full Sample)", tgt),
        subtitle = "Coeficientes variantes no tempo (2SRR Coulombe)",
        x = NULL, y = expression(beta[t])
      ) +
      theme(legend.position = "none")

    fname <- sprintf("fig_G_tvp_betas_%s", tgt)
    ggsave(file.path(results_dir, paste0(fname, ".png")), p,
           width = 10, height = 10, dpi = 200)
    cat(sprintf("  [OK] %s.png\n", fname))
  }
} else {
  cat("  TVP Betas nao disponiveis.\n")
}

###############################################################################
# PARTE H -- Grafico MSFE Relativo (barras por horizonte)
###############################################################################
cat("\n=======================================================\n")
cat(" PARTE H -- Grafico MSFE Relativo (barras)\n")

# Selecionar modelos principais para visualizacao (senao fica poluido)
modelos_barra <- intersect(c("AR", "Ridge", "LASSO", "AdaLASSO", "ElNET",
                             "RF", "2SRR", "Factor", "CSR", "Bagging"),
                           methods)

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
    labs(
      title    = sprintf("MSFE Relativo ao AR -- %s", tgt),
      subtitle = "Abaixo de 1 = modelo supera o benchmark AR",
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
cat("\n=======================================================\n")
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
# PARTE J -- Resumo Final no Console
###############################################################################
cat("\n\n")
cat("==============================================================\n")
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
cat("--------------------------------------------------------------\n")

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

cat(" Analise concluida com sucesso.\n")
