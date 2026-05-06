rm(list = ls())
set.seed(2024)

cat("
==============================================================
 10_results_analysis.R — Analise Completa de Resultados
==============================================================
\n")

# ---- 0. Pacotes ----
pkgs <- c("ggplot2", "dplyr", "tidyr", "scales", "gridExtra",
          "forecast", "viridis", "reshape2", "patchwork")
need <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(need)) install.packages(need, repos = "https://cran.r-project.org")
invisible(lapply(pkgs, library, character.only = TRUE))

# MCS
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
results_dir <- file.path("40_results", paste0("run_", timestamp_str))
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("Pasta de resultados: %s\n\n", results_dir))

# ---- 2. Carregar dados de forecast ----
# Procurar o RDS de resultados mais recente
cat("=== Carregando resultados de forecast ===\n")

# Tentar carregar forecast_results.rda / .rds
found_results <- FALSE

# Procurar em varios locais possiveis
search_paths <- c(
  "30_output/forecast_results.rda",
  "30_output/forecast_results.rds",
  "40_results/forecast_results.rda",
  "40_results/forecast_results.rds",
  "forecast_results.rda",
  "forecast_results.rds"
)

# Tambem procurar em subpastas de 10_data/
data_dirs <- list.dirs("10_data", recursive = TRUE, full.names = TRUE)
for (d in data_dirs) {
  search_paths <- c(search_paths,
                    file.path(d, "forecast_results.rda"),
                    file.path(d, "forecast_results.rds"))
}

for (sp in search_paths) {
  if (file.exists(sp)) {
    cat(sprintf("  Carregando: %s\n", sp))
    tryCatch({
      if (grepl("\\.rda$", sp)) {
        load(sp)
      } else {
        forecast_output <- readRDS(sp)
      }
      found_results <- TRUE
      cat("  [OK] Resultados carregados!\n")
      break
    }, error = function(e) {
      cat(sprintf("  Erro: %s\n", e$message))
    })
  }
}

# Se nao encontrou resultados salvos, procurar CSVs de previsoes individuais
if (!found_results) {
  cat("\n  Procurando CSVs de resultados individuais...\n")

  # Procurar em todas as pastas
  all_rds <- list.files(".", pattern = "\\.(rds|rda|RDS|RDA)$",
                        recursive = TRUE, full.names = TRUE)
  all_csv <- list.files(".", pattern = "\\.csv$",
                        recursive = TRUE, full.names = TRUE)

  cat(sprintf("  Encontrados: %d RDS/RDA, %d CSV\n",
              length(all_rds), length(all_csv)))

  # Listar todos para debug
  if (length(all_rds) > 0) {
    cat("  RDS/RDA encontrados:\n")
    for (f in all_rds) cat(sprintf("    %s (%s bytes)\n", f, file.size(f)))
  }
  if (length(all_csv) > 0) {
    cat("  CSV encontrados:\n")
    for (f in all_csv) cat(sprintf("    %s (%s bytes)\n", f, file.size(f)))
  }
}

# ---- 3. Se nao encontrou resultados, gerar dados de demonstracao ----
# (isso permite testar o script mesmo sem ter rodado o forecast)
if (!found_results || !exists("forecast_output")) {
  cat("\n[AVISO] Resultados de forecast nao encontrados.\n")
  cat("  Gerando dados de demonstracao para validar o pipeline...\n")
  cat("  (Rode o forecast.R primeiro para obter resultados reais)\n\n")

  # Dados simulados realistas
  set.seed(2024)
  n_eval <- 100
  dates_eval <- seq.Date(as.Date("2015-01-01"), by = "month", length.out = n_eval)

  # Simular serie real + previsoes de 4 modelos para 4 targets e 4 horizontes
  targets_list <- c("GDP", "INFLATION", "IR", "SPREAD")
  horizons_list <- c(1, 3, 6, 12)
  methods_list <- c("AR", "Ridge", "TVP-Ridge", "RF")

  forecast_output <- list(
    results = list(),
    targets_br = as.list(setNames(targets_list, paste0("V", 1:4))),
    horizons = horizons_list,
    methods = methods_list,
    dates = dates_eval
  )

  for (vi in seq_along(targets_list)) {
    tgt <- targets_list[vi]
    for (h in horizons_list) {
      key <- sprintf("V%d_h%d", vi, h)

      # Serie real
      actual <- cumsum(rnorm(n_eval, 0, 0.5)) +
        sin(seq(0, 4 * pi, length.out = n_eval))

      # Previsoes com diferentes niveis de erro
      preds_mat <- matrix(NA, n_eval, length(methods_list))
      colnames(preds_mat) <- methods_list
      noise_levels <- c(AR = 1.0, Ridge = 0.8, `TVP-Ridge` = 0.6, RF = 0.7)

      for (m in methods_list) {
        noise <- noise_levels[m] * (1 + h / 12)
        preds_mat[, m] <- actual + rnorm(n_eval, 0, noise)
      }

      forecast_output$results[[key]] <- list(
        target = tgt,
        horizon = h,
        dates = dates_eval,
        actuals = actual,
        preds = preds_mat
      )
    }
  }
  cat("  [OK] Dados de demonstracao gerados.\n\n")
}

# ---- Extrair objetos ----
res      <- forecast_output$results
horizons <- forecast_output$horizons
methods  <- forecast_output$methods
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
            theme(plot.title    = element_text(face = "bold", size = 14),
                  plot.subtitle = element_text(color = "grey40", size = 10),
                  legend.position = "bottom",
                  legend.title    = element_blank(),
                  panel.grid.minor = element_blank(),
                  strip.text = element_text(face = "bold")))

cores_modelos <- c("Realizado" = "black",
                   "AR"        = "grey55",
                   "Ridge"     = "steelblue",
                   "TVP-Ridge" = "#D32F2F",
                   "RF"        = "#2E7D32")

cat("=== Configuracao ===\n")
cat(sprintf("  Targets:    %s\n", paste(targets, collapse = ", ")))
cat(sprintf("  Horizontes: %s\n", paste(horizons, collapse = ", ")))
cat(sprintf("  Metodos:    %s\n", paste(methods, collapse = ", ")))
cat(sprintf("  Chaves:     %d combinacoes\n\n", length(res)))

###############################################################################
# PARTE A — Serie Real + Previsoes por Horizonte
###############################################################################
cat("=======================================================\n")
cat(" PARTE A — Graficos: Serie Real vs Previsoes\n")
cat("=======================================================\n\n")

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
  for (m in methods) lw_vals[m] <- 0.6

  p <- ggplot(df_long, aes(Date, Valor, color = Serie, linewidth = Serie)) +
    geom_line(na.rm = TRUE, alpha = 0.85) +
    scale_color_manual(values = cores_modelos) +
    scale_linewidth_manual(values = lw_vals, guide = "none") +
    labs(
      title    = sprintf("%s — Horizonte h = %d", r$target, r$horizon),
      subtitle = "Serie realizada vs previsoes pseudo-out-of-sample",
      x = NULL, y = "Valor"
    )

  fname <- sprintf("fig_A_forecast_%s_h%02d", r$target, r$horizon)
  ggsave(file.path(results_dir, paste0(fname, ".png")), p,
         width = 11, height = 5.5, dpi = 250)
  ggsave(file.path(results_dir, paste0(fname, ".pdf")), p,
         width = 11, height = 5.5)

  cat(sprintf("  [OK] %s.png/pdf\n", fname))
}

###############################################################################
# PARTE B — Tabela MSFE Relativo ao AR (benchmark)
###############################################################################
cat("\n=======================================================\n")
cat(" PARTE B — Tabelas de Metricas (MSFE / RMSE)\n")
cat("=======================================================\n\n")

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

# Salvar
write.csv(df_rel,  file.path(results_dir, "tab_B_msfe_relative.csv"),
          row.names = FALSE)
write.csv(df_rmse, file.path(results_dir, "tab_B_rmse_absolute.csv"),
          row.names = FALSE)
cat("\n  [OK] tab_B_msfe_relative.csv\n  [OK] tab_B_rmse_absolute.csv\n")

###############################################################################
# PARTE C — Teste Diebold-Mariano
###############################################################################
cat("\n=======================================================\n")
cat(" PARTE C — Teste Diebold-Mariano (DM)\n")
cat("=======================================================\n\n")

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
        pval <- dm$p.value
        stars <- ifelse(pval < 0.01, "***",
                        ifelse(pval < 0.05, "**",
                               ifelse(pval < 0.10, "*", "")))
        row[[paste0("DM_pval_", m)]]  <- round(pval, 4)
        row[[paste0("DM_stars_", m)]] <- stars
      }, error = function(e) {
        row[[paste0("DM_pval_", m)]]  <<- NA
        row[[paste0("DM_stars_", m)]] <<- ""
      })
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
# PARTE D — CSFE (Cumulative Squared Forecast Errors)
###############################################################################
cat("\n=======================================================\n")
cat(" PARTE D — CSFE (Cumulative Squared Forecast Errors)\n")
cat("=======================================================\n\n")

for (key in names(res)) {
  r   <- res[[key]]
  n_e <- length(r$actuals)
  if (n_e < 20) next

  e2_ar   <- (r$preds[1:n_e, "AR"] - r$actuals)^2
  df_csfe <- data.frame(Date = r$dates[1:n_e])

  for (m in methods) {
    if (m == "AR") next
    e2_m <- (r$preds[1:n_e, m] - r$actuals)^2
    # CSFE = cumsum(e2_AR - e2_modelo) -> sobe quando modelo ganha
    csfe <- cumsum(ifelse(is.na(e2_ar) | is.na(e2_m), 0, e2_ar - e2_m))
    df_csfe[[m]] <- csfe
  }

  df_csfe_long <- df_csfe |>
    pivot_longer(-Date, names_to = "Modelo", values_to = "CSFE")

  p <- ggplot(df_csfe_long, aes(Date, CSFE, color = Modelo)) +
    geom_line(linewidth = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = cores_modelos[-1]) +
    labs(
      title    = sprintf("CSFE — %s, h = %d", r$target, r$horizon),
      subtitle = "Acima de 0 = modelo supera o AR; abaixo = AR e melhor",
      x = NULL, y = "CSFE (acumulado)"
    )

  fname <- sprintf("fig_D_csfe_%s_h%02d", r$target, r$horizon)
  ggsave(file.path(results_dir, paste0(fname, ".png")), p,
         width = 10, height = 5, dpi = 250)
  cat(sprintf("  [OK] %s.png\n", fname))
}

###############################################################################
# PARTE E — Analise de Residuos (histograma + QQ-plot + ACF)
###############################################################################
cat("\n=======================================================\n")
cat(" PARTE E — Analise de Residuos\n")
cat("=======================================================\n\n")

for (key in names(res)) {
  r   <- res[[key]]
  n_e <- length(r$actuals)
  if (n_e < 20) next

  for (m in methods) {
    residuos <- r$actuals - r$preds[1:n_e, m]
    residuos <- residuos[!is.na(residuos)]
    if (length(residuos) < 20) next

    # Histograma
    df_res <- data.frame(Residuo = residuos)
    p_hist <- ggplot(df_res, aes(Residuo)) +
      geom_histogram(bins = 25, fill = "steelblue", alpha = 0.7,
                     color = "white") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
      labs(
        title = sprintf("Residuos: %s — %s (h=%d)", m, r$target, r$horizon),
        subtitle = sprintf("Media=%.4f, SD=%.4f, Skew=%.3f",
                           mean(residuos), sd(residuos),
                           mean(((residuos - mean(residuos)) /
                                   sd(residuos))^3)),
        x = "Residuo (Real - Previsto)", y = "Frequencia"
      )

    # QQ-plot
    p_qq <- ggplot(df_res, aes(sample = Residuo)) +
      stat_qq(color = "steelblue", alpha = 0.6) +
      stat_qq_line(color = "red", linewidth = 0.8) +
      labs(title = "QQ-Plot vs Normal",
           x = "Quantis Teoricos", y = "Quantis Amostrais")

    # ACF manual
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

    # Combinar
    p_combined <- p_hist + p_qq + p_acf + plot_layout(ncol = 3) +
      plot_annotation(
        title = sprintf("Diagnostico de Residuos: %s — %s (h=%d)",
                        m, r$target, r$horizon)
      )

    fname <- sprintf("fig_E_residuos_%s_%s_h%02d", r$target, m, r$horizon)
    ggsave(file.path(results_dir, paste0(fname, ".png")), p_combined,
           width = 15, height = 5, dpi = 200)
  }
  cat(sprintf("  [OK] Residuos %s h=%d (todos os modelos)\n",
              r$target, r$horizon))
}

###############################################################################
# PARTE F — Model Confidence Set (MCS)
###############################################################################
cat("\n=======================================================\n")
cat(" PARTE F — Model Confidence Set (MCS)\n")
cat("=======================================================\n\n")

if (has_mcs) {
  mcs_results <- list()
  for (key in names(res)) {
    r   <- res[[key]]
    n_e <- length(r$actuals)
    if (n_e < 20) next

    # Montar matriz de perdas (squared errors)
    loss_mat <- matrix(NA, n_e, length(methods))
    colnames(loss_mat) <- methods
    for (m in methods) {
      loss_mat[, m] <- (r$preds[1:n_e, m] - r$actuals)^2
    }

    # Remover linhas com NA
    ok             <- complete.cases(loss_mat)
    loss_mat_clean <- loss_mat[ok, , drop = FALSE]
    if (nrow(loss_mat_clean) < 20) next

    tryCatch({
      mcs_out   <- MCSprocedure(as.data.frame(loss_mat_clean),
                                alpha = 0.15, B = 5000,
                                statistic = "Tmax", cl = NULL)
      surviving <- names(which(
        mcs_out@show[, "Rank_M"] <= nrow(mcs_out@show)
      ))
      cat(sprintf("  %s h=%d: MCS = {%s}\n",
                  r$target, r$horizon, paste(surviving, collapse = ", ")))
      mcs_results[[key]] <- data.frame(
        Target     = r$target,
        Horizon    = r$horizon,
        MCS_Models = paste(surviving, collapse = ", ")
      )
    }, error = function(e) {
      cat(sprintf("  %s h=%d: MCS erro — %s\n",
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
# PARTE G — TVP Betas ao Longo do Tempo (se disponiveis)
###############################################################################
cat("\n=======================================================\n")
cat(" PARTE G — TVP Betas (Full Sample)\n")
cat("=======================================================\n\n")

if (exists("forecast_output") && !is.null(forecast_output$tvp_full)) {
  for (tgt in names(forecast_output$tvp_full)) {
    betas <- forecast_output$tvp_full[[tgt]]$betas
    d     <- forecast_output$tvp_full[[tgt]]$dates

    beta_var <- apply(betas, 2, var)
    top_vars <- names(sort(beta_var, decreasing = TRUE))[
      1:min(8, ncol(betas))
    ]

    df_beta <- data.frame(Date = d, betas[, top_vars, drop = FALSE]) |>
      pivot_longer(-Date, names_to = "Variable", values_to = "Beta")

    p <- ggplot(df_beta, aes(Date, Beta, color = Variable)) +
      geom_line(linewidth = 0.6) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      facet_wrap(~ Variable, scales = "free_y", ncol = 2) +
      labs(
        title    = sprintf("TVP Betas — Target: %s (Full Sample)", tgt),
        subtitle = "Coeficientes variantes no tempo (Coulombe TVP-2SRR)",
        x = NULL, y = expression(beta[t])
      ) +
      theme(legend.position = "none")

    fname <- sprintf("fig_G_tvp_betas_%s", tgt)
    ggsave(file.path(results_dir, paste0(fname, ".png")), p,
           width = 10, height = 10, dpi = 200)
    cat(sprintf("  [OK] %s.png\n", fname))
  }
} else {
  cat("  TVP Betas nao disponiveis (rode o forecast com TVP-Ridge primeiro)\n")
}

###############################################################################
# PARTE H — Grafico MSFE Relativo (barras por horizonte)
###############################################################################
cat("\n=======================================================\n")
cat(" PARTE H — Grafico MSFE Relativo (barras)\n")
cat("=======================================================\n\n")

df_rel_long <- df_rel |>
  pivot_longer(cols = all_of(methods), names_to = "Modelo",
               values_to = "MSFE_rel") |>
  filter(Modelo != "AR") |>
  mutate(Horizonte = factor(paste0("h=", Horizon)))

for (tgt in unique(df_rel_long$Target)) {
  df_tgt <- df_rel_long |> filter(Target == tgt)

  p <- ggplot(df_tgt, aes(Horizonte, MSFE_rel, fill = Modelo)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6,
             alpha = 0.85) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red",
               linewidth = 0.7) +
    scale_fill_manual(values = cores_modelos[-c(1, 2)]) +
    labs(
      title    = sprintf("MSFE Relativo ao AR — %s", tgt),
      subtitle = "Abaixo de 1 = modelo supera o benchmark AR",
      x = NULL, y = "MSFE / MSFE(AR)"
    ) +
    coord_cartesian(
      ylim = c(0, max(df_tgt$MSFE_rel, na.rm = TRUE) * 1.1)
    )

  fname <- sprintf("fig_H_msfe_bars_%s", tgt)
  ggsave(file.path(results_dir, paste0(fname, ".png")), p,
         width = 8, height = 5, dpi = 250)
  cat(sprintf("  [OK] %s.png\n", fname))
}

###############################################################################
# PARTE I — Tabela LaTeX-Ready
###############################################################################
cat("\n=======================================================\n")
cat(" PARTE I — Tabela LaTeX\n")
cat("=======================================================\n\n")

# Gerar tabela LaTeX
latex_file <- file.path(results_dir, "tab_I_latex_msfe.tex")
sink(latex_file)

cat("\\begin{table}[ht]\n")
cat("\\centering\n")
cat("\\caption{MSFE relativo ao AR --- Dados Americanos}\n")
cat("\\label{tab:msfe_relative}\n")

# Colunas
methods_no_ar <- methods[methods != "AR"]
cat(sprintf("\\begin{tabular}{ll%s}\n",
            paste(rep("c", length(methods)), collapse = "")))
cat("\\hline\\hline\n")
cat(sprintf("Target & h & %s \\\\\n", paste(methods, collapse = " & ")))
cat("\\hline\n")

for (i in seq_len(nrow(df_rel))) {
  vals <- sapply(methods, function(m) {
    v <- df_rel[i, m]
    if (is.na(v)) return("---")

    # Adicionar negrito se for o melhor (menor) para esse target/horizonte
    all_vals <- sapply(methods[-1], function(mm) df_rel[i, mm])
    is_best  <- (!is.na(v) && v == min(all_vals, na.rm = TRUE) && m != "AR")
    formatted <- sprintf("%.3f", v)
    if (is_best) formatted <- sprintf("\\textbf{%s}", formatted)

    # Adicionar estrelas do DM (se disponivel)
    dm_col <- paste0("DM_stars_", m)
    if (dm_col %in% names(df_dm)) {
      stars <- df_dm[df_dm$Target == df_rel$Target[i] &
                       df_dm$Horizon == df_rel$Horizon[i], dm_col]
      if (length(stars) > 0 && !is.na(stars) && nchar(stars) > 0) {
        formatted <- paste0(formatted, "$^{", stars, "}$")
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
# PARTE J — Resumo Final no Console
###############################################################################
cat("\n\n")
cat("==============================================================\n")
cat(" RESUMO FINAL DOS RESULTADOS\n")
cat("==============================================================\n")

# Contar arquivos gerados
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

# Resumo dos resultados
cat("\n  MELHOR MODELO POR TARGET x HORIZONTE (MSFE relativo):\n\n")
for (i in seq_len(nrow(df_rel))) {
  vals <- sapply(methods[-1], function(m) {
    v <- df_rel[i, m]
    if (is.na(v)) return(Inf)
    v
  })
  best_m <- methods[-1][which.min(vals)]
  best_v <- min(vals, na.rm = TRUE)
  cat(sprintf("    %s h=%2d: %s (%.3f)\n",
              df_rel$Target[i], df_rel$Horizon[i], best_m, best_v))
}

cat("\n==============================================================\n")
cat(" Analise concluida com sucesso.\n")
cat("==============================================================\n")
