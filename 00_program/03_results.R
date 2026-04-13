# ============================================================
# 03_results.R
#
# Avaliação POOS — estilo Coulombe (2SRR)
# Lê TVPfcst_V*_H*_M*.RData de 30_output/outputs_MM_DD_YYYY/
#
# Produz em 40_results/:
#   tab_rmsfe_relativo.csv
#   tab_rmsfe_dm_formatado.csv    <- tabela publicável (TCC)
#   tab_diebold_mariano.csv
#   tab_melhor_modelo.csv
#   fig_fcst_<VAR>_H<h>.png
#   fig_heatmap_H<h>.png
#   fig_rmsfe_barras_H<h>.png
# ============================================================

rm(list = ls())

# ============================================================
# 0. PATHS
# ============================================================

wd <- normalizePath(file.path(dirname(sys.frame(1)$ofile), ".."),
                    mustWork = FALSE)
if (!dir.exists(wd)) wd <- getwd()
setwd(wd)

paths <- list(
  data    = "10_data",
  output  = "30_output",
  results = "40_results"
)

if (!dir.exists(paths$results)) dir.create(paths$results, recursive = TRUE)

run_folder <- sort(
  list.dirs(paths$data, full.names = TRUE, recursive = FALSE)
) |> tail(1)

output_dir <- sort(
  list.dirs(paths$output, full.names = TRUE, recursive = FALSE)
) |> tail(1)

cat(sprintf("Dados de:   %s\n", run_folder))
cat(sprintf("Outputs de: %s\n", output_dir))

# ============================================================
# 1. PACOTES
# ============================================================

myPKGs <- c("dplyr", "tidyr", "ggplot2", "purrr",
            "stringr", "sandwich", "lmtest", "scales",
            "forcats", "ggtext")
need   <- myPKGs[!myPKGs %in% names(installed.packages()[, "Package"])]
if (length(need) > 0)
  install.packages(need, repos = "http://cran.us.r-project.org")
invisible(lapply(myPKGs, library, character.only = TRUE))

# ============================================================
# 2. CARREGA SÉRIES
# ============================================================

load(file.path(run_folder, "df_targets.rda"))
load(file.path(run_folder, "df_model.rda"))
load(file.path(run_folder, "targets_br.rda"))

target_names <- unname(unlist(targets_br))
szv          <- length(target_names)
dates_vec    <- df_model$date
bigt         <- length(dates_vec)

mat_targets <- df_targets |>
  dplyr::select(dplyr::all_of(target_names)) |>
  as.matrix()

cat(sprintf("T = %d  |  %s a %s\n", bigt,
            format(min(dates_vec), "%b/%Y"),
            format(max(dates_vec), "%b/%Y")))

# ============================================================
# 3. ARQUIVOS DE FORECAST
# ============================================================

rdata_files <- list.files(
  output_dir,
  pattern    = "^TVPfcst_V\\d+_H\\d+_M\\d+\\.RData$",
  full.names = TRUE
)

if (length(rdata_files) == 0)
  stop(sprintf("Nenhum arquivo TVPfcst encontrado em: %s", output_dir))

cat(sprintf("%d arquivo(s) encontrado(s).\n", length(rdata_files)))

file_meta <- tibble::tibble(path = rdata_files) |>
  dplyr::mutate(
    fname  = basename(path),
    V      = as.integer(stringr::str_extract(fname, "(?<=_V)\\d+")),
    H      = as.integer(stringr::str_extract(fname, "(?<=_H)\\d+")),
    M_spec = as.integer(stringr::str_extract(fname, "(?<=_M)\\d+"))
  )

horizons <- sort(unique(file_meta$H))
vars_idx <- sort(unique(file_meta$V))
specs    <- sort(unique(file_meta$M_spec))

# ============================================================
# 4. RÓTULOS
# ============================================================

est_labels <- c("1" = "Ridge", "2" = "2SRR", "3" = "MSRRs", "4" = "MSRRd")
spec_labels <- c("1" = "AR", "2" = "AR+F", "3" = "AR+Tgt", "4" = "AR+Pan")

# ============================================================
# 5. CARREGA TODOS OS FORECASTS → data.frame longo
# ============================================================

load_one <- function(meta_row) {
  e <- new.env(parent = emptyenv())
  load(meta_row$path, envir = e)
  fc_arr <- e$forecast
  v      <- meta_row$V
  h      <- meta_row$H
  m      <- meta_row$M_spec
  n_est  <- dim(fc_arr)[4]

  purrr::map_dfr(seq_len(n_est), function(est) {
    tibble::tibble(
      t         = seq_len(bigt),
      date      = dates_vec,
      V         = v, H = h, M_spec = m,
      estimador = est,
      forecast  = fc_arr[, h, v, est],
      realizado = mat_targets[, v]
    )
  })
}

df_long <- purrr::map_dfr(purrr::transpose(file_meta), load_one) |>
  dplyr::mutate(
    var_name   = target_names[V],
    est_label  = est_labels[as.character(estimador)],
    spec_label = spec_labels[as.character(M_spec)]
  )

cat(sprintf("df_long: %d linhas\n", nrow(df_long)))

# ============================================================
# 6. ERROS
# ============================================================

df_erros          <- df_long
df_erros$erro     <- df_erros$realizado - df_erros$forecast
df_erros$se       <- df_erros$erro^2
df_erros$ae       <- abs(df_erros$erro)
df_erros          <- df_erros[!is.na(df_erros$erro), ]

# ============================================================
# 7. MSFE e RMSFE RELATIVO
# ============================================================

df_msfe <- df_erros |>
  dplyr::group_by(var_name, V, H, M_spec, estimador, est_label, spec_label) |>
  dplyr::summarise(n_obs = dplyr::n(), MSFE = mean(se), MAE = mean(ae),
                   .groups = "drop")

benchmark <- df_msfe |>
  dplyr::filter(estimador == 1, M_spec == 1) |>
  dplyr::transmute(var_name, V, H, MSFE_bench = MSFE, n_bench = n_obs)

df_rmsfe <- df_msfe |>
  dplyr::left_join(benchmark, by = c("var_name", "V", "H")) |>
  dplyr::mutate(RMSFE_rel = sqrt(MSFE / MSFE_bench))

# ============================================================
# 8. TABELA RMSFE
# ============================================================

tab_main <- df_rmsfe |>
  dplyr::filter(!(estimador == 1 & M_spec == 1)) |>
  dplyr::mutate(col_label = paste0(spec_label, "_", est_label)) |>
  dplyr::select(var_name, H, col_label, RMSFE_rel) |>
  tidyr::pivot_wider(names_from = col_label, values_from = RMSFE_rel) |>
  dplyr::arrange(var_name, H)

cat("\n=== RMSFE relativo ao Ridge-AR ===\n")
print(tab_main, n = 100, width = 200)
write.csv(tab_main, file.path(paths$results, "tab_rmsfe_relativo.csv"),
          row.names = FALSE)

# ============================================================
# 9. DIEBOLD-MARIANO (HAC, Newey-West)
# ============================================================

dm_test_safe <- function(e_model, e_bench) {
  d <- e_model^2 - e_bench^2
  d <- d[!is.na(d)]
  n <- length(d)
  if (n < 8)
    return(data.frame(DM_stat = NA_real_, p_value = NA_real_,
                      DM_pval_onesided = NA_real_))
  tryCatch({
    fit <- lm(d ~ 1)
    ct  <- lmtest::coeftest(fit,
                            vcov = sandwich::NeweyWest(fit, lag = 4,
                                                      prewhite = FALSE))
    dm_s <- ct[1, "t value"]
    pv2s <- ct[1, "Pr(>|t|)"]
    pv1s <- pt(dm_s, df = n - 1, lower.tail = TRUE)
    data.frame(DM_stat = dm_s, p_value = pv2s, DM_pval_onesided = pv1s)
  }, error = function(e)
    data.frame(DM_stat = NA_real_, p_value = NA_real_,
               DM_pval_onesided = NA_real_)
  )
}

bench_erros <- df_erros[df_erros$estimador == 1 & df_erros$M_spec == 1,
                        c("var_name", "V", "H", "t", "erro")]
names(bench_erros)[names(bench_erros) == "erro"] <- "erro_bench"

df_erros_dm <- merge(
  df_erros[!(df_erros$estimador == 1 & df_erros$M_spec == 1), ],
  bench_erros, by = c("var_name", "V", "H", "t")
)
df_erros_dm <- df_erros_dm[!is.na(df_erros_dm$erro_bench), ]

grupos   <- unique(df_erros_dm[, c("var_name", "H", "M_spec",
                                    "estimador", "est_label", "spec_label")])
dm_list  <- vector("list", nrow(grupos))

for (i in seq_len(nrow(grupos))) {
  g   <- grupos[i, ]
  sub <- df_erros_dm[
    df_erros_dm$var_name  == g$var_name  &
    df_erros_dm$H         == g$H         &
    df_erros_dm$M_spec    == g$M_spec    &
    df_erros_dm$estimador == g$estimador, ]
  dm_list[[i]] <- cbind(g, dm_test_safe(sub$erro, sub$erro_bench))
}

df_dm <- do.call(rbind, dm_list)
rownames(df_dm) <- NULL
df_dm$sig_2s <- ifelse(df_dm$p_value < 0.01, "***",
               ifelse(df_dm$p_value < 0.05, "**",
               ifelse(df_dm$p_value < 0.10, "*", "")))
df_dm$sig_1s <- ifelse(df_dm$DM_pval_onesided < 0.01, "†††",
               ifelse(df_dm$DM_pval_onesided < 0.05, "††",
               ifelse(df_dm$DM_pval_onesided < 0.10, "†", "")))

write.csv(df_dm,
          file.path(paths$results, "tab_diebold_mariano.csv"),
          row.names = FALSE)

# ============================================================
# 10. TABELA PUBLICÁVEL — RMSFE + asteriscos DM
# ============================================================

tab_pub_base <- merge(
  df_rmsfe[!(df_rmsfe$estimador == 1 & df_rmsfe$M_spec == 1), ],
  df_dm[, c("var_name", "H", "M_spec", "estimador", "sig_2s")],
  by = c("var_name", "H", "M_spec", "estimador"), all.x = TRUE
)
tab_pub_base$sig_2s[is.na(tab_pub_base$sig_2s)] <- ""
tab_pub_base$col_label <- paste0(tab_pub_base$spec_label, "_", tab_pub_base$est_label)
tab_pub_base$cell      <- sprintf("%.3f%s", tab_pub_base$RMSFE_rel, tab_pub_base$sig_2s)

tab_pub <- tab_pub_base |>
  dplyr::select(var_name, H, col_label, cell) |>
  tidyr::pivot_wider(names_from = col_label, values_from = cell) |>
  dplyr::arrange(var_name, H)

write.csv(tab_pub,
          file.path(paths$results, "tab_rmsfe_dm_formatado.csv"),
          row.names = FALSE)
cat("\n=== Tabela publicável salva ===\n")
print(tab_pub, n = 100, width = 200)

# ============================================================
# 11. MELHOR MODELO
# ============================================================

best_models <- df_rmsfe[!(df_rmsfe$estimador == 1 & df_rmsfe$M_spec == 1), ] |>
  dplyr::group_by(var_name, H) |>
  dplyr::slice_min(RMSFE_rel, n = 1, with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::select(var_name, H, spec_label, est_label, RMSFE_rel, n_obs)

write.csv(best_models,
          file.path(paths$results, "tab_melhor_modelo.csv"),
          row.names = FALSE)
cat("\n=== Melhor modelo ===\n")
print(best_models)

# ============================================================
# 12. FIGURAS — previsão vs. realizado
# ============================================================

df_plot_best <- merge(
  df_erros,
  best_models[, c("var_name", "H", "spec_label", "est_label")],
  by = c("var_name", "H", "spec_label", "est_label")
)

for (vname in target_names) {
  for (h_val in horizons) {
    df_v <- df_plot_best[df_plot_best$var_name == vname &
                           df_plot_best$H == h_val, ]
    df_v <- df_v[order(df_v$date), ]
    if (nrow(df_v) < 5) next

    rmsfe_txt <- sprintf("%.3f",
      best_models$RMSFE_rel[best_models$var_name == vname &
                               best_models$H == h_val])

    p <- ggplot2::ggplot(df_v, ggplot2::aes(x = date)) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.4) +
      ggplot2::geom_line(ggplot2::aes(y = realizado, colour = "Realizado"),
                         linewidth = 0.75) +
      ggplot2::geom_line(ggplot2::aes(y = forecast, colour = "Previsão"),
                         linewidth = 0.75, linetype = "dashed") +
      ggplot2::scale_colour_manual(
        values = c("Realizado" = "#01696f", "Previsão" = "#964219"), name = NULL
      ) +
      ggplot2::scale_x_date(date_labels = "%Y", date_breaks = "3 years") +
      ggplot2::labs(
        title    = sprintf("%s — H=%d | %s/%s",
                           vname, h_val,
                           unique(df_v$spec_label), unique(df_v$est_label)),
        subtitle = sprintf("RMSFE relativo = %s", rmsfe_txt),
        x = NULL, y = "Variação"
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(legend.position = "top",
                     panel.grid.minor = ggplot2::element_blank())

    ggplot2::ggsave(
      file.path(paths$results, sprintf("fig_fcst_%s_H%d.png", vname, h_val)),
      p, width = 9, height = 4, dpi = 150
    )
  }
}

# ============================================================
# 13. FIGURA — heatmap RMSFE
# ============================================================

for (h_val in horizons) {
  df_heat <- df_rmsfe[df_rmsfe$H == h_val, ]
  df_heat$modelo_full <- paste0(df_heat$spec_label, "\n", df_heat$est_label)
  df_heat$var_name    <- factor(df_heat$var_name, levels = rev(target_names))
  col_order <- unique(df_heat[order(df_heat$M_spec, df_heat$estimador),
                               "modelo_full"])
  df_heat$modelo_full <- factor(df_heat$modelo_full, levels = col_order)

  p_heat <- ggplot2::ggplot(
    df_heat, ggplot2::aes(x = modelo_full, y = var_name, fill = RMSFE_rel)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.3f", RMSFE_rel),
                   colour = RMSFE_rel > 1.05),
      size = 2.8, fontface = "bold", show.legend = FALSE
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#01696f", mid = "#f7f6f2", high = "#a12c7b",
      midpoint = 1, limits = c(0.6, 1.4),
      oob = scales::squish, name = "RMSFE\nrelativo"
    ) +
    ggplot2::scale_colour_manual(
      values = c("FALSE" = "grey20", "TRUE" = "grey60")
    ) +
    ggplot2::labs(
      title    = sprintf("RMSFE relativo ao Ridge-AR — H = %d", h_val),
      subtitle = "< 1 = melhor que benchmark | * p<0.10 ** p<0.05 *** p<0.01 (DM bilateral)",
      x = NULL, y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid  = ggplot2::element_blank(),
      plot.title  = ggplot2::element_text(face = "bold"),
      axis.text.x = ggplot2::element_text(size = 8, lineheight = 1.2)
    )

  ggplot2::ggsave(
    file.path(paths$results, sprintf("fig_heatmap_H%d.png", h_val)),
    p_heat, width = 12, height = 5, dpi = 150
  )
}

# ============================================================
# 14. FIGURA — barras RMSFE por variável
# ============================================================

for (h_val in horizons) {
  df_bar <- df_rmsfe[df_rmsfe$H == h_val &
                       !(df_rmsfe$estimador == 1 & df_rmsfe$M_spec == 1), ]
  df_bar$modelo_full <- paste0(df_bar$spec_label, " / ", df_bar$est_label)
  df_bar$var_name    <- factor(df_bar$var_name, levels = target_names)
  df_bar$melhor      <- df_bar$RMSFE_rel < 1

  p_bar <- ggplot2::ggplot(
    df_bar, ggplot2::aes(x = reorder(modelo_full, RMSFE_rel),
                         y = RMSFE_rel, fill = melhor)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed",
                        colour = "grey40", linewidth = 0.5) +
    ggplot2::scale_fill_manual(
      values = c("TRUE" = "#01696f", "FALSE" = "#964219"), guide = "none"
    ) +
    ggplot2::scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
    ggplot2::facet_wrap(~var_name, scales = "free_y", ncol = 2) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title    = sprintf("RMSFE relativo ao Ridge-AR — H = %d", h_val),
      subtitle = "Verde = melhor que benchmark | Laranja = pior",
      x = NULL, y = "RMSFE relativo"
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      strip.text         = ggplot2::element_text(face = "bold")
    )

  ggplot2::ggsave(
    file.path(paths$results, sprintf("fig_rmsfe_barras_H%d.png", h_val)),
    p_bar, width = 11, height = 3 * ceiling(szv / 2), dpi = 150
  )
}

# ============================================================
# 15. SUMÁRIO
# ============================================================

cat("\n============================================================\n")
cat(sprintf("Resultados em: %s\n", paths$results))
cat("  tab_rmsfe_relativo.csv\n")
cat("  tab_rmsfe_dm_formatado.csv   <- TCC\n")
cat("  tab_diebold_mariano.csv\n")
cat("  tab_melhor_modelo.csv\n")
cat("  fig_fcst_<VAR>_H<h>.png\n")
cat("  fig_heatmap_H<h>.png\n")
cat("  fig_rmsfe_barras_H<h>.png\n")
cat("============================================================\n")
