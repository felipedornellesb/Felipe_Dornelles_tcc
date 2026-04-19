# ============================================================
# 03_03_results.R
#
# Avaliação POOS completa — estilo Coulombe (2SRR)
# Lê todos os TVPfcst_V*_H*_M*.RData já gerados
#
# Correções aplicadas:
#   1) df_erros construído via $ em vez de pipe/mutate
#   2) filter sem referência interna a colunas novas
#   3) todas as dependências em ordem sequencial
#
# Produz em 40_results/:
#   tab_rmsfe_relativo.csv
#   tab_rmsfe_dm_formatado.csv   <- tabela publicável (TCC)
#   tab_diebold_mariano.csv
#   tab_melhor_modelo.csv
#   fig_fcst_<VAR>_H<h>.png
#   fig_heatmap_H<h>.png
#   fig_rmsfe_barras_H<h>.png
# ============================================================

rm(list = ls())

# ============================================================
# 1. PACOTES
# ============================================================

myPKGs <- c("dplyr", "tidyr", "ggplot2", "purrr",
            "stringr", "sandwich", "lmtest", "scales",
            "forcats", "ggtext")

InstalledPKGs    <- names(installed.packages()[, "Package"])
InstallThesePKGs <- myPKGs[!myPKGs %in% InstalledPKGs]
if (length(InstallThesePKGs) > 0)
  install.packages(InstallThesePKGs, repos = "http://cran.us.r-project.org")

invisible(lapply(myPKGs, library, character.only = TRUE))

# ============================================================
# 2. CAMINHOS
# ============================================================

wd <- "C:/Users/felip/OneDrive/Documentos/tcc/Felipe_Dornelles_tcc/"
setwd(wd)

paths <- list(
  data    = "10_data",
  output  = "30_output",
  results = "40_results"
)

if (!dir.exists(paths$results)) dir.create(paths$results, recursive = TRUE)

run_folder <- list.dirs(paths$data,   full.names = TRUE, recursive = FALSE) |>
  sort() |> tail(1)
output_dir <- list.dirs(paths$output, full.names = TRUE, recursive = FALSE) |>
  sort() |> tail(1)

cat(sprintf("Dados de:   %s\n", run_folder))
cat(sprintf("Outputs de: %s\n", output_dir))

# ============================================================
# 3. CARREGA SÉRIES ORIGINAIS
# ============================================================

load(file.path(run_folder, "df_targets.rda"))
load(file.path(run_folder, "df_model.rda"))
load(file.path(run_folder, "targets_br.rda"))

target_names <- unname(unlist(targets_br))
szv          <- length(target_names)
dates_vec    <- df_model$date
bigt         <- length(dates_vec)

mat_targets <- df_targets |>
  select(all_of(target_names)) |>
  as.matrix()

cat(sprintf("T = %d  |  %s a %s\n", bigt,
            format(min(dates_vec), "%b/%Y"),
            format(max(dates_vec), "%b/%Y")))
cat(sprintf("Variáveis alvo: %s\n", paste(target_names, collapse = " | ")))

# ============================================================
# 4. DESCOBRE ARQUIVOS
# ============================================================

rdata_files <- list.files(
  output_dir,
  pattern    = "^TVPfcst_V\\d+_H\\d+_M\\d+\\.RData$",
  full.names = TRUE
)

if (length(rdata_files) == 0)
  stop(sprintf("Nenhum arquivo TVPfcst encontrado em: %s", output_dir))

cat(sprintf("\n%d arquivo(s) encontrado(s).\n", length(rdata_files)))

file_meta <- tibble(path = rdata_files) |>
  mutate(
    fname  = basename(path),
    V      = as.integer(str_extract(fname, "(?<=_V)\\d+")),
    H      = as.integer(str_extract(fname, "(?<=_H)\\d+")),
    M_spec = as.integer(str_extract(fname, "(?<=_M)\\d+"))
  )

horizons <- sort(unique(file_meta$H))
vars_idx <- sort(unique(file_meta$V))
specs    <- sort(unique(file_meta$M_spec))

cat(sprintf("Horizontes H : %s\n", paste(horizons, collapse = ", ")))
cat(sprintf("Variáveis  V : %s\n", paste(vars_idx,  collapse = ", ")))
cat(sprintf("Specs      M : %s\n", paste(specs,     collapse = ", ")))

# ============================================================
# 5. RÓTULOS
# ============================================================

est_labels <- c(
  "1" = "Ridge",
  "2" = "2SRR",
  "3" = "MSRRs",
  "4" = "MSRRd"
)

spec_labels <- c(
  "1" = "AR",
  "2" = "AR+F2",
  "3" = "AR+Tgt",
  "4" = "AR+Pan"
)

# ============================================================
# 6. CARREGA TODOS OS FORECASTS → data.frame longo
# ============================================================

cat("\nCarregando forecasts...\n")

load_one <- function(meta_row) {
  e <- new.env(parent = emptyenv())
  load(meta_row$path, envir = e)

  fc_arr <- e$forecast   # [bigt, max_H, szv, 4]
  v      <- meta_row$V
  h      <- meta_row$H
  m      <- meta_row$M_spec
  n_est  <- dim(fc_arr)[4]

  purrr::map_dfr(seq_len(n_est), function(est) {
    tibble(
      t         = seq_len(bigt),
      date      = dates_vec,
      V         = v,
      H         = h,
      M_spec    = m,
      estimador = est,
      forecast  = fc_arr[, h, v, est],
      realizado = mat_targets[, v]
    )
  })
}

df_long <- purrr::map_dfr(
  purrr::transpose(file_meta),
  load_one
) |>
  mutate(
    var_name   = target_names[V],
    est_label  = est_labels[as.character(estimador)],
    spec_label = spec_labels[as.character(M_spec)]
  )

cat(sprintf("df_long: %d linhas  |  %d combinações únicas\n",
            nrow(df_long),
            n_distinct(df_long |> select(V, H, M_spec, estimador))))

# ============================================================
# 7. ERROS DE PREVISÃO
# CORREÇÃO: usa $ direto para evitar bug de referência no dplyr
# ============================================================

df_erros        <- df_long
df_erros$erro   <- df_erros$realizado - df_erros$forecast
df_erros$se     <- df_erros$erro^2
df_erros$ae     <- abs(df_erros$erro)
df_erros        <- df_erros[!is.na(df_erros$erro), ]

cat(sprintf("df_erros: %d linhas com previsão válida\n", nrow(df_erros)))

# ============================================================
# 8. MSFE, MAE E RMSFE RELATIVO
# ============================================================

df_msfe <- df_erros |>
  group_by(var_name, V, H, M_spec, estimador, est_label, spec_label) |>
  summarise(
    n_obs = n(),
    MSFE  = mean(se),
    MAE   = mean(ae),
    .groups = "drop"
  )

# Benchmark: Ridge (est=1) + especificação AR pura (M_spec=1)
benchmark <- df_msfe |>
  filter(estimador == 1, M_spec == 1) |>
  transmute(var_name, V, H, MSFE_bench = MSFE, n_bench = n_obs)

df_rmsfe <- df_msfe |>
  left_join(benchmark, by = c("var_name", "V", "H")) |>
  mutate(RMSFE_rel = sqrt(MSFE / MSFE_bench))

cat(sprintf("df_rmsfe: %d linhas\n", nrow(df_rmsfe)))

# ============================================================
# 9. TABELA PRINCIPAL — RMSFE relativo
# ============================================================

tab_main <- df_rmsfe |>
  filter(!(estimador == 1 & M_spec == 1)) |>
  mutate(col_label = paste0(spec_label, "_", est_label)) |>
  select(var_name, H, col_label, RMSFE_rel) |>
  pivot_wider(names_from = col_label, values_from = RMSFE_rel) |>
  arrange(var_name, H)

cat("\n=== RMSFE relativo ao Ridge-AR ===\n")
print(tab_main, n = 100, width = 200)

write.csv(tab_main,
          file.path(paths$results, "tab_rmsfe_relativo.csv"),
          row.names = FALSE)

# ============================================================
# 10. TESTE DIEBOLD-MARIANO (HAC, Newey-West)
# ============================================================

dm_test_safe <- function(e_model, e_bench) {
  d <- e_model^2 - e_bench^2
  d <- d[!is.na(d)]
  n <- length(d)

  if (n < 8)
    return(data.frame(DM_stat = NA_real_,
                      p_value = NA_real_,
                      DM_pval_onesided = NA_real_))

  tryCatch({
    df_d <- data.frame(d = d)
    fit  <- lm(d ~ 1, data = df_d)
    ct   <- coeftest(fit,
                     vcov = NeweyWest(fit, lag = 4, prewhite = FALSE))
    dm_s  <- ct[1, "t value"]
    pv_2s <- ct[1, "Pr(>|t|)"]
    pv_1s <- pt(dm_s, df = n - 1, lower.tail = TRUE)
    data.frame(DM_stat = dm_s,
               p_value = pv_2s,
               DM_pval_onesided = pv_1s)
  }, error = function(e)
    data.frame(DM_stat = NA_real_,
               p_value = NA_real_,
               DM_pval_onesided = NA_real_)
  )
}

# Erros do benchmark por (var, H, t)
bench_erros <- df_erros[df_erros$estimador == 1 & df_erros$M_spec == 1,
                        c("var_name", "V", "H", "t", "erro")]
names(bench_erros)[names(bench_erros) == "erro"] <- "erro_bench"

df_erros_dm <- merge(
  df_erros[!(df_erros$estimador == 1 & df_erros$M_spec == 1), ],
  bench_erros,
  by = c("var_name", "V", "H", "t")
)
df_erros_dm <- df_erros_dm[!is.na(df_erros_dm$erro_bench), ]

# Calcula DM por grupo
grupos <- unique(df_erros_dm[, c("var_name", "H", "M_spec",
                                  "estimador", "est_label", "spec_label")])

dm_list <- vector("list", nrow(grupos))

for (i in seq_len(nrow(grupos))) {
  g <- grupos[i, ]

  sub <- df_erros_dm[
    df_erros_dm$var_name  == g$var_name  &
    df_erros_dm$H         == g$H         &
    df_erros_dm$M_spec    == g$M_spec    &
    df_erros_dm$estimador == g$estimador, ]

  dm_res <- dm_test_safe(sub$erro, sub$erro_bench)

  dm_list[[i]] <- cbind(g, dm_res)
}

df_dm <- do.call(rbind, dm_list)
rownames(df_dm) <- NULL

df_dm$sig_2s <- ifelse(df_dm$p_value < 0.01, "***",
                ifelse(df_dm$p_value < 0.05, "**",
                ifelse(df_dm$p_value < 0.10, "*", "")))

df_dm$sig_1s <- ifelse(df_dm$DM_pval_onesided < 0.01, "†††",
                ifelse(df_dm$DM_pval_onesided < 0.05, "††",
                ifelse(df_dm$DM_pval_onesided < 0.10, "†", "")))

cat("\n=== Diebold-Mariano (bilateral) vs. Ridge-AR ===\n")
print(df_dm[order(df_dm$var_name, df_dm$H),
            c("var_name", "H", "spec_label", "est_label",
              "DM_stat", "p_value", "sig_2s")],
      row.names = FALSE)

write.csv(df_dm,
          file.path(paths$results, "tab_diebold_mariano.csv"),
          row.names = FALSE)

# ============================================================
# 11. TABELA PUBLICÁVEL — RMSFE + asteriscos DM
# ============================================================

tab_pub_base <- merge(
  df_rmsfe[!(df_rmsfe$estimador == 1 & df_rmsfe$M_spec == 1), ],
  df_dm[, c("var_name", "H", "M_spec", "estimador", "sig_2s")],
  by = c("var_name", "H", "M_spec", "estimador"),
  all.x = TRUE
)

tab_pub_base$sig_2s[is.na(tab_pub_base$sig_2s)] <- ""
tab_pub_base$col_label <- paste0(tab_pub_base$spec_label, "_",
                                  tab_pub_base$est_label)
tab_pub_base$cell      <- sprintf("%.3f%s",
                                   tab_pub_base$RMSFE_rel,
                                   tab_pub_base$sig_2s)

tab_pub <- tab_pub_base |>
  select(var_name, H, col_label, cell) |>
  pivot_wider(names_from = col_label, values_from = cell) |>
  arrange(var_name, H)

cat("\n=== Tabela publicável (RMSFE + DM) ===\n")
print(tab_pub, n = 100, width = 200)

write.csv(tab_pub,
          file.path(paths$results, "tab_rmsfe_dm_formatado.csv"),
          row.names = FALSE)

# ============================================================
# 12. MELHOR MODELO POR VARIÁVEL × HORIZONTE
# ============================================================

df_rmsfe_nobench <- df_rmsfe[!(df_rmsfe$estimador == 1 &
                                 df_rmsfe$M_spec == 1), ]

best_models <- df_rmsfe_nobench |>
  group_by(var_name, H) |>
  slice_min(RMSFE_rel, n = 1, with_ties = FALSE) |>
  ungroup() |>
  select(var_name, H, spec_label, est_label, RMSFE_rel, n_obs)

cat("\n=== Melhor modelo por variável × horizonte ===\n")
print(best_models)

write.csv(best_models,
          file.path(paths$results, "tab_melhor_modelo.csv"),
          row.names = FALSE)

# ============================================================
# 13. FIGURAS — previsão vs. realizado (melhor modelo)
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

    spec_txt <- unique(df_v$spec_label)
    est_txt  <- unique(df_v$est_label)

    p <- ggplot(df_v, aes(x = date)) +
      geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.4) +
      geom_line(aes(y = realizado, colour = "Realizado"),
                linewidth = 0.75, alpha = 0.85) +
      geom_line(aes(y = forecast, colour = "Previsão"),
                linewidth = 0.75, linetype = "dashed") +
      scale_colour_manual(
        values = c("Realizado" = "#01696f", "Previsão" = "#964219"),
        name   = NULL
      ) +
      scale_x_date(date_labels = "%Y", date_breaks = "3 years") +
      labs(
        title    = sprintf("%s  —  H = %d mês(es)  |  %s / %s",
                           vname, h_val, spec_txt, est_txt),
        subtitle = sprintf("RMSFE relativo = %s", rmsfe_txt),
        x = NULL, y = "Variação"
      ) +
      theme_minimal(base_size = 11) +
      theme(
        legend.position  = "top",
        panel.grid.minor = element_blank(),
        plot.title  = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 9, colour = "grey40")
      )

    fname <- file.path(paths$results,
                       sprintf("fig_fcst_%s_H%d.png", vname, h_val))
    ggsave(fname, p, width = 9, height = 4, dpi = 150)
    cat(sprintf("  Salvo: %s\n", basename(fname)))
  }
}

# ============================================================
# 14. FIGURA — heatmap RMSFE
# ============================================================

for (h_val in horizons) {
  df_heat <- df_rmsfe[df_rmsfe$H == h_val, ]
  df_heat$modelo_full <- paste0(df_heat$spec_label, "\n", df_heat$est_label)
  df_heat$var_name    <- factor(df_heat$var_name, levels = rev(target_names))

  col_order <- unique(df_heat[order(df_heat$M_spec, df_heat$estimador),
                               "modelo_full"])
  df_heat$modelo_full <- factor(df_heat$modelo_full, levels = col_order)

  p_heat <- ggplot(df_heat,
                   aes(x = modelo_full, y = var_name, fill = RMSFE_rel)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(
      aes(label  = sprintf("%.3f", RMSFE_rel),
          colour = RMSFE_rel > 1.05),
      size = 2.8, fontface = "bold", show.legend = FALSE
    ) +
    scale_fill_gradient2(
      low      = "#01696f",
      mid      = "#f7f6f2",
      high     = "#a12c7b",
      midpoint = 1,
      limits   = c(0.6, 1.4),
      oob      = scales::squish,
      name     = "RMSFE\nrelativo"
    ) +
    scale_colour_manual(values = c("FALSE" = "grey20", "TRUE" = "grey60")) +
    labs(
      title    = sprintf("RMSFE relativo ao Ridge-AR  —  H = %d", h_val),
      subtitle = "< 1 = melhor que benchmark  |  * p<0.10  ** p<0.05  *** p<0.01 (DM bilateral)",
      x = NULL, y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid    = element_blank(),
      plot.title    = element_text(face = "bold"),
      plot.subtitle = element_text(size = 8, colour = "grey40"),
      axis.text.x   = element_text(size = 8, lineheight = 1.2)
    )

  fname <- file.path(paths$results, sprintf("fig_heatmap_H%d.png", h_val))
  ggsave(fname, p_heat, width = 12, height = 5, dpi = 150)
  cat(sprintf("  Salvo: %s\n", basename(fname)))
}

# ============================================================
# 15. FIGURA — barras RMSFE por variável
# ============================================================

for (h_val in horizons) {
  df_bar <- df_rmsfe[df_rmsfe$H == h_val &
                       !(df_rmsfe$estimador == 1 & df_rmsfe$M_spec == 1), ]
  df_bar$modelo_full <- paste0(df_bar$spec_label, " / ", df_bar$est_label)
  df_bar$var_name    <- factor(df_bar$var_name, levels = target_names)
  df_bar$melhor      <- df_bar$RMSFE_rel < 1

  p_bar <- ggplot(df_bar,
                  aes(x = reorder(modelo_full, RMSFE_rel),
                      y = RMSFE_rel,
                      fill = melhor)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept = 1, linetype = "dashed",
               colour = "grey40", linewidth = 0.5) +
    scale_fill_manual(
      values = c("TRUE" = "#01696f", "FALSE" = "#964219"),
      guide  = "none"
    ) +
    scale_y_continuous(labels = number_format(accuracy = 0.01)) +
    facet_wrap(~var_name, scales = "free_y", ncol = 2) +
    coord_flip() +
    labs(
      title    = sprintf("RMSFE relativo ao Ridge-AR  —  H = %d", h_val),
      subtitle = "Verde = melhor que benchmark  |  Laranja = pior",
      x = NULL, y = "RMSFE relativo"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major.y = element_blank(),
      strip.text         = element_text(face = "bold"),
      plot.title         = element_text(face = "bold")
    )

  fname <- file.path(paths$results,
                     sprintf("fig_rmsfe_barras_H%d.png", h_val))
  ggsave(fname, p_bar,
         width  = 11,
         height = 3 * ceiling(szv / 2),
         dpi    = 150)
  cat(sprintf("  Salvo: %s\n", basename(fname)))
}

# ============================================================
# 16. SUMÁRIO FINAL
# ============================================================

cat("\n============================================================\n")
cat(sprintf("Arquivos salvos em: %s\n", paths$results))
cat("------------------------------------------------------------\n")
cat("CSVs:\n")
cat("  tab_rmsfe_relativo.csv\n")
cat("  tab_rmsfe_dm_formatado.csv   <- tabela TCC (publicável)\n")
cat("  tab_diebold_mariano.csv\n")
cat("  tab_melhor_modelo.csv\n")
cat("Figuras:\n")
cat("  fig_fcst_<VAR>_H<h>.png      <- previsão vs. realizado\n")
cat("  fig_heatmap_H<h>.png         <- heatmap RMSFE\n")
cat("  fig_rmsfe_barras_H<h>.png    <- barras por variável\n")
cat("============================================================\n")