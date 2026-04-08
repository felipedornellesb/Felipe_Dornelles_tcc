# ============================================================
# 03_02_results_nathalia_comparison.R
# Avaliação de previsão: RMSPE, DM, CSFE, MCS
# ============================================================

rm(list = ls())

library(tidyverse)
library(forecast)
library(MCS)
library(xtable)
library(ggplot2)
library(ggsci)
library(ggpubr)

# ============================================================
# 0. CAMINHOS E PARÂMETROS
# ============================================================

paths <- list(
  forecasts = "forecasts/",
  results   = "40_results/40_nathalia_comparative/"   # <-- ALTERADO
)
dir.create(paths$results, showWarnings = FALSE, recursive = TRUE)  # recursive = TRUE garante que cria subpastas se não existirem

UNRATE <- "SEADE12_TDTGSP12"
IPCA   <- "PRECOS12_IPCA12"
SPREAD <- "JPM366_EMBI366"
variable     <- c(UNRATE, IPCA, SPREAD)
horizon_list <- c(1, 3, 6, 9, 12)

var_labels <- c(
  "SEADE12_TDTGSP12" = "Desemprego (UNRATE)",
  "PRECOS12_IPCA12"  = "Inflação (IPCA)",
  "JPM366_EMBI366"   = "Risco-País (SPREAD)"
)

# ============================================================
# 1. CARREGAR PREVISÕES
# ============================================================

actual_values <- readRDS(file.path(paths$forecasts, "actual_values.rda"))
test_matrix   <- as.matrix(actual_values)

ar_bic_rmspe <- readRDS(file.path(paths$forecasts, "ar_bic_results.rda"))
model_ar_bic <- readRDS(file.path(paths$forecasts, "ar_bic.rda"))

tvp_2srr <- readRDS(file.path(paths$forecasts, "tvp_2srr.rda"))
ar        <- readRDS(file.path(paths$forecasts, "ar.rda"))
ols       <- readRDS(file.path(paths$forecasts, "ols.rda"))
ridge_1s  <- readRDS(file.path(paths$forecasts, "ridge_1srr.rda"))

model_list  <- c(tvp_2srr, ar, ols, ridge_1s)
models_list <- c("ar_aic", "ar_cv", "ols_bic", "ridge_1srr_cv", "tvp_2srr_cv")

# ============================================================
# 2. CALCULAR RMSPE RELATIVO, DM E CSFE
# ============================================================

get_stars <- function(p) {
  ifelse(is.na(p), "",
  ifelse(p < 0.01, "***",
  ifelse(p < 0.05, "**",
  ifelse(p < 0.10, "*", ""))))
}

rmspe_data         <- data.frame()
matrix_list        <- list()
CSFE               <- list()
CSFE_list          <- list()
squared_error_list <- list()

for (v in variable) {
  for (h in horizon_list) {

    bic_filtered     <- model_ar_bic[sapply(model_ar_bic,
                          function(x) x$horizon == h & x$variable == v &
                                      x$model == "ar_bic_bic")]
    model_matrix_bic <- sapply(bic_filtered, function(x) x$pred)
    pred_error_bic   <- (test_matrix - model_matrix_bic)
    bic_finite       <- pred_error_bic[is.finite(pred_error_bic)]

    for (model in models_list) {

      mf           <- model_list[sapply(model_list,
                        function(x) x$horizon == h & x$variable == v &
                                    x$model == model)]
      model_matrix <- sapply(mf, function(x) x$pred)
      pred_error   <- (test_matrix - model_matrix)
      mod_finite   <- pred_error[is.finite(pred_error)]

      matrix_list[[model]] <- mod_finite^2
      CSFE[[model]]        <- cumsum(bic_finite^2 - mod_finite^2)

      rmspe_rel <- sqrt(mean(mod_finite^2, na.rm = TRUE)) /
                   ar_bic_rmspe$RMSPE[ar_bic_rmspe$variable == v &
                                      ar_bic_rmspe$horizon  == h]

      if (all(pred_error == pred_error_bic, na.rm = TRUE)) {
        dm_pval <- 1
      } else {
        dm_res  <- forecast::dm.test(pred_error, pred_error_bic,
                                     h = h, varestimator = "bartlett")
        dm_pval <- dm_res$p.value
      }

      rmspe_data <- rbind(rmspe_data,
        data.frame(variable = v, horizon = h, model = model,
                   RMSPE   = rmspe_rel,
                   DM_pval = format(dm_pval, scientific = FALSE),
                   stringsAsFactors = FALSE))
    }

    squared_error_list[[paste(v, h, sep = "_")]] <- do.call(cbind, matrix_list)
    CSFE_list[[paste(v, h, sep = "_")]]          <- do.call(cbind, CSFE)
  }
}

results <- list(
  rmspe_data    = rmspe_data,
  squared_error = squared_error_list,
  CSFE          = CSFE_list
)

# ============================================================
# 3. TABELAS RMSPE COM ESTRELAS
# ============================================================

rmspe_full <- bind_rows(
  ar_bic_rmspe %>% mutate(DM_pval = "1", model = "ar_bic"),
  rmspe_data
) %>%
  arrange(variable, horizon) %>%
  mutate(
    stars = get_stars(as.numeric(DM_pval)),
    RMSPE = round(RMSPE, 4)
  ) %>%
  unite("values", RMSPE, stars, sep = "", remove = FALSE)

make_rmspe_table <- function(v) {
  rmspe_full %>%
    filter(variable == v) %>%
    select(model, horizon, values) %>%
    pivot_wider(names_from   = horizon,
                names_prefix = "h=",
                values_from  = values,
                values_fill  = NULL)
}

tab_SPREAD <- make_rmspe_table(SPREAD)
tab_IPCA   <- make_rmspe_table(IPCA)
tab_UNRATE <- make_rmspe_table(UNRATE)

print(xtable(tab_SPREAD, caption = "RMSPE Relativo — SPREAD"), type = "latex")
print(xtable(tab_IPCA,   caption = "RMSPE Relativo — IPCA"),   type = "latex")
print(xtable(tab_UNRATE, caption = "RMSPE Relativo — UNRATE"), type = "latex")

write.csv(tab_SPREAD, file.path(paths$results, "TAB_rmspe_SPREAD.csv"), row.names = FALSE)
write.csv(tab_IPCA,   file.path(paths$results, "TAB_rmspe_IPCA.csv"),   row.names = FALSE)
write.csv(tab_UNRATE, file.path(paths$results, "TAB_rmspe_UNRATE.csv"), row.names = FALSE)

# ============================================================
# 4. MCS
# ============================================================

run_mcs <- function(v, h) {
  mat <- results$squared_error[[paste(v, h, sep = "_")]]
  mat <- mat[, colSums(is.na(mat)) == 0, drop = FALSE]
  if (ncol(mat) < 2) return(NULL)
  MCSprocedure(mat)
}

mcs_results <- list()
for (v in variable) {
  for (h in horizon_list) {
    key <- paste(v, h, sep = "_")
    mcs_results[[key]] <- run_mcs(v, h)
    cat(sprintf("\nMCS — %s | h=%d\n", var_labels[v], h))
    print(mcs_results[[key]])
  }
}

# ============================================================
# 5. GRÁFICOS CSFE
# ============================================================

date_start  <- as.Date("2009-01-01")
n_forecasts <- nrow(test_matrix)

make_csfe_df <- function(v, h, top_models) {
  mat <- results$CSFE[[paste(v, h, sep = "_")]]
  as.data.frame(mat) %>%
    mutate(dates = seq(date_start, by = "month", length.out = n_forecasts)) %>%
    select(all_of(c(top_models, "dates"))) %>%
    pivot_longer(-dates, names_to = "Model", values_to = "value")
}

make_csfe_plot <- function(df, h_label) {
  ggplot(df, aes(x = dates, y = value, color = Model)) +
    geom_line(linetype = "solid", linewidth = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    scale_color_npg() +
    labs(x = "Tempo", y = NULL, title = sprintf("h = %s", h_label)) +
    theme_bw() +
    theme(
      axis.text.x  = element_text(size = 10),
      legend.title = element_text(size = 10),
      legend.text  = element_text(size = 9),
      plot.title   = element_text(size = 14)
    )
}

top_models_all <- c("ar_aic", "ar_cv", "ols_bic", "ridge_1srr_cv", "tvp_2srr_cv")

make_var_panel <- function(v, title_label) {
  plots <- lapply(horizon_list, function(h) {
    df <- make_csfe_df(v, h, top_models_all)
    make_csfe_plot(df, h_label = as.character(h))
  })

  panel <- ggarrange(plotlist = plots, ncol = 2, nrow = 3,
                     common.legend = TRUE, legend = "bottom")

  panel_titled <- annotate_figure(panel,
    top = text_grob(
      sprintf("CSFE Cumulativo vs AR-BIC — %s", title_label),
      face = "bold", size = 14
    )
  )

  ggsave(
    filename = file.path(paths$results, sprintf("FIG_CSFE_%s.png", v)),
    plot     = panel_titled,
    width    = 12, height = 14, dpi = 300
  )

  panel_titled
}

plot_SPREAD <- make_var_panel(SPREAD, var_labels[SPREAD])
plot_IPCA   <- make_var_panel(IPCA,   var_labels[IPCA])
plot_UNRATE <- make_var_panel(UNRATE, var_labels[UNRATE])

# ============================================================
# 6. GRÁFICO RMSPE POR HORIZONTE
# ============================================================

p_rmspe <- rmspe_full %>%
  filter(model %in% top_models_all) %>%
  mutate(horizon   = as.integer(horizon),
         var_label = var_labels[variable]) %>%
  ggplot(aes(x = horizon, y = RMSPE, color = model, group = model)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  facet_wrap(~ var_label, scales = "free_y", ncol = 1) +
  scale_color_npg() +
  scale_x_continuous(breaks = horizon_list) +
  labs(
    title    = "RMSPE Relativo ao AR-BIC por Horizonte",
    subtitle = "< 1 = melhor que o benchmark | linha tracejada = benchmark",
    x        = "Horizonte (meses)", y = "RMSPE Relativo",
    color    = "Modelo"
  ) +
  theme_bw() +
  theme(strip.text = element_text(face = "bold", size = 11),
        legend.position = "bottom")

ggsave(
  filename = file.path(paths$results, "FIG_RMSPE_por_horizonte.png"),
  plot     = p_rmspe,
  width    = 10, height = 12, dpi = 300
)

cat(sprintf("\n✅ 03_results.R concluído. Resultados em: %s\n", paths$results))