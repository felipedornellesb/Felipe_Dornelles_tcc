# ============================================================
# 02_FORECAST.R
#
# Forecast TVP-2SRR adaptado ao 01_data_prep.R do projeto.
# Versão enxuta para a base brasileira:
#   - usa df_model.rda (base final limpa)
#   - usa df_targets.rda e df_panel_pca.rda quando disponíveis
#   - sem dataprep exagerado
#   - lags curtos (1:3)
#   - poucos fatores (K = 1:3)
#   - janela inicial definida pela amostra disponível
# ============================================================

rm(list = ls())

# ============================================================
# SETUP
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

run_date   <- format(Sys.Date(), "%m_%d_%Y")
output_run <- file.path(paths$output, paste0("outputs_", run_date))

for (p in c(paths$output, output_run, paths$results)) {
  if (!dir.exists(p)) dir.create(p, recursive = TRUE)
}

cat(sprintf("Output folder: %s\n", output_run))

# ============================================================
# PACKAGES
# ============================================================

myPKGs <- c("dplyr", "pracma", "forecast", "MTS", "matrixcalc", "fGarch")

InstalledPKGs    <- names(installed.packages()[, "Package"])
InstallThesePKGs <- myPKGs[!myPKGs %in% InstalledPKGs]
if (length(InstallThesePKGs) > 0)
  install.packages(InstallThesePKGs, repos = "http://cran.us.r-project.org")

invisible(lapply(myPKGs, library, character.only = TRUE))

# ============================================================
# LOAD FUNCTIONS
# ============================================================

source(file.path(paths$functions, "MV2SRR_v221103.R"))

# ============================================================
# LOAD LATEST DATA RUN
# ============================================================

data_subfolders <- list.dirs(paths$data, recursive = FALSE, full.names = TRUE)
data_subfolders <- data_subfolders[grepl("data_\\d{2}_\\d{2}_\\d{4}$",
                                         basename(data_subfolders))]

if (length(data_subfolders) == 0) {
  stop("Nenhuma pasta data_MM_DD_YYYY encontrada em 10_data/. Rode o 01_data_prep.R primeiro.")
}

run_folder <- data_subfolders[length(data_subfolders)]
cat(sprintf("Loading data from: %s\n", run_folder))

load(file.path(run_folder, "df_model.rda"))

if (file.exists(file.path(run_folder, "df_targets.rda"))) {
  load(file.path(run_folder, "df_targets.rda"))
} else {
  df_targets <- NULL
}

if (file.exists(file.path(run_folder, "df_panel_pca.rda"))) {
  load(file.path(run_folder, "df_panel_pca.rda"))
} else {
  df_panel_pca <- NULL
}

if (file.exists(file.path(run_folder, "targets_br.rda"))) {
  load(file.path(run_folder, "targets_br.rda"))
} else {
  targets_br <- list(
    V1 = "PIB",
    V2 = "IPCA",
    V3 = "SELIC",
    V4 = "CAMBIO",
    V5 = "DESEMPREGO"
  )
}

# ============================================================
# PARAMETERS
# ============================================================

variable_list <- unname(unlist(targets_br))
variable_list <- variable_list[variable_list %in% names(df_model)]

horizon_list <- c(1, 2, 4)
lag_orders   <- 1:3
K_grid       <- 1:3
model_name   <- "coulombe_2srr_br"
window_type  <- "expanding"
reopt_every  <- 12L

# ============================================================
# PREPARE MODEL BASE
# ============================================================

if (!"date" %in% names(df_model)) {
  stop("df_model precisa conter a coluna date.")
}

model_df <- df_model |>
  dplyr::select(-date)

# Mantém apenas colunas finitas
keep_cols <- names(model_df)[sapply(model_df, function(x) all(is.finite(x)))]
model_df  <- model_df[, keep_cols, drop = FALSE]

# Confirma alvos presentes
ausentes <- variable_list[!variable_list %in% colnames(model_df)]
if (length(ausentes) > 0) {
  stop(paste("Variáveis-alvo ausentes em df_model:", paste(ausentes, collapse = ", ")))
}

# Padronização do painel
scaled_mat <- scale(model_df)
df_scaled  <- as.data.frame(scaled_mat)

n_total <- nrow(df_scaled)

# Janela inicial: 70% da amostra, com mínimo de 60 observações
start_window <- max(60L, floor(n_total * 0.70))
if (start_window >= n_total) start_window <- n_total - 1L

end_window <- n_total
nwindows   <- end_window - start_window + 1L

cat(sprintf("Panel: %d obs x %d vars | OOS windows: %d (index %d to %d)\n",
            n_total, ncol(df_scaled), nwindows, start_window, end_window))

# ============================================================
# HELPER FUNCTIONS
# ============================================================

remove_high_corr <- function(x_mat, cutoff = 0.98) {
  if (ncol(x_mat) <= 1L) return(x_mat)
  cor_mat <- suppressWarnings(cor(x_mat, use = "pairwise.complete.obs"))
  keep <- rep(TRUE, ncol(x_mat))
  for (j in 2:ncol(x_mat)) {
    if (any(abs(cor_mat[j, 1:(j - 1)]) > cutoff, na.rm = TRUE)) {
      keep[j] <- FALSE
    }
  }
  x_mat[, keep, drop = FALSE]
}

prep_window <- function(df, variable, horizon, lag_y = 2L, K = 2L) {
  y <- as.numeric(df[[variable]])
  x <- as.matrix(df[, setdiff(names(df), variable), drop = FALSE])

  if (any(!is.finite(y))) stop(sprintf("y inválido para %s", variable))

  # remove colunas ruins
  if (ncol(x) > 0L) {
    ok_fin <- apply(x, 2, function(col) all(is.finite(col)))
    x <- x[, ok_fin, drop = FALSE]
  }

  if (ncol(x) > 0L) {
    ok_var <- apply(x, 2, function(col) stats::sd(col) > 1e-8)
    x <- x[, ok_var, drop = FALSE]
  }

  if (ncol(x) > 1L) {
    x <- remove_high_corr(x, cutoff = 0.98)
  }

  if (ncol(x) < 1L) stop(sprintf("Sem regressoras válidas para %s", variable))

  K_actual <- min(K, ncol(x), nrow(x) - 1L)
  if (K_actual < 1L) stop(sprintf("K inválido para %s", variable))

  pca <- prcomp(x, scale. = FALSE, center = FALSE)
  factors <- pca$x[, 1:K_actual, drop = FALSE]

  max_lag <- lag_y
  n <- nrow(df)

  start_idx <- max_lag + 1L
  end_idx   <- n - horizon

  if (end_idx < start_idx) {
    stop(sprintf("Janela pequena demais: n=%d, h=%d, lag=%d", n, horizon, lag_y))
  }

  Xin <- NULL
  yin <- NULL

  for (t in start_idx:end_idx) {
    y_tgt <- y[t + horizon]
    y_lags <- y[(t - 1L):(t - lag_y)]
    f_t <- factors[t, , drop = TRUE]
    row_x <- c(y_lags, f_t)
    Xin <- rbind(Xin, row_x)
    yin <- c(yin, y_tgt)
  }

  t_last <- n
  y_lags_last <- y[(t_last):(t_last - lag_y + 1L)]
  f_last <- factors[n, , drop = TRUE]
  Xout <- c(y_lags_last, f_last)

  Xin <- as.matrix(Xin)
  yin <- as.numeric(yin)
  Xout <- as.numeric(Xout)

  good_rows <- apply(Xin, 1, function(z) all(is.finite(z))) & is.finite(yin)
  Xin <- Xin[good_rows, , drop = FALSE]
  yin <- yin[good_rows]

  if (nrow(Xin) < 25L) {
    stop(sprintf("Poucas observações úteis: %d", nrow(Xin)))
  }

  list(Xin = Xin, yin = yin, Xout = Xout)
}

run_tvp_2srr <- function(df, variable, horizon, lag_y, K, best = list(), reopt = FALSE) {
  prep <- prep_window(df = df, variable = variable, horizon = horizon, lag_y = lag_y, K = K)

  lambdavec <- exp(pracma::linspace(-4, 10, n = 12))
  if (is.null(best$lambda)) reopt <- TRUE

  lambda_use <- if (reopt) lambdavec else best$lambda

  out <- tvp.ridge(
    X                 = prep$Xin,
    Y                 = matrix(prep$yin, ncol = 1),
    lambda.candidates = lambda_use,
    oosX              = prep$Xout,
    kfold             = 5L,
    CV.2SRR           = TRUE,
    CV.plot           = FALSE,
    sig.eps.param     = 0.75,
    sig.u.param       = 0.75
  )

  if (reopt) {
    best <- list(lambda = out$lambdas, lag = lag_y, K = K)
  }

  list(pred = as.numeric(out$forecast), best = best)
}

# ============================================================
# FORECAST LOOP
# ============================================================

set.seed(1234)
forecasts_list <- list()

for (v in variable_list) {
  for (h in horizon_list) {

    cat(sprintf("\n[%s] var=%-12s | h=%d\n", model_name, v, h))

    best <- list()
    model_list <- vector("list", nwindows)

    for (i in seq_len(nwindows)) {
      train_end <- start_window - h - 1L + i

      if (train_end < 40L) {
        model_list[[i]] <- list(pred = NA_real_)
        next
      }

      Df <- if (window_type == "expanding") {
        df_scaled[1:train_end, , drop = FALSE]
      } else {
        df_scaled[max(1L, train_end - 119L):train_end, , drop = FALSE]
      }

      reopt <- (i == 1L) || (i %% reopt_every == 1L)
      if (reopt) cat(sprintf("  re-optimizing lambda (window %d/%d)\n", i, nwindows))

      model <- tryCatch({
        run_tvp_2srr(
          df       = Df,
          variable = v,
          horizon  = h,
          lag_y    = 2L,
          K        = 2L,
          best     = best,
          reopt    = reopt
        )
      }, error = function(e) {
        message(sprintf("  [ERROR tvp.ridge] var=%s h=%d win=%d: %s", v, h, i, e$message))
        list(pred = NA_real_, best = best)
      })

      if (is.na(model$pred)) {
        model$pred <- if (i > 1L && !is.na(model_list[[i - 1L]]$pred)) {
          model_list[[i - 1L]]$pred
        } else {
          0
        }
        cat(sprintf("  [fallback] window %d uses previous pred\n", i))
      }

      model_list[[i]] <- list(pred = model$pred)
      best <- model$best
    }

    preds <- vapply(model_list, function(x) x$pred, numeric(1L))

    actual_idx_end <- min(end_window, start_window + nwindows - 1L)
    actual_vec <- as.numeric(df_scaled[start_window:actual_idx_end, v])

    forecasts_list[[length(forecasts_list) + 1L]] <- list(
      variable = v,
      horizon  = h,
      model    = model_name,
      pred     = preds,
      actual   = actual_vec
    )

    cat(sprintf("  -> %d forecasts | NAs: %d\n", length(preds), sum(is.na(preds))))
  }
}

# ============================================================
# SAVE FORECASTS
# ============================================================

saveRDS(forecasts_list, file = file.path(output_run, paste0(model_name, ".rds")))

cat(sprintf("\nForecast completo. %d combinações salvas em %s/\n",
            length(forecasts_list), output_run))

# ============================================================
# CHECKUPS
# ============================================================

results <- do.call(rbind, lapply(forecasts_list, function(f) {
  n <- min(length(f$pred), length(f$actual))
  err <- f$pred[1:n] - f$actual[1:n]
  data.frame(
    variable = f$variable,
    horizon  = f$horizon,
    MSFE     = mean(err^2, na.rm = TRUE),
    RMSFE    = sqrt(mean(err^2, na.rm = TRUE)),
    MAE      = mean(abs(err), na.rm = TRUE)
  )
}))

rw_results <- do.call(rbind, lapply(forecasts_list, function(f) {
  act <- f$actual
  h   <- f$horizon
  n   <- length(act)
  if (n <= h) return(NULL)
  rw_pred <- act[1:(n - h)]
  rw_act  <- act[(1 + h):n]
  rw_err  <- rw_pred - rw_act
  data.frame(
    variable = f$variable,
    horizon  = h,
    MSFE_RW  = mean(rw_err^2, na.rm = TRUE)
  )
}))

results_rel <- merge(results, rw_results, by = c("variable", "horizon"), all.x = TRUE)
results_rel$MSFE_rel <- results_rel$MSFE / results_rel$MSFE_RW

dm_results <- do.call(rbind, lapply(forecasts_list, function(f) {
  act <- f$actual
  h   <- f$horizon
  n   <- min(length(f$pred), length(act))

  e1 <- f$pred[1:n] - act[1:n]
  e2 <- rep(NA_real_, n)
  if (n > h) {
    for (i in (h + 1L):n) e2[i] <- act[i - h] - act[i]
  }

  keep <- !is.na(e1) & !is.na(e2)
  if (sum(keep) < 10L) {
    return(data.frame(variable = f$variable,
                      horizon  = h,
                      DM_stat  = NA_real_,
                      p_value  = NA_real_,
                      MSFE_rel = NA_real_))
  }

  dm_test <- tryCatch(
    forecast::dm.test(e1[keep], e2[keep], alternative = "less", h = h, power = 2),
    error = function(e) NULL
  )

  data.frame(
    variable = f$variable,
    horizon  = h,
    DM_stat  = if (!is.null(dm_test)) as.numeric(dm_test$statistic) else NA_real_,
    p_value  = if (!is.null(dm_test)) dm_test$p.value else NA_real_,
    MSFE_rel = mean(e1[keep]^2, na.rm = TRUE) / mean(e2[keep]^2, na.rm = TRUE)
  )
}))

results_summary <- list(
  fit_stats = results,
  msfe_rel  = results_rel,
  dm_vs_rw  = dm_results
)

saveRDS(results_summary, file = file.path(paths$results, "results_summary.rds"))
write.csv(results,      file = file.path(paths$results, "TAB01_fit_stats.csv"), row.names = FALSE)
write.csv(results_rel,  file = file.path(paths$results, "TAB02_msfe_relativo.csv"), row.names = FALSE)
write.csv(dm_results,   file = file.path(paths$results, "TAB03_diebold_mariano.csv"), row.names = FALSE)

cat("\n=== RESUMO GERAL ===\n")
cat(sprintf("Modelo      : %s\n", model_name))
cat(sprintf("Janela      : %s\n", window_type))
cat(sprintf("Variáveis   : %s\n", paste(variable_list, collapse = ", ")))
cat(sprintf("Horizontes  : %s\n", paste(horizon_list, collapse = ", ")))
cat(sprintf("Combinações : %d\n", length(forecasts_list)))
cat(sprintf("OOS windows : %d\n", nwindows))
cat(sprintf("Período OOS : índices %d a %d\n", start_window, end_window))
cat("\n02_FORECAST.R finalizado.\n")

