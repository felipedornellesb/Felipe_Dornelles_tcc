# ============================================================
# 02_03_forecast.R
#
# TVP-2SRR forecast pipeline revisado para o projeto Brasil.
# Ajustes desta versão:
#   - sem fallback enganoso: falha vira NA
#   - log limpo de execução salvo em arquivo
#   - contagem de falhas por variável/horizonte
#   - uso explícito de grid lag x K
#   - seleção de hiperparâmetros por OOS validation dentro da janela inicial
#   - painel PCA separado dos targets quando disponível
#   - checagem explícita de dependências do MV2SRR_v221103.R
#   - avaliação final apenas sobre previsões realmente estimadas
#
# Observação metodológica:
#   Esta versão aproxima melhor a lógica de forecasting recursivo do tipo
#   Coulombe: fatores extraídos do painel, horizonte explícito, escolha de
#   hiperparâmetros em tempo real e comparação com random walk.
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

log_file <- file.path(output_run, "02_03_forecast_log.txt")
if (file.exists(log_file)) file.remove(log_file)

log_line <- function(...) {
  txt <- paste0(...)
  cat(txt, "\n")
  cat(txt, "\n", file = log_file, append = TRUE)
}

log_line(sprintf("Output folder: %s", output_run))

# ============================================================
# PACKAGES
# ============================================================

myPKGs <- c(
  "dplyr", "pracma", "forecast", "MTS", "matrixcalc", "fGarch",
  "tibble", "stats"
)

InstalledPKGs    <- names(installed.packages()[, "Package"])
InstallThesePKGs <- myPKGs[!myPKGs %in% InstalledPKGs]
if (length(InstallThesePKGs) > 0) {
  install.packages(InstallThesePKGs, repos = "https://cran.r-project.org")
}

invisible(lapply(myPKGs, library, character.only = TRUE))

# Compatibilidade explícita para scripts legados que chamam vec() sem namespace
if (!exists("vec", mode = "function")) {
  if ("package:MTS" %in% search() && exists("vec", where = as.environment("package:MTS"), mode = "function")) {
    vec <- get("vec", envir = as.environment("package:MTS"))
  }
}

if (!exists("garchFit", mode = "function")) {
  stop("garchFit() não disponível. Verifique instalação/carregamento do pacote fGarch.")
}

# ============================================================
# LOAD FUNCTIONS
# ============================================================

mv2_file <- file.path(paths$functions, "MV2SRR_v221103.R")
if (!file.exists(mv2_file)) {
  stop(sprintf("Arquivo não encontrado: %s", mv2_file))
}

source(mv2_file)

if (!exists("tvp.ridge", mode = "function")) {
  stop("tvp.ridge() não foi carregada a partir de MV2SRR_v221103.R")
}

log_line("MV2SRR_v221103.R carregado com sucesso.")

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
log_line(sprintf("Loading data from: %s", run_folder))

load(file.path(run_folder, "df_model.rda"))
if (file.exists(file.path(run_folder, "df_targets.rda"))) load(file.path(run_folder, "df_targets.rda")) else df_targets <- NULL
if (file.exists(file.path(run_folder, "df_panel_pca.rda"))) load(file.path(run_folder, "df_panel_pca.rda")) else df_panel_pca <- NULL
if (file.exists(file.path(run_folder, "targets_br.rda"))) {
  load(file.path(run_folder, "targets_br.rda"))
} else {
  targets_br <- list(V1 = "PIB", V2 = "IPCA", V3 = "SELIC", V4 = "CAMBIO", V5 = "DESEMPREGO")
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
min_train_n  <- 60L
min_eval_n   <- 25L
corr_cutoff  <- 0.98
lambda_grid  <- exp(pracma::linspace(-4, 10, n = 12))

# ============================================================
# PREPARE MODEL BASE
# ============================================================

if (!"date" %in% names(df_model)) stop("df_model precisa conter a coluna date.")

if (!is.null(df_panel_pca) && "date" %in% names(df_panel_pca)) {
  panel_df <- df_panel_pca |>
    dplyr::select(-date)
} else {
  panel_candidates <- setdiff(names(df_model), c("date", variable_list))
  panel_df <- df_model[, panel_candidates, drop = FALSE]
}

target_df <- df_model[, variable_list, drop = FALSE]

keep_panel <- names(panel_df)[sapply(panel_df, function(x) all(is.finite(x)) && stats::sd(x) > 1e-8)]
panel_df   <- panel_df[, keep_panel, drop = FALSE]

keep_target <- names(target_df)[sapply(target_df, function(x) all(is.finite(x)))]
target_df   <- target_df[, keep_target, drop = FALSE]
variable_list <- intersect(variable_list, names(target_df))

if (length(variable_list) == 0) stop("Nenhuma variável-alvo válida encontrada.")
if (ncol(panel_df) == 0) stop("Painel PCA vazio após limpeza. O modelo Coulombe precisa de regressoras de painel.")

scaled_target <- scale(target_df)
scaled_panel  <- scale(panel_df)

df_scaled <- as.data.frame(cbind(scaled_target, scaled_panel))

n_total <- nrow(df_scaled)
start_window <- max(min_train_n, floor(n_total * 0.70))
if (start_window >= n_total) start_window <- n_total - 1L
end_window <- n_total
nwindows   <- end_window - start_window + 1L

log_line(sprintf("Panel: %d obs x %d vars | OOS windows: %d (index %d to %d)",
                 n_total, ncol(df_scaled), nwindows, start_window, end_window))
log_line(sprintf("Targets: %s", paste(variable_list, collapse = ", ")))
log_line(sprintf("Panel vars for PCA: %d", ncol(panel_df)))

# ============================================================
# HELPER FUNCTIONS
# ============================================================

remove_high_corr <- function(x_mat, cutoff = 0.98) {
  if (is.null(dim(x_mat)) || ncol(x_mat) <= 1L) return(x_mat)
  cor_mat <- suppressWarnings(cor(x_mat, use = "pairwise.complete.obs"))
  keep <- rep(TRUE, ncol(x_mat))
  for (j in 2:ncol(x_mat)) {
    if (any(abs(cor_mat[j, 1:(j - 1)]) > cutoff, na.rm = TRUE)) keep[j] <- FALSE
  }
  x_mat[, keep, drop = FALSE]
}

make_design <- function(df, variable, horizon, lag_y = 2L, K = 2L, corr_cutoff = 0.98, panel_names = NULL) {
  y <- as.numeric(df[[variable]])

  if (is.null(panel_names)) panel_names <- setdiff(names(df), variable)
  panel_names <- intersect(panel_names, names(df))
  x <- as.matrix(df[, panel_names, drop = FALSE])

  if (any(!is.finite(y))) stop(sprintf("y inválido para %s", variable))
  if (ncol(x) < 1L) stop(sprintf("Painel vazio para %s", variable))

  ok_fin <- apply(x, 2, function(col) all(is.finite(col)))
  x <- x[, ok_fin, drop = FALSE]
  if (ncol(x) < 1L) stop(sprintf("Sem regressoras finitas para %s", variable))

  ok_var <- apply(x, 2, function(col) stats::sd(col) > 1e-8)
  x <- x[, ok_var, drop = FALSE]
  if (ncol(x) < 1L) stop(sprintf("Sem regressoras com variância para %s", variable))

  if (ncol(x) > 1L) x <- remove_high_corr(x, cutoff = corr_cutoff)
  if (ncol(x) < 1L) stop(sprintf("Sem regressoras após filtro de correlação para %s", variable))

  K_actual <- min(K, ncol(x), nrow(x) - 1L)
  if (K_actual < 1L) stop(sprintf("K inválido para %s", variable))

  pca <- prcomp(x, scale. = FALSE, center = FALSE)
  factors <- pca$x[, 1:K_actual, drop = FALSE]

  n <- nrow(df)
  start_idx <- lag_y + 1L
  end_idx   <- n - horizon

  if (end_idx < start_idx) {
    stop(sprintf("Janela pequena demais: n=%d, h=%d, lag=%d", n, horizon, lag_y))
  }

  Xin <- vector("list", end_idx - start_idx + 1L)
  yin <- numeric(end_idx - start_idx + 1L)

  pos <- 1L
  for (t in start_idx:end_idx) {
    y_tgt  <- y[t + horizon]
    y_lags <- y[(t - 1L):(t - lag_y)]
    f_t    <- factors[t, , drop = TRUE]
    Xin[[pos]] <- c(y_lags, f_t)
    yin[pos]   <- y_tgt
    pos <- pos + 1L
  }

  Xin <- do.call(rbind, Xin)
  Xout <- c(y[n:(n - lag_y + 1L)], factors[n, , drop = TRUE])

  Xin  <- as.matrix(Xin)
  yin  <- as.numeric(yin)
  Xout <- as.numeric(Xout)

  keep <- apply(Xin, 1, function(z) all(is.finite(z))) & is.finite(yin)
  Xin <- Xin[keep, , drop = FALSE]
  yin <- yin[keep]

  list(Xin = Xin, yin = yin, Xout = Xout, K_actual = K_actual)
}

fit_one_tvp <- function(df, variable, horizon, lag_y, K, lambda_use, panel_names) {
  des <- make_design(
    df = df,
    variable = variable,
    horizon = horizon,
    lag_y = lag_y,
    K = K,
    corr_cutoff = corr_cutoff,
    panel_names = panel_names
  )

  if (nrow(des$Xin) < min_eval_n) {
    stop(sprintf("Poucas observações úteis: %d", nrow(des$Xin)))
  }

  out <- tvp.ridge(
    X                 = des$Xin,
    Y                 = matrix(des$yin, ncol = 1),
    lambda.candidates = lambda_use,
    oosX              = des$Xout,
    kfold             = 5L,
    CV.2SRR           = TRUE,
    CV.plot           = FALSE,
    sig.eps.param     = 0.75,
    sig.u.param       = 0.75
  )

  list(
    pred   = as.numeric(out$forecast),
    lambda = out$lambdas,
    K_used = des$K_actual
  )
}

select_hyperparams <- function(df, variable, horizon, lag_grid, K_grid, lambda_grid, panel_names) {
  combos <- expand.grid(lag_y = lag_grid, K = K_grid)
  scores <- rep(NA_real_, nrow(combos))
  lambdas <- vector("list", nrow(combos))

  n <- nrow(df)
  val_size <- max(12L, horizon + 6L)
  split_pt <- n - val_size
  if (split_pt <= (max(lag_grid) + horizon + 10L)) {
    stop("Janela insuficiente para seleção de hiperparâmetros.")
  }

  train_df <- df[1:split_pt, , drop = FALSE]
  valid_df <- df[1:(split_pt + val_size), , drop = FALSE]

  for (ii in seq_len(nrow(combos))) {
    lag_y <- combos$lag_y[ii]
    K     <- combos$K[ii]

    tmp <- tryCatch({
      base_fit <- fit_one_tvp(train_df, variable, horizon, lag_y, K, lambda_grid, panel_names)
      roll_preds <- c()
      roll_actual <- c()

      for (tt in seq(split_pt, split_pt + val_size - 1L)) {
        sub_df <- valid_df[1:tt, , drop = FALSE]
        fit_tt <- fit_one_tvp(sub_df, variable, horizon, lag_y, K, base_fit$lambda, panel_names)
        idx_act <- tt + horizon
        if (idx_act <= nrow(valid_df)) {
          roll_preds <- c(roll_preds, fit_tt$pred)
          roll_actual <- c(roll_actual, as.numeric(valid_df[idx_act, variable]))
        }
      }

      err <- roll_preds - roll_actual
      list(score = mean(err^2, na.rm = TRUE), lambda = base_fit$lambda)
    }, error = function(e) {
      NULL
    })

    if (!is.null(tmp)) {
      scores[ii] <- tmp$score
      lambdas[[ii]] <- tmp$lambda
    }
  }

  if (all(is.na(scores))) {
    stop(sprintf("Falha ao selecionar hiperparâmetros para %s h=%d", variable, horizon))
  }

  best_i <- which.min(scores)
  list(
    lag_y = combos$lag_y[best_i],
    K = combos$K[best_i],
    lambda = lambdas[[best_i]],
    score = scores[best_i],
    grid = cbind(combos, msfe = scores)
  )
}

# ============================================================
# FORECAST LOOP
# ============================================================

set.seed(1234)
forecasts_list <- list()
failure_log <- list()
hyperparam_log <- list()

panel_names <- colnames(panel_df)

for (v in variable_list) {
  for (h in horizon_list) {
    log_line(sprintf("\n[%s] var=%-12s | h=%d", model_name, v, h))

    best <- NULL
    model_list <- vector("list", nwindows)
    fail_count <- 0L

    for (i in seq_len(nwindows)) {
      train_end <- start_window - h - 1L + i

      if (train_end < min_train_n) {
        model_list[[i]] <- list(pred = NA_real_)
        fail_count <- fail_count + 1L
        failure_log[[length(failure_log) + 1L]] <- data.frame(
          variable = v, horizon = h, window = i,
          train_end = train_end, type = "short_train",
          message = "train_end below minimum"
        )
        next
      }

      Df <- if (window_type == "expanding") {
        df_scaled[1:train_end, , drop = FALSE]
      } else {
        df_scaled[max(1L, train_end - 119L):train_end, , drop = FALSE]
      }

      reopt <- is.null(best) || (i == 1L) || (i %% reopt_every == 1L)
      if (reopt) {
        log_line(sprintf("  selecting hyperparameters (window %d/%d)", i, nwindows))
      }

      model <- tryCatch({
        if (reopt) {
          best <- select_hyperparams(
            df = Df,
            variable = v,
            horizon = h,
            lag_grid = lag_orders,
            K_grid = K_grid,
            lambda_grid = lambda_grid,
            panel_names = panel_names
          )

          hyperparam_log[[length(hyperparam_log) + 1L]] <- data.frame(
            variable = v,
            horizon = h,
            window = i,
            lag_y = best$lag_y,
            K = best$K,
            tuning_msfe = best$score
          )

          log_line(sprintf("    chosen lag=%d | K=%d | tuning MSFE=%.6f",
                           best$lag_y, best$K, best$score))
        }

        fit <- fit_one_tvp(
          df = Df,
          variable = v,
          horizon = h,
          lag_y = best$lag_y,
          K = best$K,
          lambda_use = best$lambda,
          panel_names = panel_names
        )

        list(pred = fit$pred)
      }, error = function(e) {
        fail_count <<- fail_count + 1L
        failure_log[[length(failure_log) + 1L]] <<- data.frame(
          variable = v, horizon = h, window = i,
          train_end = train_end, type = "model_error",
          message = conditionMessage(e)
        )
        log_line(sprintf("  [ERROR] var=%s h=%d win=%d: %s", v, h, i, conditionMessage(e)))
        list(pred = NA_real_)
      })

      model_list[[i]] <- list(pred = model$pred)
    }

    preds <- vapply(model_list, function(x) x$pred, numeric(1L))
    actual_idx_end <- min(end_window, start_window + nwindows - 1L)
    actual_vec <- as.numeric(df_scaled[start_window:actual_idx_end, v])

    forecasts_list[[length(forecasts_list) + 1L]] <- list(
      variable = v,
      horizon  = h,
      model    = model_name,
      pred     = preds,
      actual   = actual_vec,
      failures = fail_count
    )

    log_line(sprintf("  -> %d forecasts | NAs: %d | failures: %d",
                     length(preds), sum(is.na(preds)), fail_count))
  }
}

# ============================================================
# SAVE FORECASTS
# ============================================================

saveRDS(forecasts_list, file = file.path(output_run, paste0(model_name, ".rds")))

# ============================================================
# CHECKUPS
# ============================================================

results <- do.call(rbind, lapply(forecasts_list, function(f) {
  n <- min(length(f$pred), length(f$actual))
  keep <- !is.na(f$pred[1:n]) & !is.na(f$actual[1:n])
  if (sum(keep) == 0) {
    return(data.frame(
      variable = f$variable,
      horizon  = f$horizon,
      n_eval   = 0,
      fail_windows = f$failures,
      MSFE     = NA_real_,
      RMSFE    = NA_real_,
      MAE      = NA_real_
    ))
  }
  err <- f$pred[1:n][keep] - f$actual[1:n][keep]
  data.frame(
    variable = f$variable,
    horizon  = f$horizon,
    n_eval   = sum(keep),
    fail_windows = f$failures,
    MSFE     = mean(err^2),
    RMSFE    = sqrt(mean(err^2)),
    MAE      = mean(abs(err))
  )
}))

rw_results <- do.call(rbind, lapply(forecasts_list, function(f) {
  act <- f$actual
  h   <- f$horizon
  n   <- length(act)
  if (n <= h) return(NULL)
  rw_pred <- act[1:(n - h)]
  rw_act  <- act[(1 + h):n]
  keep <- !is.na(rw_pred) & !is.na(rw_act)
  data.frame(
    variable = f$variable,
    horizon  = h,
    MSFE_RW  = if (sum(keep) > 0) mean((rw_pred[keep] - rw_act[keep])^2) else NA_real_
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
                      n_eval   = sum(keep),
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
    n_eval   = sum(keep),
    DM_stat  = if (!is.null(dm_test)) as.numeric(dm_test$statistic) else NA_real_,
    p_value  = if (!is.null(dm_test)) dm_test$p.value else NA_real_,
    MSFE_rel = mean(e1[keep]^2) / mean(e2[keep]^2)
  )
}))

failure_table <- if (length(failure_log) > 0) do.call(rbind, failure_log) else data.frame(
  variable = character(), horizon = integer(), window = integer(),
  train_end = integer(), type = character(), message = character()
)

failure_summary <- if (nrow(failure_table) > 0) {
  aggregate(window ~ variable + horizon + type, data = failure_table, FUN = length) |>
    dplyr::rename(n_fail = window)
} else {
  data.frame(variable = character(), horizon = integer(), type = character(), n_fail = integer())
}

hyperparam_table <- if (length(hyperparam_log) > 0) do.call(rbind, hyperparam_log) else data.frame(
  variable = character(), horizon = integer(), window = integer(),
  lag_y = integer(), K = integer(), tuning_msfe = numeric()
)

results_summary <- list(
  fit_stats = results,
  msfe_rel  = results_rel,
  dm_vs_rw  = dm_results,
  failure_table = failure_table,
  failure_summary = failure_summary,
  hyperparams = hyperparam_table
)

saveRDS(results_summary, file = file.path(paths$results, "results_summary.rds"))
write.csv(results,           file = file.path(paths$results, "TAB01_fit_stats.csv"), row.names = FALSE)
write.csv(results_rel,       file = file.path(paths$results, "TAB02_msfe_relativo.csv"), row.names = FALSE)
write.csv(dm_results,        file = file.path(paths$results, "TAB03_diebold_mariano.csv"), row.names = FALSE)
write.csv(failure_table,     file = file.path(paths$results, "TAB04_failure_log.csv"), row.names = FALSE)
write.csv(failure_summary,   file = file.path(paths$results, "TAB05_failure_summary.csv"), row.names = FALSE)
write.csv(hyperparam_table,  file = file.path(paths$results, "TAB06_hyperparams.csv"), row.names = FALSE)

log_line("\n=== RESUMO GERAL ===")
log_line(sprintf("Modelo      : %s", model_name))
log_line(sprintf("Janela      : %s", window_type))
log_line(sprintf("Variáveis   : %s", paste(variable_list, collapse = ", ")))
log_line(sprintf("Horizontes  : %s", paste(horizon_list, collapse = ", ")))
log_line(sprintf("Combinações : %d", length(forecasts_list)))
log_line(sprintf("OOS windows : %d", nwindows))
log_line(sprintf("Período OOS : índices %d a %d", start_window, end_window))
log_line(sprintf("Arquivo de log: %s", log_file))
log_line("\n02_03_forecast.R finalizado.")

