# ============================================================
# 02_forecast.R
#
# POOS para modelos Ridge / TVP-2SRR (2SRR) / MSRRs / MSRRd
# Lê dados de 10_data/data_MM_DD_YYYY/ (pasta mais recente)
# Salva TVPfcst_V*_H*_M*.RData em 30_output/outputs_MM_DD_YYYY/
#
# Especificações (M):
#   1 = AR puro
#   2 = AR + fatores PCA
#   3 = AR + outros alvos
#   4 = AR + fatores PCA + outros alvos
#
# Estimadores (dim 4 do array forecast):
#   [1] Ridge   [2] 2SRR   [3] MSRRs   [4] MSRRd
# ============================================================

rm(list = ls())

# ============================================================
# 0. PATHS
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
for (p in c(paths$output, output_run, paths$results))
  if (!dir.exists(p)) dir.create(p, recursive = TRUE)

cat(sprintf("Output folder: %s\n", output_run))

# ============================================================
# MODO TESTE — mude para FALSE para rodar completo
# ============================================================
TEST_MODE    <- TRUE
TEST_V       <- 1L    # 1=PIB, 2=IPCA, 3=SELIC, 4=CAMBIO
TEST_H       <- 1L
TEST_M       <- 1L
TEST_JANELAS <- 24L   # janelas OOS no modo teste

# ============================================================
# 1. PACOTES
# ============================================================

myPKGs <- c("pracma", "fGarch", "matrixcalc", "dplyr",
            "glmnet", "zoo", "timeSeries")
need   <- myPKGs[!myPKGs %in% names(installed.packages()[, "Package"])]
if (length(need) > 0)
  install.packages(need, repos = "http://cran.us.r-project.org")
invisible(lapply(myPKGs, library, character.only = TRUE))

# ============================================================
# 2. FUNÇÕES
# ============================================================

source(file.path(paths$tools,     "EM_sw.R"))
source(file.path(paths$tools,     "ICp2.R"))
source(file.path(paths$tools,     "Xgenerators_v190127.R"))
source(file.path(paths$functions, "dualGRRmdA_v190215.R"))
source(file.path(paths$functions, "CVGSBHK_v181127.R"))
source(file.path(paths$functions, "zfun_v190304.R"))
source(file.path(paths$functions, "factor.R"))
source(file.path(paths$functions, "TVPRRcosso_v181120.R"))
source(file.path(paths$functions, "TVPRRcossoF_v190125.R"))
source(file.path(paths$functions, "TVPRR_v181111.R"))
source(file.path(paths$functions, "fastZrot_v181125.R"))
source(file.path(paths$functions, "CVKFMV_v190214.R"))
source(file.path(paths$functions, "TVPRR_VARF_v190304.R"))
source(file.path(paths$functions, "MV2SRR_v221103.R"))

# ============================================================
# 3. CARREGA DADOS (pasta mais recente em 10_data/)
# ============================================================

data_dirs <- sort(
  list.dirs(paths$data, recursive = FALSE, full.names = TRUE)
)
data_dirs <- data_dirs[grepl("data_\\d{2}_\\d{2}_\\d{4}$",
                              basename(data_dirs))]
if (length(data_dirs) == 0)
  stop("Nenhuma pasta data_MM_DD_YYYY encontrada. Rode 01_data_prep.R primeiro.")

run_folder <- tail(data_dirs, 1)
cat(sprintf("Dados de: %s\n", run_folder))

load(file.path(run_folder, "df_model.rda"))
load(file.path(run_folder, "df_targets.rda"))
load(file.path(run_folder, "df_panel_pca.rda"))
load(file.path(run_folder, "targets_br.rda"))
load(file.path(run_folder, "all_options.rda"))

target_names <- unname(unlist(targets_br))
szv          <- length(target_names)
bigt         <- nrow(df_model)
dates_vec    <- df_model$date

cat(sprintf("T = %d  |  %s a %s\n", bigt,
            format(min(dates_vec), "%b/%Y"),
            format(max(dates_vec), "%b/%Y")))
cat(sprintf("Alvos: %s\n", paste(target_names, collapse = " | ")))

mat_y   <- as.matrix(df_targets[, target_names])
pan_cols <- setdiff(names(df_panel_pca), "date")
mat_pan  <- as.matrix(df_panel_pca[, pan_cols])

cat(sprintf("mat_pan: %d x %d\n", nrow(mat_pan), ncol(mat_pan)))

# Checagem de NAs
na_check <- colSums(is.na(mat_y))
if (any(na_check > 0)) {
  cat("ATENÇÃO — NAs em mat_y:\n"); print(na_check[na_check > 0])
}

# ============================================================
# 4. PARÂMETROS
# ============================================================

tau        <- 156L
lags_y     <- 4L
n_factors  <- min(6L, ncol(mat_pan) - 1L)
lambda_vec <- exp(pracma::linspace(-6, 20, n = 15))
reopt_freq <- 12L
scree_ts   <- 0.05
maxf       <- min(3L, n_factors)
alpha_m4   <- 0.15
sv.param   <- 0
silenceplz <- 1

horizons_all <- sort(unique(all_options$H))
max_H        <- max(horizons_all)

cat(sprintf("tau=%d (%s)  lags_y=%d  n_fac=%d  max_H=%d\n",
            tau, format(dates_vec[tau], "%b/%Y"),
            lags_y, n_factors, max_H))

# ============================================================
# 5. EM_sw — IMPUTAÇÃO DO PAINEL
# ============================================================

col_var       <- apply(mat_pan, 2, function(x) var(x, na.rm = TRUE))
mat_pan_clean <- mat_pan[, col_var > 1e-10, drop = FALSE]
cat(sprintf("Colunas com variância > 0: %d\n", ncol(mat_pan_clean)))

if (ncol(mat_pan_clean) >= (n_factors + 1L)) {
  panel_full <- tryCatch({
    EM_sw(data = mat_pan_clean, n = n_factors, it_max = 500)$data
  }, error = function(e) {
    warning(sprintf("EM_sw falhou: %s — usando mat_pan_clean.", e$message))
    mat_pan_clean
  })
} else {
  panel_full <- mat_pan_clean
  n_factors  <- max(1L, ncol(mat_pan_clean) - 1L)
  maxf       <- min(maxf, n_factors)
}

if (any(is.na(panel_full))) panel_full[is.na(panel_full)] <- 0
cat(sprintf("panel_full: %d x %d\n", nrow(panel_full), ncol(panel_full)))

# ============================================================
# 6. FUNÇÕES DO LOOP POOS
# ============================================================

build_regressors <- function(y_full, mat_y_full, panel_imput,
                              v_idx, h, spec, lags_y, n_factors,
                              train_end) {
  y_train <- y_full[1:train_end]
  if (sum(is.finite(y_train)) < (lags_y + h + 10L)) return(NULL)

  lag_y <- embed(y_train, lags_y)
  T_lag <- nrow(lag_y)
  y_dep <- y_full[(lags_y + h):train_end]
  T_eff <- length(y_dep)

  if (T_eff > T_lag || T_eff < 20L) return(NULL)

  X_ar  <- lag_y[1:T_eff, , drop = FALSE]
  Xo_ar <- as.numeric(lag_y[T_lag, ])

  X_fac <- NULL; Xo_fac <- NULL
  if (spec %in% c(2L, 4L) && ncol(panel_imput) >= 2L) {
    pan_w  <- panel_imput[1:train_end, , drop = FALSE]
    col_ok <- apply(pan_w, 2, function(x)
      all(is.finite(x)) && var(x, na.rm = TRUE) > 1e-10)
    pan_w  <- pan_w[, col_ok, drop = FALSE]
    if (ncol(pan_w) >= 2L) {
      K_f  <- min(n_factors, ncol(pan_w), nrow(pan_w) - 1L)
      pc   <- prcomp(pan_w, center = TRUE, scale. = TRUE)
      fac  <- pc$x[, seq_len(K_f), drop = FALSE]
      fl   <- embed(fac, lags_y)
      T_f  <- nrow(fl)
      if (T_f >= T_eff) {
        X_fac  <- fl[1:T_eff, , drop = FALSE]
        Xo_fac <- as.numeric(fl[T_f, ])
      }
    }
  }

  X_tgt <- NULL; Xo_tgt <- NULL
  if (spec %in% c(3L, 4L)) {
    tp <- list(); to <- list()
    for (vi in setdiff(seq_len(ncol(mat_y_full)), v_idx)) {
      yi <- mat_y_full[1:train_end, vi]
      if (any(!is.finite(yi))) next
      yl  <- embed(yi, lags_y)
      T_y <- nrow(yl)
      if (T_y < T_eff) next
      tp[[length(tp) + 1L]] <- yl[1:T_eff, , drop = FALSE]
      to[[length(to) + 1L]] <- as.numeric(yl[T_y, ])
    }
    if (length(tp) > 0) {
      X_tgt  <- do.call(cbind, tp)
      Xo_tgt <- unlist(to)
    }
  }

  Xin  <- X_ar;  Xout <- Xo_ar
  if (!is.null(X_fac)) { Xin <- cbind(Xin, X_fac); Xout <- c(Xout, Xo_fac) }
  if (!is.null(X_tgt)) { Xin <- cbind(Xin, X_tgt); Xout <- c(Xout, Xo_tgt) }

  col_ok <- apply(Xin, 2, function(x)
    all(is.finite(x)) && var(x, na.rm = TRUE) > 1e-10)
  Xin  <- Xin[, col_ok, drop = FALSE]
  Xout <- Xout[col_ok]

  row_ok <- is.finite(y_dep) & apply(Xin, 1, function(r) all(is.finite(r)))
  Xin    <- Xin[row_ok, , drop = FALSE]
  y_dep  <- y_dep[row_ok]

  if (nrow(Xin) < 20L || ncol(Xin) < 1L) return(NULL)
  list(Xin = Xin, yin = y_dep, Xout = Xout)
}

OF <- function(pred, y_ref, go.to) {
  tol <- 2
  mx  <- (max(y_ref) - mean(y_ref)) * tol + mean(y_ref)
  mn  <- mean(y_ref) - (mean(y_ref) - min(y_ref)) * tol
  ifelse(pred > mx | pred < mn, go.to, pred)
}

run_estimators <- function(Xin, yin, Xout, lambda_vec,
                            lv_fixed, reoptimize,
                            silenceplz, scree_ts, maxf,
                            alpha_m4, sv.param) {
  preds  <- rep(NA_real_, 4L)
  lv_new <- lv_fixed

  # [1] Ridge
  pred_lin <- NA_real_
  tryCatch({
    CV        <- cv.glmnet(x = Xin, y = yin, family = "gaussian",
                           alpha = 0, nfolds = 5)
    mdl       <- glmnet(x = Xin, y = yin, family = "gaussian",
                        alpha = 0, lambda = CV$lambda.min)
    pred_lin  <- as.numeric(predict(mdl, newx = matrix(Xout, nrow = 1)))
    preds[1L] <- pred_lin
    lv_new    <- CV$lambda.min
  }, error = function(e) message(sprintf("  [Ridge] %s", e$message)))

  lv2 <- if (!is.null(lv_new) && is.finite(lv_new)) lv_new else 1.0

  # [2] 2SRR
  tryCatch({
    aa2 <- TVPRR_cosso(
      y = yin, X = Xin, lambdavec = lambda_vec,
      sweigths = 1, type = 2, alpha = 0.01,
      silent = silenceplz, kfold = 5,
      lambda2 = lv2, tol = 1e-6, maxit = 10, oosX = Xout
    )
    p2 <- as.numeric(aa2$fcast)
    preds[2L] <- OF(p2, yin, if (is.finite(pred_lin)) pred_lin else p2)
  }, error = function(e) message(sprintf("  [2SRR] %s", e$message)))

  # [3] MSRRs
  tryCatch({
    aa3 <- TVPRR(
      y = yin, X = Xin, lambdavec = lambda_vec,
      sweigths = 1, type = 3, alpha = 0.001,
      silent = silenceplz, kfold = 5,
      lambda2 = lv2, tol = 1e-5, maxit = 15, oosX = Xout
    )
    p3 <- as.numeric(aa3$fcast)
    preds[3L] <- OF(p3, yin, if (is.finite(pred_lin)) pred_lin else p3)
  }, error = function(e) message(sprintf("  [MSRRs] %s", e$message)))

  # [4] MSRRd
  tryCatch({
    aa4 <- TVPRR_VARF(
      Y = yin, X = Xin, orthoFac = TRUE,
      lambdavec = lambda_vec[4:15],
      sweigths = 1, type = 2, fp.model = 1,
      sv.param = sv.param, alpha = alpha_m4,
      silent = silenceplz, kfold = 5,
      lambda2 = lv2, max.step.cv = 8,
      adaptive = 1, aparam = -0.5,
      tol = 1e-10, maxit = 20, lambdabooster = 1,
      var.share = scree_ts, override = maxf,
      id = 1, oosX = Xout
    )
    p4 <- as.numeric(aa4$fcast)
    preds[4L] <- OF(p4, yin, if (is.finite(pred_lin)) pred_lin else p4)
  }, error = function(e) message(sprintf("  [MSRRd] %s", e$message)))

  list(preds = preds, lv_fixed = lv_new)
}

# ============================================================
# 7. COMBINAÇÕES
# ============================================================

set.seed(1234L)

if (TEST_MODE) {
  combos    <- data.frame(V = TEST_V, H = TEST_H, M = TEST_M)
  tau_test  <- max(tau, bigt - TEST_JANELAS)
  cat(sprintf("\n*** MODO TESTE: V=%d (%s) H=%d M=%d | janelas OOS: %d ***\n",
              TEST_V, target_names[TEST_V], TEST_H, TEST_M,
              bigt - tau_test))
  fname_test <- file.path(output_run,
    sprintf("TVPfcst_V%d_H%d_M%d.RData", TEST_V, TEST_H, TEST_M))
  if (file.exists(fname_test)) { file.remove(fname_test); cat("Arquivo anterior removido.\n") }
} else {
  combos   <- unique(all_options[, c("V", "H", "M")])
  tau_test <- tau
}

cat(sprintf("%d combinação(ões) a processar.\n", nrow(combos)))

# ============================================================
# 8. LOOP POOS
# ============================================================

t0_total <- proc.time()

for (row_i in seq_len(nrow(combos))) {
  v_idx <- combos$V[row_i]
  h     <- combos$H[row_i]
  m     <- combos$M[row_i]

  fname <- file.path(output_run,
                     sprintf("TVPfcst_V%d_H%d_M%d.RData", v_idx, h, m))

  if (file.exists(fname) && !TEST_MODE) {
    cat(sprintf("  [%d/%d] V%d H%d M%d — já existe.\n",
                row_i, nrow(combos), v_idx, h, m))
    next
  }

  cat(sprintf("\n[%d/%d] V%d (%s) | H=%d | M=%d\n",
              row_i, nrow(combos), v_idx, target_names[v_idx], h, m))

  forecast <- array(NA_real_, dim = c(bigt, max_H, szv, 4L))
  y_full   <- mat_y[, v_idx]
  lv_fixed <- NULL
  n_valid  <- 0L
  n_skip   <- 0L
  t0_combo <- proc.time()
  t_start  <- if (TEST_MODE) tau_test else tau

  for (t in t_start:bigt) {
    train_end <- t - h
    if (train_end < 30L || !is.finite(y_full[t])) { n_skip <- n_skip + 1L; next }

    reopt <- is.null(lv_fixed) || ((t - t_start) %% reopt_freq == 0L)

    prep <- tryCatch(
      build_regressors(y_full, mat_y, panel_full,
                       v_idx, h, m, lags_y, n_factors, train_end),
      error = function(e) { message(sprintf("  build t=%d: %s", t, e$message)); NULL }
    )
    if (is.null(prep)) { n_skip <- n_skip + 1L; next }

    res <- run_estimators(prep$Xin, prep$yin, prep$Xout,
                          lambda_vec, lv_fixed, reopt,
                          silenceplz, scree_ts, maxf, alpha_m4, sv.param)

    if (!is.null(res$lv_fixed)) lv_fixed <- res$lv_fixed
    forecast[t, h, v_idx, ] <- res$preds
    n_valid <- n_valid + 1L

    freq_log <- if (TEST_MODE) 1L else 24L
    if (n_valid %% freq_log == 0L) {
      el <- (proc.time() - t0_combo)["elapsed"]
      cat(sprintf("  t=%d (%s)  ok=%d skip=%d  [%.0fs]\n",
                  t, format(dates_vec[t], "%b/%Y"), n_valid, n_skip, el))
    }
  }

  n_fill <- sum(!is.na(forecast[, h, v_idx, ]))
  cat(sprintf("  -> %d previsões salvas\n", n_fill))
  if (n_fill == 0L)
    warning(sprintf("NENHUMA previsão válida para V%d H%d M%d!", v_idx, h, m))

  save(forecast, file = fname)
  cat(sprintf("  Salvo: %s\n", basename(fname)))
}

cat(sprintf("\n=== POOS finalizado em %.1f min ===\n",
            (proc.time() - t0_total)["elapsed"] / 60))

# ============================================================
# 9. DIAGNÓSTICO (MODO TESTE)
# ============================================================

if (TEST_MODE) {
  fname_diag <- file.path(output_run,
    sprintf("TVPfcst_V%d_H%d_M%d.RData", TEST_V, TEST_H, TEST_M))

  cat("\n========== DIAGNÓSTICO ==========\n")
  if (!file.exists(fname_diag)) {
    cat("!!! Arquivo de previsão NÃO criado — cheque erros acima !!!\n")
  } else {
    load(fname_diag)
    preds_mat <- forecast[, TEST_H, TEST_V, ]
    janelas   <- which(!is.na(preds_mat[, 1]))
    cat(sprintf("Janelas com Ridge válido: %d (esperado: %d)\n",
                length(janelas), bigt - tau_test + 1L))
    if (length(janelas) > 0) {
      y_real <- mat_y[, TEST_V]
      df_val <- data.frame(
        date = dates_vec[janelas],
        real = y_real[janelas],
        ridge = preds_mat[janelas, 1],
        srr2  = preds_mat[janelas, 2],
        msrrs = preds_mat[janelas, 3],
        msrrd = preds_mat[janelas, 4]
      )
      msfe  <- function(p, r) mean((p - r)^2, na.rm = TRUE)
      base  <- msfe(df_val$ridge, df_val$real)
      cat(sprintf("  Ridge  : 1.000 (MSFE=%.6f)\n", base))
      cat(sprintf("  2SRR   : %.3f\n", msfe(df_val$srr2,  df_val$real) / base))
      cat(sprintf("  MSRRs  : %.3f\n", msfe(df_val$msrrs, df_val$real) / base))
      cat(sprintf("  MSRRd  : %.3f\n", msfe(df_val$msrrd, df_val$real) / base))
    }
  }
  cat("=================================\n")
}
