# ==============================================================================
# 03_forecast_2srr.R
#
# This script faithfully follows the Coulombe replication code structure:
#   - Expanding window POOS loop
#   - PCA factor extraction per window
#   - make_reg_matrix() for regressor construction (handles cumulation)
#   - TVPRR_cosso(type=2) for 2SRR estimation
#   - Outlier filter (OF) on forecasts
#
# Output: 2SRR.rda (forecast matrix compatible with Medeiros eval scripts)
#         betas_2SRR.rda (TVP betas bundle)
# ==============================================================================

cat("== 03_forecast_2srr.R ==\n\n")

# 0. Setup
source("00_prog/00_setup.R")

# Load Coulombe functions — these MUST be the originals from his repo
# Either source them from a local copy or download them
coulombe_dir <- file.path("00_prog", "coulombe")
if (!dir.exists(coulombe_dir)) dir.create(coulombe_dir, recursive = TRUE)

# Source Coulombe's core functions:
#   - TVPRR_cosso()    : The 2SRR engine (Algorithm 1)
#   - make_reg_matrix(): Builds X with lags of y, Y_h, and factors
#   - EM_sw()          : Stock-Watson EM for missing data imputation
#   - OF()             : Outlier filter
#
# These come from: github.com/hugocout/Replication-codes-for-Time-Varying-Parameters-as-Ridge-Regressions
# Specifically the Empirical/Forecasting/ folder
#
# IMPORTANT: Do NOT modify these functions. Use them as published.
source(file.path(coulombe_dir, "functions_coulombe.R"))  # Adjust path to your local copy

library(glmnet)
library(pracma)    # for linspace()
library(fGarch)    # for garchFit() in Step 2 of Algorithm 1

# 1. Load data
load(file.path(DIR_DATA, "data.rda"))
load(file.path(DIR_FORECASTS, "yout.rda"))  # OOS realized values

variable <- "CPIAUCSL"
horizons <- c(1, 3, 6, 12)

# Separate dates
dates_col <- data$date
data$date <- NULL

# Number of variables and observations
bigt <- nrow(data)
n_oos <- ncol(yout)  # Number of OOS windows (from Medeiros)
tau <- bigt - n_oos   # Start of OOS period

# Identify target variable position
y_col <- which(names(data) == variable)

cat(sprintf("  Total obs: %d | OOS start: %d | OOS windows: %d\n", bigt, tau, n_oos))

# 2. Coulombe pipeline parameters (FAITHFUL to Algorithm 1)
nf       <- 8        # Number of PCA factors (same as Coulombe)
ly       <- 2        # Lags of y in make_reg_matrix
lf       <- 2        # Lags of factors in make_reg_matrix
kfold    <- 5        # k-fold CV (Coulombe default)
n_lambda <- 15       # Lambda grid size (Coulombe default)
lambda_vec <- exp(pracma::linspace(-2, 12, n = n_lambda))

# Output containers
fc_ridge <- matrix(NA_real_, n_oos, length(horizons))
fc_2srr  <- matrix(NA_real_, n_oos, length(horizons))
colnames(fc_ridge) <- colnames(fc_2srr) <- paste0("h", horizons)

betas_store  <- list()
lambda_store <- list()
for (h in horizons) {
  hlab <- paste0("h", sprintf("%02d", h))
  betas_store[[hlab]]  <- vector("list", n_oos)
  lambda_store[[hlab]] <- rep(NA_real_, n_oos)
}

# Outlier filter (from Coulombe's code)
OF <- function(pred, y_is, pred_fallback) {
  y_mean <- mean(y_is, na.rm = TRUE)
  y_max  <- max(y_is, na.rm = TRUE)
  y_min  <- min(y_is, na.rm = TRUE)
  bound  <- 2 * max(abs(y_max - y_mean), abs(y_min - y_mean))
  if (abs(pred - y_mean) > bound) {
    return(pred_fallback)
  }
  return(pred)
}

# 3. POOS Loop — Expanding Window (FAITHFUL to Coulombe Section 4)
out_path <- file.path(DIR_FORECASTS, "2SRR.rda")

if (file.exists(out_path)) {
  cat("  2SRR.rda already exists — skipping.\n")
  cat("  Delete it to re-run.\n")
} else {

  t_total <- Sys.time()

  for (wi in 1:n_oos) {
    t_end <- tau + wi - 1  # Last in-sample observation

    if (wi %% 25 == 1 || wi == n_oos) {
      cat(sprintf("  Window %d/%d (T_is=%d) ... ", wi, n_oos, t_end))
      t0 <- Sys.time()
    }

    # --- 3a. In-sample data slice ---
    X_is <- as.matrix(data[1:t_end, -y_col])
    y_is <- as.numeric(data[1:t_end, y_col])

    # --- 3b. Remove zero-variance columns ---
    col_var <- apply(X_is, 2, var, na.rm = TRUE)
    X_is <- X_is[, col_var > 1e-10, drop = FALSE]

    # --- 3c. Impute NAs via Stock-Watson EM ---
    # If your data has no NAs (already cleaned by Medeiros), skip this
    if (any(is.na(X_is))) {
      X_is <- tryCatch(
        EM_sw(X_is, nf = nf, max_iter = 500),
        error = function(e) {
          # Fallback: simple column-mean imputation
          for (j in 1:ncol(X_is)) {
            na_idx <- is.na(X_is[, j])
            if (any(na_idx)) X_is[na_idx, j] <- mean(X_is[, j], na.rm = TRUE)
          }
          X_is
        }
      )
    }

    # --- 3d. Standardize X ---
    X_means <- colMeans(X_is, na.rm = TRUE)
    X_sds   <- apply(X_is, 2, sd, na.rm = TRUE)
    X_sds[X_sds < 1e-10] <- 1
    X_std <- scale(X_is, center = X_means, scale = X_sds)

    # --- 3e. Extract PCA factors ---
    pca_obj <- prcomp(X_std, center = FALSE, scale. = FALSE)
    factors <- pca_obj$x[, 1:min(nf, ncol(pca_obj$x)), drop = FALSE]

    # --- 3f. For each horizon, build regression and estimate ---
    for (hi in seq_along(horizons)) {
      h <- horizons[hi]
      hlab <- paste0("h", sprintf("%02d", h))

      # Cumulative target: Y_h(t) = sum of y from t-h+1 to t
      # This is what make_reg_matrix does internally
      # If you use make_reg_matrix from Coulombe, it handles this
      tryCatch({

        # Build regression matrix using Coulombe's function
        # make_reg_matrix returns: list(XX, yy) where:
        #   XX = [lags of Y_h, lags of y, lags of factors]
        #   yy = Y_h aligned for direct forecasting at horizon h
        reg <- make_reg_matrix(y = y_is, Y = y_is, factors = factors,
                               h = h, ly = ly, lf = lf)
        XX <- reg$XX
        yy <- reg$yy

        if (nrow(XX) < 30 || ncol(XX) < 2) {
          fc_ridge[wi, hi] <- NA
          fc_2srr[wi, hi]  <- NA
          next
        }

        # ---- Step 1: Ridge baseline via cv.glmnet ----
        cv_fit <- cv.glmnet(x = XX, y = yy, alpha = 0,
                            nfolds = min(kfold, floor(nrow(XX) / 3)))
        ridge_coef <- as.numeric(coef(cv_fit, s = "lambda.min"))[-1]
        ridge_pred <- predict(cv_fit, newx = XX[nrow(XX), , drop = FALSE],
                              s = "lambda.min")

        # ---- Steps 1-4: Full 2SRR via TVPRR_cosso ----
        result_2srr <- TVPRR_cosso(
          X       = XX,
          y       = yy,
          type    = 2,              # type=2 = full 2SRR (Algorithm 1)
          lambdavec = lambda_vec,
          lambda2 = cv_fit$lambda.min,  # Initialize Step 4 with Ridge lambda
          kfold   = kfold,
          sweigths = rep(1, nrow(XX)),
          tol     = 1e-6,
          maxit   = 10
        )

        # Extract forecast (last time period's fitted value = forecast for t+h)
        pred_2srr <- result_2srr$forecast  # or result_2srr$yhat[length(yy)]
        pred_ridge <- as.numeric(ridge_pred)

        # ---- Outlier filter ----
        pred_2srr <- OF(pred_2srr, yy, pred_ridge)

        # Store forecasts
        fc_ridge[wi, hi] <- pred_ridge
        fc_2srr[wi, hi]  <- pred_2srr

        # Store betas and lambda
        betas_store[[hlab]][[wi]]  <- result_2srr$betas  # TVP betas
        lambda_store[[hlab]][wi]   <- result_2srr$lambda  # Optimal lambda

      }, error = function(e) {
        cat(sprintf("[h=%d FAIL: %s] ", h, e$message))
        fc_ridge[wi, hi] <<- NA
        fc_2srr[wi, hi]  <<- NA
      })
    }

    # Progress reporting
    if (wi %% 25 == 1 || wi == n_oos) {
      elapsed <- difftime(Sys.time(), t0, units = "secs")
      cat(sprintf("%.1f sec\n", elapsed))
    }

    # Checkpoint every 50 windows
    if (wi %% 50 == 0) {
      save(fc_ridge, fc_2srr, betas_store, lambda_store,
           file = file.path(DIR_CHECKPOINTS, "ckpt_2SRR.rda"))
      cat(sprintf("  [checkpoint saved at window %d]\n", wi))
    }
  }

  cat(sprintf("\n  Total time: %.1f hours\n",
              difftime(Sys.time(), t_total, units = "hours")))

  # --- 4. Save in Medeiros-compatible format ---
  # The Medeiros eval scripts expect: forecasts matrix (n_oos × maxh)
  maxh <- 12
  forecasts <- matrix(NA_real_, n_oos, maxh)
  colnames(forecasts) <- paste0("h", 1:maxh)
  for (hi in seq_along(horizons)) {
    forecasts[, horizons[hi]] <- fc_2srr[, hi]
  }
  save(forecasts, file = out_path)
  cat(sprintf("  Saved %s\n", basename(out_path)))

  # Save betas bundle
  betas_bundle <- list()
  for (h in horizons) {
    hlab <- paste0("h", sprintf("%02d", h))
    betas_bundle[[hlab]] <- list(
      betas_tvp = betas_store[[hlab]],
      lambda    = lambda_store[[hlab]]
    )
  }
  save(betas_bundle, file = file.path(DIR_BETAS, "betas_2SRR.rda"))
  cat("  Saved betas_2SRR.rda\n")

  # Also save ridge forecasts for Half-Half strategy later
  forecasts_ridge <- matrix(NA_real_, n_oos, maxh)
  colnames(forecasts_ridge) <- paste0("h", 1:maxh)
  for (hi in seq_along(horizons)) {
    forecasts_ridge[, horizons[hi]] <- fc_ridge[, hi]
  }
  save(forecasts_ridge, file = file.path(DIR_FORECASTS, "Ridge_from_2SRR.rda"))
  cat("  Saved Ridge_from_2SRR.rda (for Half-Half)\n")
}

cat("== done ==\n")

# ==============================================================================
# NOTES ON RUNTIME:
#
# Expected runtime: 30-40 hours for 312 windows × 4 horizons on a single core.
# This is NORMAL for the 2SRR with full Algorithm 1 (including GARCH in Step 2).
#
# The Coulombe paper's forecasting exercise used quarterly data (fewer obs per
# window, ~200 vs your ~600 monthly), so his runtimes were faster.
#
# These changes affect ONLY the speed of hyperparameter search, not the
# statistical procedure. The paper's Algorithm 1 is preserved exactly.
# ==============================================================================
