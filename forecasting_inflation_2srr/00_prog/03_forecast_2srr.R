# ==============================================================================
# 03_forecast_2srr.R
#
# Runs the 2SRR via rolling_window() from Medeiros.
# Saves 2SRR.rda (forecast matrix) and betas_2SRR.rda (beta bundle).
# ==============================================================================
cat("== 03_forecast_2srr.R ==\n\n")
source("00_prog/00_setup.R")

load(file.path(DIR_DATA, "data.rda"))

variable <- "CPIAUCSL"
nwindows <- 180
maxh     <- 12

# Remove date column
dates_col <- data$date
data$date <- NULL
rownames(data) <- as.character(dates_col)

out_path <- file.path(DIR_FORECASTS, "2SRR.rda")

if (file.exists(out_path)) {
  cat("  2SRR.rda exists, skipping.\n")
} else {
  model_list <- vector("list", maxh)

  for (h in 1:maxh) {
    cat(sprintf("  h=%2d: %d windows...", h, nwindows))
    t0 <- Sys.time()

    tryCatch({
      model_list[[h]] <- rolling_window(
        fn       = run2srr,
        df       = data,
        nwindow  = nwindows,
        horizon  = h,
        variable = variable,
        kfold    = 5,
        n_lags   = 4
      )
      cat(sprintf(" %.1f min\n", difftime(Sys.time(), t0, units = "mins")))
    }, error = function(e) {
      cat(sprintf(" FAILED: %s\n", e$message))
      model_list[[h]] <<- list(forecast = rep(NA_real_, nwindows),
                                outputs = vector("list", nwindows))
    })

    # Checkpoint
    save(model_list, file = file.path(DIR_CHECKPOINTS, "ckpt_2SRR.rda"))
  }

  # Assemble forecast matrix
  forecasts <- matrix(NA_real_, nwindows, maxh)
  for (h in 1:maxh) {
    fc <- model_list[[h]]$forecast
    forecasts[1:min(length(fc), nwindows), h] <- fc[1:min(length(fc), nwindows)]
  }
  colnames(forecasts) <- paste0("h", 1:maxh)

  if (exists("accumulate_model", mode = "function"))
    forecasts <- accumulate_model(forecasts)

  save(forecasts, file = out_path)
  cat(sprintf("  Saved %s\n", basename(out_path)))

  # Extract betas
  betas_bundle <- list()
  for (h in 1:maxh) {
    ml <- model_list[[h]]
    if (is.null(ml$outputs)) next
    betas_h  <- vector("list", length(ml$outputs))
    omega_h  <- vector("list", length(ml$outputs))
    lambda_h <- rep(NA_real_, length(ml$outputs))
    for (i in seq_along(ml$outputs)) {
      out <- ml$outputs[[i]]
      if (is.null(out)) next
      betas_h[[i]]  <- out$betas_tvp
      omega_h[[i]]  <- out$omega
      lambda_h[i]   <- ifelse(is.null(out$lambda), NA, out$lambda)
    }
    betas_bundle[[paste0("h", h)]] <- list(
      betas_tvp = betas_h, omega = omega_h, lambda = lambda_h)
  }
  save(betas_bundle, file = file.path(DIR_BETAS, "betas_2SRR.rda"))
  cat("  Saved betas_2SRR.rda\n")
}

cat("== done ==\n")
