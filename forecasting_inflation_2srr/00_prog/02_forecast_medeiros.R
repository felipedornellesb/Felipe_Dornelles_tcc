# ==============================================================================
# 02_forecast_medeiros.R
# ==============================================================================
cat("== 02_forecast_medeiros.R ==\n\n")
source("00_prog/00_setup.R")

load(file.path(DIR_DATA, "data.rda"))

variable <- "CPIAUCSL"
horizons <- c(1, 3, 6, 12)
maxh     <- 12

# CRITICAL: remove date column (embed() and model functions require numeric)
dates_col <- data$date
data$date <- NULL
rownames(data) <- as.character(dates_col)
cat("Data prepared:", nrow(data), "obs x", ncol(data), "numeric vars\n")

# List available model functions
all_models <- list(
  list(name = "LASSO",    fn = "runlasso"),
  list(name = "Ridge",    fn = "runridge"),
  list(name = "ElNET",    fn = "runelnet"),
  list(name = "AdaLASSO", fn = "runadaptlasso"),
  list(name = "AdaElNET", fn = "runadaptelnet"),
  list(name = "RF",       fn = "runrf"),
  list(name = "Bagging",  fn = "runbag"),
  list(name = "Factor",   fn = "runfactor"),
  list(name = "T.Factor", fn = "runtargetfactor"),
  list(name = "CSR",      fn = "runcsr"),
  list(name = "AR",       fn = "runar"),
  list(name = "AR_BIC",   fn = "runarbic")
)

cat("\nAvailable functions: ")
avail <- sapply(all_models, function(m) exists(m$fn, mode = "function"))
cat(paste(sapply(all_models[avail], "[[", "name"), collapse = ", "), "\n\n")

for (m in all_models) {
  out_path <- file.path(DIR_FORECASTS, paste0(m$name, ".rda"))
  if (file.exists(out_path)) { cat(sprintf("  %-12s exists\n", m$name)); next }
  if (!exists(m$fn, mode = "function")) { next }

  cat(sprintf("  %-12s running...\n", m$name))
  t0 <- Sys.time()
  
  forecasts_mat <- matrix(NA_real_, nwindows, maxh)
  betas_bundle <- list()
  
  for (h in horizons) {
    cat(sprintf("    h=%2d...", h))
    
    # Direct forecasting: create cumulative target for horizon h
    data_h <- data
    y_h <- as.numeric(stats::filter(data[[variable]], rep(1, h), sides = 1))
    if (h > 1) y_h[1:(h-1)] <- y_h[h]
    data_h[[variable]] <- y_h
    
    tryCatch({
      result <- rolling_window(get(m$fn), data_h, nwindows, h, variable)
      forecasts_mat[, h] <- result$forecast
      betas_bundle[[paste0("h", h)]] <- result$outputs
      cat(" done.\n")
    }, error = function(e) cat(sprintf(" FAILED: %s\n", e$message)))
  }
  
  tryCatch({
    forecasts <- forecasts_mat
    colnames(forecasts) <- paste0("h", 1:maxh)
    save(forecasts, file = out_path)
    
    # Save the betas and lambdas bundle
    if (length(betas_bundle) > 0) {
      save(betas_bundle, file = file.path(DIR_BETAS, paste0("betas_", m$name, ".rda")))
    }
    
    cat(sprintf("  %-12s %.1f min\n", m$name, difftime(Sys.time(), t0, units = "mins")))
  }, error = function(e) cat(sprintf(" FAILED to save: %s\n", e$message)))
}

cat("== done ==\n")
