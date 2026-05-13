# ==============================================================================
# 02_forecast_medeiros.R
# ==============================================================================
cat("== 02_forecast_medeiros.R ==\n\n")
source("00_prog/00_setup.R")

load(file.path(DIR_DATA, "data.rda"))

variable <- "CPIAUCSL"
nwindows <- 180
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
  for (h in 1:maxh) {
    cat(sprintf("    h=%2d...", h))
    tryCatch({
      result <- rolling_window(get(m$fn), data, nwindows, h, variable)
      forecasts_mat[, h] <- result$forecast
      cat(" done.\n")
    }, error = function(e) cat(sprintf(" FAILED: %s\n", e$message)))
  }
  
  tryCatch({
    if (exists("accumulate_model", mode = "function"))
      forecasts <- accumulate_model(forecasts_mat)
    else
      forecasts <- forecasts_mat
      
    save(forecasts, file = out_path)
    cat(sprintf("  %-12s %.1f min\n", m$name, difftime(Sys.time(), t0, units = "mins")))
  }, error = function(e) cat(sprintf(" FAILED to accumulate/save: %s\n", e$message)))
}

cat("== done ==\n")
