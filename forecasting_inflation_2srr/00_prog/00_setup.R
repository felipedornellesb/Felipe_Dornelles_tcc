# ==============================================================================
# 00_setup.R
# ==============================================================================

cat("== 00_setup.R ==\n\n")

ROOT <- getwd()
stopifnot("Run from project root (folder containing 00_prog/)." =
            dir.exists(file.path(ROOT, "00_prog")))

DIR_DATA        <- file.path(ROOT, "10_data")
DIR_COULOMBE    <- file.path(ROOT, "20_tools", "21_coulombe")
DIR_MEDEIROS    <- file.path(ROOT, "20_tools", "22_medeiros")
DIR_ADAPTED     <- file.path(ROOT, "20_tools", "23_adapted")
DIR_FORECASTS   <- file.path(ROOT, "30_output", "forecasts")
DIR_BETAS       <- file.path(ROOT, "30_output", "betas")
DIR_CHECKPOINTS <- file.path(ROOT, "30_output", "checkpoints")
DIR_TABLES      <- file.path(ROOT, "40_results", "tables")
DIR_FIGURES     <- file.path(ROOT, "40_results", "figures")

for (d in c(DIR_DATA, DIR_COULOMBE, DIR_MEDEIROS, DIR_ADAPTED,
            DIR_FORECASTS, DIR_BETAS, DIR_CHECKPOINTS, DIR_TABLES, DIR_FIGURES))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# -- Packages ------------------------------------------------------------------
pkgs <- c("glmnet", "pracma", "randomForest",
          "forecast", "lmtest", "sandwich", "ggplot2", "reshape2", "xtable")
new_pkgs <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(new_pkgs)) install.packages(new_pkgs, repos = "https://cran.r-project.org")
suppressPackageStartupMessages(invisible(lapply(pkgs, library, character.only = TRUE)))

if (!"HDeconometrics" %in% installed.packages()[, "Package"]) {
  if (!"devtools" %in% installed.packages()[, "Package"])
    install.packages("devtools", repos = "https://cran.r-project.org")
  devtools::install_github("gabrielrvsc/HDeconometrics")
}
suppressPackageStartupMessages(library(HDeconometrics))

for (gp in c("rugarch", "fGarch"))
  if (!gp %in% installed.packages()[, "Package"])
    tryCatch(install.packages(gp, repos = "https://cran.r-project.org"),
             error = function(e) NULL)

# ==============================================================================
# CRITICAL: factor() PCA override BEFORE sourcing Coulombe code.
# EM_sw.R calls factor(X, n_fac=n). Without this, base::factor errors.
# ==============================================================================
base_factor_backup <- base::factor
factor <- function(X, n_fac = NULL, ...) {
  if (!is.null(n_fac)) {
    X[is.na(X)] <- 0
    pc  <- prcomp(X, center = TRUE, scale. = TRUE)
    nf  <- min(n_fac, ncol(pc$x))
    fac <- pc$x[, 1:nf, drop = FALSE]
    lam <- pc$rotation[, 1:nf, drop = FALSE]
    mse <- mean((X - fac %*% t(lam))^2, na.rm = TRUE)
    return(list(factors = fac, lambda = lam, mse = mse))
  }
  base_factor_backup(X, ...)
}
cat("[OK] factor() PCA defined\n")

# -- Download Coulombe functions (ALL 10 files) --------------------------------
base_tools <- paste0("https://raw.githubusercontent.com/hugocout/",
                      "Replication-codes-for-Time-Varying-Parameters-",
                      "as-Ridge-Regressions/main/Empirical/20_tools")
base_func  <- paste0(base_tools, "/functions")

coulombe_files <- list(
  list(url = paste0(base_tools, "/EM_sw.R"),              f = "EM_sw.R"),
  list(url = paste0(base_tools, "/ICp2.R"),               f = "ICp2.R"),
  list(url = paste0(base_func, "/zfun_v190304.R"),         f = "zfun_v190304.R"),
  list(url = paste0(base_func, "/dualGRRmdA_v190215.R"),  f = "dualGRRmdA_v190215.R"),
  list(url = paste0(base_func, "/CVGSBHK_v181127.R"),     f = "CVGSBHK_v181127.R"),
  list(url = paste0(base_func, "/CVKFMV_v190214.R"),      f = "CVKFMV_v190214.R"),
  list(url = paste0(base_func, "/TVPRR_v181111.R"),       f = "TVPRR_v181111.R"),
  list(url = paste0(base_func, "/TVPRRcosso_v181120.R"),  f = "TVPRRcosso_v181120.R"),
  list(url = paste0(base_func, "/fastZrot_v181125.R"),     f = "fastZrot_v181125.R"),
  list(url = paste0(base_func, "/Xgenerators_v190127.R"), f = "Xgenerators_v190127.R")
)

cat("\nCoulombe functions:\n")
for (cf in coulombe_files) {
  dest <- file.path(DIR_COULOMBE, cf$f)
  if (!file.exists(dest)) {
    tryCatch({
      download.file(cf$url, dest, quiet = TRUE)
      cat(sprintf("  downloaded %-30s (%d bytes)\n", cf$f, file.size(dest)))
    }, error = function(e) cat(sprintf("  FAILED %-30s %s\n", cf$f, e$message)))
  } else {
    cat(sprintf("  found      %-30s (%d bytes)\n", cf$f, file.size(dest)))
  }
}

# -- Download Medeiros functions -----------------------------------------------
medeiros_base <- paste0("https://raw.githubusercontent.com/gabrielrvsc/",
                         "ForecastingInflation/main")
cat("\nMedeiros functions:\n")
for (f in c("functions/functions.R", "functions/rolling_window.R")) {
  dest <- file.path(DIR_MEDEIROS, basename(f))
  if (!file.exists(dest)) {
    tryCatch({
      download.file(paste0(medeiros_base, "/", f), dest, quiet = TRUE)
      cat("  downloaded", basename(f), "\n")
    }, error = function(e) cat("  FAILED", basename(f), "\n"))
  } else {
    cat("  found", basename(f), "\n")
  }
}

# -- data.rda ------------------------------------------------------------------
if (!file.exists(file.path(DIR_DATA, "data.rda")))
  cat("\n*** Place data.rda in 10_data/ before proceeding. ***\n\n")

# -- Source Coulombe (ORDER MATTERS) -------------------------------------------
cat("\nSourcing Coulombe:\n")
for (f in c("EM_sw.R", "ICp2.R", "zfun_v190304.R",
            "dualGRRmdA_v190215.R", "CVGSBHK_v181127.R", "CVKFMV_v190214.R",
            "TVPRR_v181111.R", "TVPRRcosso_v181120.R",
            "fastZrot_v181125.R", "Xgenerators_v190127.R")) {
  fp <- file.path(DIR_COULOMBE, f)
  if (file.exists(fp)) {
    tryCatch({ source(fp, local = FALSE); cat(sprintf("  sourced %s\n", f)) },
             error = function(e) cat(sprintf("  WARN %s: %s\n", f, e$message)))
  }
}

# Verify critical functions
cat("\nFunction check:\n")
for (fn in c("dualGRR", "TVPRR_cosso", "TVPRR", "cvgs.bhk2015",
             "Zfun", "EM_sw", "make_reg_matrix", "fastZrot"))
  cat(sprintf("  %-20s %s\n", fn,
              ifelse(exists(fn, mode = "function"), "OK", "*** MISSING ***")))

# -- Source Medeiros -----------------------------------------------------------
cat("\nSourcing Medeiros:\n")
for (f in c("functions.R", "rolling_window.R")) {
  fp <- file.path(DIR_MEDEIROS, f)
  if (file.exists(fp)) { source(fp, local = FALSE); cat("  sourced", f, "\n") }
}

# -- Source adapted ------------------------------------------------------------
for (fp in list.files(DIR_ADAPTED, "\\.R$", full.names = TRUE)) {
  source(fp, local = FALSE); cat("  sourced", basename(fp), "\n")
}

cat("\n== Setup complete ==\n")
