# ==============================================================================
# 05_descriptive_analysis.R
#
# Descriptive analysis of 2SRR betas on the FULL sample, pre-forecasting.
# Compares TVP betas with static Ridge betas as requested by the professor.
# ==============================================================================
cat("== 05_descriptive_analysis.R ==\n\n")
source("00_prog/00_setup.R")

load(file.path(DIR_DATA, "data.rda"))

variable <- "CPIAUCSL"
horizon  <- 1
n_lags   <- 4
kfold    <- 5

# Data preparation
dates_col <- data$date
data$date <- NULL
rownames(data) <- as.character(dates_col)

df_num <- as.matrix(data[, sapply(data, is.numeric)])
y_raw <- as.numeric(data[[variable]])

# Embed with lags
X_embed <- embed(df_num, n_lags)
n_align <- nrow(X_embed) - horizon
X_in    <- X_embed[1:n_align, , drop = FALSE]
y_in    <- y_raw[(n_lags + horizon):(n_lags + horizon + n_align - 1)]

# Clean up variables
cv <- apply(X_in, 2, var, na.rm = TRUE)
good <- which(is.finite(cv) & cv > 1e-10)
X_in <- X_in[, good, drop = FALSE]

ok <- complete.cases(X_in, y_in)
X_in <- X_in[ok, , drop = FALSE]; y_in <- y_in[ok]

# Standardize
X_mu <- colMeans(X_in); X_sd <- apply(X_in, 2, sd); X_sd[X_sd < 1e-10] <- 1
X_sc <- sweep(sweep(X_in, 2, X_mu, "-"), 2, X_sd, "/")

# ==============================================================================
# 1. RUN 2SRR (FULL SAMPLE)
# ==============================================================================
cat("\nRunning 2SRR on full sample...\n")
t0 <- Sys.time()
fit_2srr <- tvp_2srr_standalone(X_sc, y_in, kfold = kfold)
cat(sprintf("Done in %.1f mins.\n", difftime(Sys.time(), t0, units = "mins")))

betas_tvp <- fit_2srr$beta
lam2_opt  <- fit_2srr$lambda

# ==============================================================================
# 2. RUN RIDGE (FULL SAMPLE)
# ==============================================================================
cat("Running Ridge on full sample...\n")
library(glmnet)
cv_ridge <- cv.glmnet(X_sc, y_in, alpha = 0, nfolds = kfold)
beta_ridge <- as.numeric(coef(cv_ridge, s = "lambda.min"))[-1] # remove intercept

# ==============================================================================
# 3. ANALYSIS AND EXPORT
# ==============================================================================

# Create comparison table
K <- ncol(X_sc)
var_names <- paste0("X", 1:K)
colnames(betas_tvp) <- var_names

comp <- data.frame(
  predictor    = var_names,
  ridge_beta   = beta_ridge,
  tvp_mean     = colMeans(betas_tvp, na.rm = TRUE),
  tvp_sd       = apply(betas_tvp, 2, sd, na.rm = TRUE),
  tvp_cv       = apply(betas_tvp, 2, function(x) sd(x)/abs(mean(x))),
  sigma2_u     = fit_2srr$omega
)

comp <- comp[order(-comp$sigma2_u), ]
write.csv(comp, file.path(DIR_TABLES, "pre_forecast_betas_comp.csv"), row.names = FALSE)
cat("\nSaved pre_forecast_betas_comp.csv\n")

# Top 10 varying predictors
top10 <- head(comp$predictor, 10)

# Save beta trajectory for top 10
traj_top10 <- data.frame(t = 1:nrow(betas_tvp), betas_tvp[, top10, drop=FALSE])
write.csv(traj_top10, file.path(DIR_TABLES, "pre_forecast_traj_top10.csv"), row.names = FALSE)

# Count 'active' parameters (parsimony)
med_sigma2u <- median(fit_2srr$omega)
n_active <- sum(fit_2srr$omega > med_sigma2u)
pct_active <- 100 * n_active / K

cat(sprintf("\nParsimony Analysis:\n"))
cat(sprintf("  Optimal Lambda2: %.4f\n", lam2_opt))
cat(sprintf("  Active TVP Parameters (> median sigma2u): %d out of %d (%.1f%%)\n", n_active, K, pct_active))

cat("\n== 05_descriptive_analysis.R complete ==\n")
