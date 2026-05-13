# ==============================================================================
# 01_data_prep.R
# ==============================================================================
cat("== 01_data_prep.R ==\n\n")
source("00_prog/00_setup.R")

load(file.path(DIR_DATA, "data.rda"))
cat("Data:", nrow(data), "obs x", ncol(data), "vars\n")

variable <- "CPIAUCSL"
nwindows <- 180
maxh     <- 12
y   <- data[[variable]]
n   <- length(y)
tau <- n - nwindows

yout <- matrix(NA, nwindows, maxh)
rw   <- matrix(NA, nwindows, maxh)
for (h in 1:maxh) for (i in 1:nwindows) {
  t_end <- tau + i
  if ((t_end + h) <= n) yout[i, h] <- sum(y[(t_end + 1):(t_end + h)])
  rw[i, h] <- h * y[t_end]
}
colnames(yout) <- colnames(rw) <- paste0("h", 1:maxh)

save(yout, file = file.path(DIR_FORECASTS, "yout.rda"))
save(rw,   file = file.path(DIR_FORECASTS, "rw.rda"))
cat("Saved yout and rw:", nwindows, "x", maxh, "\n== done ==\n")
