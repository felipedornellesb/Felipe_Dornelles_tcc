# ==============================================================================
# 04_analysis.R
#
# Unified analysis: RMSFE, DM, Clark-West, beta tables, parsimony,
# sub-periods, CSSED, figures, LaTeX.
# ==============================================================================
cat("== 04_analysis.R ==\n\n")
source("00_prog/00_setup.R")

load(file.path(DIR_FORECASTS, "yout.rda"))
load(file.path(DIR_FORECASTS, "rw.rda"))
load(file.path(DIR_DATA, "data.rda"))

horizons <- c(1, 3, 6, 12)
all_h    <- 1:12
nw       <- nrow(yout)
oos_dates <- tail(data$date, nw)

# ==============================================================================
# 1. LOAD ALL FORECASTS
# ==============================================================================
cat("--- Loading forecasts ---\n")

med_names <- c("Ridge", "LASSO", "ElNET", "AdaLASSO", "AdaElNET",
               "RF", "Bagging", "Factor", "T.Factor", "CSR", "AR", "AR_BIC")
all_fc <- list()
for (mn in c(med_names, "2SRR")) {
  fp <- file.path(DIR_FORECASTS, paste0(mn, ".rda"))
  if (file.exists(fp)) {
    env <- new.env(); load(fp, envir = env)
    all_fc[[mn]] <- as.matrix(get(ls(env)[1], envir = env))
    cat(sprintf("  %-12s %d x %d\n", mn, nrow(all_fc[[mn]]), ncol(all_fc[[mn]])))
  }
}

if (!is.null(all_fc[["2SRR"]]) && !is.null(all_fc[["Ridge"]])) {
  all_fc[["Half-Half"]] <- 0.5 * all_fc[["2SRR"]] + 0.5 * all_fc[["Ridge"]]
  cat(sprintf("  %-12s %d x %d\n", "Half-Half", nrow(all_fc[["Half-Half"]]), ncol(all_fc[["Half-Half"]])))
}

# ==============================================================================
# 2. RMSFE TABLE
# ==============================================================================
cat("\n--- RMSFE ---\n")

rmsfe_fn <- function(y, f) {
  ok <- complete.cases(y, f); if (sum(ok) < 5) return(NA)
  sqrt(mean((y[ok] - f[ok])^2))
}

rmsfe_rw <- sapply(all_h, function(h) rmsfe_fn(yout[, h], rw[, h]))

rmsfe_table <- data.frame(model = "RW", stringsAsFactors = FALSE)
for (h in all_h) rmsfe_table[[paste0("h", h)]] <- 1.0

for (mn in names(all_fc)) {
  row <- data.frame(model = mn, stringsAsFactors = FALSE)
  fc_mat <- all_fc[[mn]]
  for (h in all_h) {
    val <- if (h <= ncol(fc_mat)) rmsfe_fn(yout[, h], fc_mat[, h]) / rmsfe_rw[h] else NA
    row[[paste0("h", h)]] <- round(val, 4)
  }
  rmsfe_table <- rbind(rmsfe_table, row)
}

cols_show <- c("model", paste0("h", horizons))
print(rmsfe_table[, intersect(cols_show, names(rmsfe_table))], row.names = FALSE)
write.csv(rmsfe_table, file.path(DIR_TABLES, "rmsfe_all.csv"), row.names = FALSE)

# ==============================================================================
# 3. DM AND CLARK-WEST TESTS
# ==============================================================================
cat("\n--- Statistical tests ---\n")

tests <- data.frame()
if (!is.null(all_fc[["2SRR"]])) {
  for (h in horizons) {
    if (h > ncol(all_fc[["2SRR"]])) next
    y_h  <- yout[, h]
    fc_s <- all_fc[["2SRR"]][, h]
    ok   <- complete.cases(y_h, fc_s, rw[, h])
    if (sum(ok) < 20) next

    # DM: 2SRR vs RW
    dm <- tryCatch(dm.test(ts((y_h[ok]-fc_s[ok])^2), ts((y_h[ok]-rw[ok,h])^2),
                            alternative = "less", h = h),
                    error = function(e) list(statistic = NA, p.value = NA))
    tests <- rbind(tests, data.frame(h = h, test = "DM", comp = "2SRR_vs_RW",
                                      stat = as.numeric(dm$statistic), p = dm$p.value))

    # DM: 2SRR vs Ridge (if available)
    if (!is.null(all_fc[["Ridge"]]) && h <= ncol(all_fc[["Ridge"]])) {
      fc_r <- all_fc[["Ridge"]][, h]
      ok2  <- ok & complete.cases(fc_r)
      dm2 <- tryCatch(dm.test(ts((y_h[ok2]-fc_s[ok2])^2), ts((y_h[ok2]-fc_r[ok2])^2),
                               alternative = "less", h = h),
                       error = function(e) list(statistic = NA, p.value = NA))
      tests <- rbind(tests, data.frame(h = h, test = "DM", comp = "2SRR_vs_Ridge",
                                        stat = as.numeric(dm2$statistic), p = dm2$p.value))
      cw <- clark_west(y_h[ok2], fc_r[ok2], fc_s[ok2])
      tests <- rbind(tests, data.frame(h = h, test = "CW", comp = "2SRR_vs_Ridge",
                                        stat = cw$stat, p = cw$pvalue))
    }
  }
}
if (nrow(tests) > 0) {
  print(tests, row.names = FALSE, digits = 3)
  write.csv(tests, file.path(DIR_TABLES, "statistical_tests.csv"), row.names = FALSE)
}

# ==============================================================================
# 4. BETA ANALYSIS
# ==============================================================================
cat("\n--- Beta analysis ---\n")

bpath <- file.path(DIR_BETAS, "betas_2SRR.rda")
if (file.exists(bpath)) {
  load(bpath)  # betas_bundle

  for (hname in names(betas_bundle)) {
    bb <- betas_bundle[[hname]]
    bt <- bb$betas_tvp
    om <- bb$omega
    la <- bb$lambda

    # Extract last-period beta from each window
    K <- NULL
    for (b in bt) if (!is.null(b) && is.matrix(b)) { K <- ncol(b); break }
    if (is.null(K)) { cat(sprintf("  %s: no betas\n", hname)); next }

    beta_mat <- matrix(NA, length(bt), K)
    for (i in seq_along(bt)) {
      b <- bt[[i]]
      if (!is.null(b) && is.matrix(b) && ncol(b) == K) beta_mat[i, ] <- b[nrow(b), ]
    }
    colnames(beta_mat) <- paste0("X", 1:K)

    # Table A: beta trajectory
    df_traj <- data.frame(window = 1:nrow(beta_mat), beta_mat, check.names = FALSE)
    write.csv(df_traj, file.path(DIR_TABLES, paste0("betas_trajectory_", hname, ".csv")),
              row.names = FALSE)

    # Table B: sigma2_u ranking (last window with data)
    last_om <- NULL
    for (i in rev(seq_along(om))) if (!is.null(om[[i]])) { last_om <- om[[i]]; break }
    if (!is.null(last_om) && length(last_om) == K) {
      sigma_df <- data.frame(predictor = paste0("X", 1:K), sigma2_u = last_om,
                              pct = round(100 * last_om / sum(last_om), 2))
      sigma_df <- sigma_df[order(-sigma_df$sigma2_u), ]
      write.csv(sigma_df, file.path(DIR_TABLES, paste0("sigma2u_", hname, ".csv")),
                row.names = FALSE)
    }

    # Table C: TVP vs constant comparison
    comp <- data.frame(
      predictor  = paste0("X", 1:K),
      tvp_mean   = colMeans(beta_mat, na.rm = TRUE),
      tvp_sd     = apply(beta_mat, 2, sd, na.rm = TRUE),
      norm_delta = apply(beta_mat, 2, function(x) sqrt(sum(diff(x)^2, na.rm = TRUE))),
      tvp_cv     = apply(beta_mat, 2, function(x) {
        m <- mean(x, na.rm = TRUE)
        ifelse(abs(m) > 1e-10, sd(x, na.rm = TRUE) / abs(m), NA) })
    )
    comp <- comp[order(-comp$norm_delta), ]
    write.csv(comp, file.path(DIR_TABLES, paste0("tvp_vs_constant_", hname, ".csv")),
              row.names = FALSE)

    # Lambda trajectory
    write.csv(data.frame(window = seq_along(la), lambda = la),
              file.path(DIR_TABLES, paste0("lambda_", hname, ".csv")), row.names = FALSE)

    med_cv <- median(comp$tvp_cv, na.rm = TRUE)
    cat(sprintf("  %s: K=%d, median CV=%.3f %s\n", hname, K, med_cv,
                ifelse(med_cv > 0.5, "(substantial variation)", "(modest)")))
  }
}

# ==============================================================================
# 5. SUB-PERIODS
# ==============================================================================
cat("\n--- Sub-periods ---\n")

periods <- list(
  pre_GFC  = oos_dates < as.Date("2007-12-01"),
  GFC      = oos_dates >= as.Date("2007-12-01") & oos_dates <= as.Date("2009-06-01"),
  post_GFC = oos_dates > as.Date("2009-06-01") & oos_dates < as.Date("2020-02-01"),
  COVID    = oos_dates >= as.Date("2020-02-01") & oos_dates <= as.Date("2021-12-01"),
  post_COVID = oos_dates > as.Date("2021-12-01")
)

sub_table <- data.frame()
for (pn in names(periods)) {
  idx <- which(periods[[pn]])
  if (length(idx) < 5) next
  for (mn in names(all_fc)) {
    for (h in horizons) {
      if (h > ncol(all_fc[[mn]])) next
      rw_sub <- rmsfe_fn(yout[idx, h], rw[idx, h])
      if (!is.finite(rw_sub) || rw_sub < 1e-10) next
      sub_table <- rbind(sub_table, data.frame(
        period = pn, model = mn, h = h,
        rmsfe = round(rmsfe_fn(yout[idx, h], all_fc[[mn]][idx, h]) / rw_sub, 4),
        n = length(idx), stringsAsFactors = FALSE))
    }
  }
}
if (nrow(sub_table) > 0)
  write.csv(sub_table, file.path(DIR_TABLES, "subperiod_rmsfe.csv"), row.names = FALSE)

# ==============================================================================
# 6. CSSED
# ==============================================================================
cat("\n--- CSSED ---\n")
if (!is.null(all_fc[["2SRR"]])) {
  for (h in horizons) {
    if (h > ncol(all_fc[["2SRR"]])) next
    y_h <- yout[, h]; fc_s <- all_fc[["2SRR"]][, h]
    ok <- complete.cases(y_h, rw[, h], fc_s)
    cssed <- cumsum((y_h[ok] - rw[ok, h])^2 - (y_h[ok] - fc_s[ok])^2)
    write.csv(data.frame(window = which(ok), date = oos_dates[ok], cssed = cssed),
              file.path(DIR_TABLES, sprintf("cssed_h%02d.csv", h)), row.names = FALSE)
  }
}

# ==============================================================================
# 7. FIGURES
# ==============================================================================
cat("\n--- Figures ---\n")

# RMSFE bar chart
tryCatch({
  key <- rmsfe_table[rmsfe_table$model %in% c("Ridge","LASSO","RF","2SRR","RW"), ]
  if (nrow(key) > 0) {
    dl <- reshape2::melt(key, id.vars = "model",
                          measure.vars = paste0("h", horizons),
                          variable.name = "horizon", value.name = "rmsfe")
    dl$rmsfe <- as.numeric(dl$rmsfe)
    p <- ggplot(dl, aes(x = horizon, y = rmsfe, fill = model)) +
      geom_bar(stat = "identity", position = "dodge", width = 0.7) +
      geom_hline(yintercept = 1, linetype = "dashed") +
      labs(title = "RMSFE relative to RW", x = "Horizon", y = "RMSFE / RW") +
      theme_minimal(base_size = 11) + theme(legend.position = "bottom")
    ggsave(file.path(DIR_FIGURES, "rmsfe_comparison.pdf"), p, width = 10, height = 6)
    cat("  rmsfe_comparison.pdf\n")
  }
}, error = function(e) NULL)

# LaTeX
tryCatch({
  xt <- xtable(rmsfe_table[, intersect(c("model", paste0("h", horizons)),
                                         names(rmsfe_table))],
               caption = "RMSFE relative to Random Walk", label = "tab:rmsfe", digits = 4)
  print(xt, file = file.path(DIR_TABLES, "rmsfe.tex"), include.rownames = FALSE)
  cat("  rmsfe.tex\n")
}, error = function(e) NULL)

cat("\n== 04_analysis.R complete ==\n")
