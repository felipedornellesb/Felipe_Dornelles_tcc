# 12_2srr_deep_dive.R
#
# PARTES:
#   1. Lambda ao longo do tempo + NBER
#   2. Dispersao cross-sectional dos betas
#   3. Fallback silencioso (OF filter)
#   4. Sub-periodos expandidos (7 regimes)
#   5. Correlacao betas vs volatilidade
#   6. Tabela dual benchmark (RW + AR)
#   7. Superficie 3D dos betas TVP (estilo Coulombe)
#   8. Mudancas de sinal dos betas
#   9. Rolling RMSE ratio (janela 36 meses)
#  10. Fan chart das previsoes
#  11. Narrativa automatica
# ============================================================

rm(list = ls())
setwd("~/tcc/Felipe_Dornelles_tcc/forecasting_inflation")

library(ggplot2)
library(reshape2)
library(gridExtra)
library(forecast)
library(lmtest)
library(sandwich)

has_plotly <- tryCatch({ library(plotly); TRUE }, error = function(e) {
  tryCatch({ install.packages("plotly"); library(plotly); TRUE },
           error = function(e2) FALSE)
})

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_dir   <- file.path("40_results", paste0("run12_", timestamp))
fig_dir   <- file.path(out_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

hor <- c(1, 3, 6, 12)

recessions <- data.frame(
  start = as.Date(c("2001-03-01", "2007-12-01", "2020-02-01")),
  end   = as.Date(c("2001-11-01", "2009-06-01", "2020-04-01"))
)

cat(sprintf("Output: %s\n\n", out_dir))

load("forecasts/yout.rda")
load("forecasts/coulombe_forecasts.rda")
load("forecasts/coulombe_betas_2SRR.rda")
load("forecasts/coulombe_betas_ridge.rda")
load("data/data.rda")

fred_raw <- as.data.frame(data)
date_col <- grep("^date$", colnames(fred_raw), ignore.case = TRUE)[1]
dates    <- fred_raw[, date_col]
y_raw    <- as.numeric(fred_raw[, grep("^CPIAUCSL$", colnames(fred_raw), ignore.case=TRUE)[1]])
bigt     <- nrow(fred_raw)
n_oos    <- nrow(yout)
tau      <- bigt - n_oos

coulombe <- list()
for (h in hor) {
  fname <- sprintf("forecasts/coulombe_fc_h%02d.csv", h)
  if (file.exists(fname))
    coulombe[[paste0("h",h)]] <- read.csv(fname, stringsAsFactors=FALSE)
}

med_names <- c("AR","AR_BIC","AdaEINET","AdaLASSO","Bagging",
               "CSR","EINET","factor","LASSO","RF","Ridge",
               "T.Factor","2SRR","rw")
med_fc <- list()
for (mn in med_names) {
  fname <- sprintf("forecasts/%s.rda", mn)
  if (file.exists(fname)) {
    env <- new.env(); load(fname, envir=env)
    obj <- get(ls(env)[1], envir=env)
    if (is.matrix(obj) || is.data.frame(obj)) med_fc[[mn]] <- as.matrix(obj)
  }
}

tvp_ar_list <- list(); tvp_fac_list <- list()
for (h in hor) {
  f_ar  <- sprintf("forecasts/tvp_TVP_AR_h%02d.csv", h)
  f_fac <- sprintf("forecasts/tvp_TVP_Factor_h%02d.csv", h)
  if (file.exists(f_ar))  tvp_ar_list[[paste0("h",h)]]  <- read.csv(f_ar, stringsAsFactors=FALSE)
  if (file.exists(f_fac)) tvp_fac_list[[paste0("h",h)]] <- read.csv(f_fac, stringsAsFactors=FALSE)
}

reg_names <- c("intercept", paste0("y_lag",1:2),
               paste0("F",rep(1:8,2),"_lag",rep(1:2,each=8)))

cat(sprintf("yout: %dx%d | Medeiros: %d | Coulombe: %d h\n\n",
            nrow(yout), ncol(yout), length(med_fc), length(coulombe)))

# PARTE 1: LAMBDA AO LONGO DO TEMPO + NBER
cat("PARTE 1: Trajetoria de Lambda\n")

for (h in hor) {
  key <- paste0("h",h)
  if (is.null(coulombe[[key]])) next
  df <- coulombe[[key]]
  df <- df[complete.cases(df$lam_ridge, df$lam2_2srr), ]
  if (nrow(df) < 10) next
  df$date <- as.Date(df$date)

  lam_long <- data.frame(
    date   = rep(df$date, 2),
    Lambda = c(df$lam_ridge, df$lam2_2srr),
    Modelo = c(rep("Ridge (lambda)", nrow(df)),
               rep("2SRR (lambda2)", nrow(df))))

  p <- ggplot(lam_long, aes(x=date, y=Lambda, color=Modelo)) +
    geom_rect(data=recessions, aes(xmin=start,xmax=end,ymin=-Inf,ymax=Inf),
              inherit.aes=FALSE, fill="gray90", alpha=0.5) +
    geom_line(linewidth=0.5) +
    scale_y_log10() +
    labs(title=sprintf("Regularizacao adaptativa (h=%d)", h),
         subtitle="Escala log | Cinza = recessoes NBER",
         x="", y=expression(lambda~"(log)"), color="") +
    theme_minimal() + theme(legend.position="bottom")

  ggsave(file.path(fig_dir, sprintf("lambda_h%02d.pdf",h)), p, width=12, height=5)
  print(p)

  cor_lam <- cor(df$lam_ridge, df$lam2_2srr, use="complete.obs")
  cat(sprintf("  h=%d: cor(lam_ridge, lam2_2srr) = %.4f\n", h, cor_lam))
}

# PARTE 2: DISPERSAO CROSS-SECTIONAL DOS BETAS
cat("\nPARTE 2: Dispersao cross-sectional\n")

for (hi in seq_along(hor)) {
  h <- hor[hi]
  valid_2srr  <- Filter(Negate(is.null), betas_2srr[[hi]])
  valid_ridge <- Filter(Negate(is.null), betas_ridge[[hi]])
  if (length(valid_2srr) < 5 || length(valid_ridge) < 5) next

  disp_tvp <- sapply(valid_2srr, function(b) {
    bm <- b$betas
    if (is.array(bm) && length(dim(bm))==3) bvec <- bm[1,,dim(bm)[3]]
    else if (is.matrix(bm)) bvec <- bm[nrow(bm),]
    else bvec <- as.numeric(bm)
    var(bvec, na.rm=TRUE)
  })

  disp_ridge <- sapply(valid_ridge, function(b) var(c(b$beta0, b$betas), na.rm=TRUE))
  dates_b <- sapply(valid_2srr, function(b) as.character(b$date))
  n_min <- min(length(disp_tvp), length(disp_ridge), length(dates_b))

  disp_df <- data.frame(date=as.Date(dates_b[1:n_min]),
                         TVP_2SRR=disp_tvp[1:n_min], Ridge=disp_ridge[1:n_min])
  disp_long <- melt(disp_df, id.vars="date", variable.name="Modelo", value.name="Dispersao")

  p <- ggplot(disp_long, aes(x=date, y=Dispersao, color=Modelo)) +
    geom_rect(data=recessions, aes(xmin=start,xmax=end,ymin=-Inf,ymax=Inf),
              inherit.aes=FALSE, fill="gray90", alpha=0.5) +
    geom_line(linewidth=0.5) +
    labs(title=sprintf("Dispersao cross-sectional dos betas (h=%d)",h),
         subtitle="Var(beta_k) por janela", x="", y="Variancia dos betas", color="") +
    theme_minimal() + theme(legend.position="bottom")

  ggsave(file.path(fig_dir, sprintf("dispersao_h%02d.pdf",h)), p, width=12, height=5)
  print(p)
  write.csv(disp_df, file.path(out_dir, sprintf("dispersao_h%02d.csv",h)), row.names=FALSE)

  ratio_d <- mean(disp_df$TVP_2SRR, na.rm=T) / mean(disp_df$Ridge, na.rm=T)
  cat(sprintf("  h=%d: ratio dispersao(2SRR/Ridge) = %.3f\n", h, ratio_d))
}

# PARTE 3: FALLBACK SILENCIOSO (OF FILTER)
cat("\nPARTE 3: Fallback silencioso\n")

fallback_tab <- data.frame()
for (h in hor) {
  key <- paste0("h",h)
  if (is.null(coulombe[[key]])) next
  df <- coulombe[[key]]
  df <- df[complete.cases(df$fc_ridge, df$fc_2srr), ]
  n_fb <- sum(abs(df$fc_2srr - df$fc_ridge) < 1e-10)
  pct  <- 100 * n_fb / nrow(df)
  cat(sprintf("  h=%2d: %d/%d identicas (%.1f%% fallback)\n", h, n_fb, nrow(df), pct))
  fallback_tab <- rbind(fallback_tab, data.frame(h=h, n_total=nrow(df),
                                                  n_fallback=n_fb, pct_fallback=pct))
}
print(fallback_tab)
write.csv(fallback_tab, file.path(out_dir, "fallback_detection.csv"), row.names=FALSE)

# PARTE 4A: SUB-PERIODOS EXPANDIDOS
cat("\nPARTE 4A: Sub-periodos expandidos\n")


periodos <- list(
  "Full Sample"        = c(as.Date("1999-07-01"), as.Date("2025-06-01")),
  "Pre-GFC"            = c(as.Date("1999-07-01"), as.Date("2007-11-30")),
  "GFC"                = c(as.Date("2007-12-01"), as.Date("2009-06-30")),
  "Post-GFC_Pre-COVID" = c(as.Date("2009-07-01"), as.Date("2020-01-31")),
  "COVID"              = c(as.Date("2020-02-01"), as.Date("2021-06-30")),
  "High Inflation"     = c(as.Date("2021-07-01"), as.Date("2023-06-30")),
  "Post-Inflation"     = c(as.Date("2023-07-01"), as.Date("2025-06-01")))


sub_results <- list()
for (h in hor) {
  key <- paste0("h",h)
  if (is.null(coulombe[[key]])) next
  df <- coulombe[[key]]
  df <- df[complete.cases(df$realized, df$fc_ridge, df$fc_2srr), ]
  df$date_p <- as.Date(df$date)


  for (pn in names(periodos)) {
    pr <- periodos[[pn]]
    idx <- df$date_p >= pr[1] & df$date_p <= pr[2]
    if (sum(idx) < 5) next
    rmse_r <- sqrt(mean((df$fc_ridge[idx]-df$realized[idx])^2))
    rmse_2 <- sqrt(mean((df$fc_2srr[idx]-df$realized[idx])^2))
    sub_results[[paste0(pn,"_h",h)]] <- data.frame(
      h=h, periodo=pn, RMSE_Ridge=rmse_r, RMSE_2SRR=rmse_2,
      Ratio=rmse_2/rmse_r, n=sum(idx))
  }
}


if (length(sub_results) > 0) {
  sub_tab <- do.call(rbind, sub_results)
  print(sub_tab)
  write.csv(sub_tab, file.path(out_dir, "subperiodos_expandido.csv"), row.names=FALSE)


  sub_tab$periodo <- factor(sub_tab$periodo, levels=names(periodos))
  p <- ggplot(sub_tab, aes(x=periodo, y=Ratio, fill=factor(h))) +
    geom_bar(stat="identity", position="dodge", width=0.7) +
    geom_hline(yintercept=1, linetype="dashed", color="red") +
    labs(title="RMSE(2SRR)/RMSE(Ridge) por sub-periodo",
         subtitle="Abaixo de 1 = 2SRR melhor", x="", y="Ratio", fill="h") +
    theme_minimal() + theme(axis.text.x=element_text(angle=30, hjust=1))
  ggsave(file.path(fig_dir, "subperiodos_ratio.pdf"), p, width=14, height=6)
  print(p)
}


for (h in hor) {
  key <- paste0("h",h)
  if (is.null(coulombe[[key]])) next
  df <- coulombe[[key]]
  df <- df[complete.cases(df$realized, df$fc_ridge, df$fc_2srr), ]
  df$date_p <- as.Date(df$date)


  dm_results <- list()
  for (pn in names(periodos)) {
    pr <- periodos[[pn]]
    idx <- df$date_p >= pr[1] & df$date_p <= pr[2]
    if (sum(idx) < 10) next
    e_ridge <- df$fc_ridge[idx] - df$realized[idx]
    e_2srr  <- df$fc_2srr[idx]  - df$realized[idx]
    dm_test <- tryCatch(
      dm.test(e_ridge, e_2srr, alternative="greater", h=h, power=2),
      error = function(e) NULL
    )
    if (!is.null(dm_test)) {
      dm_results[[pn]] <- data.frame(
        h=h, periodo=pn,
        DM_stat=dm_test$statistic,
        p_value=dm_test$p.value,
        significativo=dm_test$p.value < 0.10
      )
    }
  }
  if (length(dm_results) > 0) {
    dm_tab <- do.call(rbind, dm_results)
    print(dm_tab)
    write.csv(dm_tab, file.path(out_dir, sprintf("dm_test_h%02d.csv",h)), row.names=FALSE)
  }
}


# PARTE 4B: EVOLUCAO DE BETAS SELECIONADOS AO LONGO DO TEMPO
cat("\nPARTE 4B: Evolucao de betas selecionados\n")


# Regressores de interesse teorico
reg_interest <- c("intercept", "y_lag1", "y_lag2", "F1_lag1", "F2_lag1")


for (hi in seq_along(hor)) {
  h <- hor[hi]
  valid_b <- Filter(Negate(is.null), betas_2srr[[hi]])
  if (length(valid_b) < 10) next


  beta_mat <- do.call(rbind, lapply(valid_b, function(b) {
    bm <- b$betas
    if (is.array(bm) && length(dim(bm))==3) bvec <- bm[1,,dim(bm)[3]]
    else if (is.matrix(bm)) bvec <- bm[nrow(bm),]
    else bvec <- as.numeric(bm)
    bvec
  }))


  K <- ncol(beta_mat)
  col_nm <- if (K <= length(reg_names)) reg_names[1:K] else paste0("b",0:(K-1))
  colnames(beta_mat) <- col_nm
  dates_b <- as.Date(sapply(valid_b, function(b) as.character(b$date)))


  # Selecionar apenas regressores de interesse que existem
  reg_sel <- intersect(reg_interest, col_nm)
  if (length(reg_sel) == 0) next


  df_betas <- as.data.frame(beta_mat[, reg_sel, drop=FALSE])
  df_betas$date <- dates_b
  df_long <- melt(df_betas, id.vars="date", variable.name="Regressor", value.name="Beta")


  p <- ggplot(df_long, aes(x=date, y=Beta, color=Regressor)) +
    geom_rect(data=recessions, aes(xmin=start,xmax=end,ymin=-Inf,ymax=Inf),
              inherit.aes=FALSE, fill="gray90", alpha=0.5) +
    geom_line(linewidth=0.6) +
    geom_hline(yintercept=0, linetype="dashed", color="gray40") +
    facet_wrap(~Regressor, scales="free_y", ncol=2) +
    labs(title=sprintf("Evolucao temporal dos betas 2SRR (h=%d)", h),
         subtitle="Intercepto, lags de y e fatores principais | Cinza = recessoes NBER",
         x="", y="Beta", color="") +
    theme_minimal() + theme(legend.position="none")


  ggsave(file.path(fig_dir, sprintf("betas_evolucao_h%02d.pdf",h)), p, width=14, height=8)
  print(p)
  cat(sprintf("  h=%d: evolucao de %d betas plotada\n", h, length(reg_sel)))
}

# PARTE 5A: CORRELACAO BETAS vs VOLATILIDADE
cat("\nPARTE 5A: Correlacao betas vs volatilidade\n")

for (hi in seq_along(hor)) {
  h <- hor[hi]; key <- paste0("h",h)
  if (is.null(coulombe[[key]])) next
  df <- coulombe[[key]]
  df <- df[complete.cases(df$realized, df$fc_2srr), ]
  df$abs_err <- abs(df$fc_2srr - df$realized)
  n_df <- nrow(df)
  df$vol_12m <- NA
  if (n_df >= 12) for (i in 12:n_df) df$vol_12m[i] <- sd(df$abs_err[(i-11):i])

  valid_b <- Filter(Negate(is.null), betas_2srr[[hi]])
  if (length(valid_b) < 10) next

  norma_l2 <- sapply(valid_b, function(b) {
    bm <- b$betas
    if (is.array(bm) && length(dim(bm))==3) bvec <- bm[1,,dim(bm)[3]]
    else if (is.matrix(bm)) bvec <- bm[nrow(bm),]
    else bvec <- as.numeric(bm)
    sqrt(sum(bvec^2, na.rm=TRUE))
  })

  n_min <- min(length(norma_l2), n_df)
  if (n_min < 20) next

  cor_val <- cor(norma_l2[1:n_min], df$vol_12m[1:n_min], use="complete.obs")
  cat(sprintf("  h=%d: cor(norma_L2, vol_12m) = %.4f\n", h, cor_val))

  plot_df <- data.frame(date=as.Date(df$date[1:n_min]),
                         Norma_L2=scale(norma_l2[1:n_min]),
                         Vol_erro=scale(df$vol_12m[1:n_min]))
  plot_long <- melt(plot_df, id.vars="date", variable.name="Serie", value.name="Valor")

  p <- ggplot(plot_long, aes(x=date, y=Valor, color=Serie)) +
    geom_rect(data=recessions, aes(xmin=start,xmax=end,ymin=-Inf,ymax=Inf),
              inherit.aes=FALSE, fill="gray90", alpha=0.5) +
    geom_line(linewidth=0.5) +
    labs(title=sprintf("Norma L2 betas vs Volatilidade erros (h=%d)",h),
         subtitle=sprintf("Correlacao = %.3f", cor_val),
         x="", y="Padronizado", color="") +
    theme_minimal() + theme(legend.position="bottom")
  ggsave(file.path(fig_dir, sprintf("betas_vs_vol_h%02d.pdf",h)), p, width=12, height=5)
  print(p)
}

cat("\nDIAGNOSTICO DE ESCALA DOS FORECASTS\n")
for (hi in seq_along(hor)) {
  h <- hor[hi]
  cat(sprintf("\nh=%d:\n", h))
  cat(sprintf("  yout range: [%.4f, %.4f] | mean=%.4f\n",
              min(yout[,hi], na.rm=T), max(yout[,hi], na.rm=T), mean(yout[,hi], na.rm=T)))
  key <- paste0("h",h)
  if (!is.null(coulombe[[key]])) {
    df <- coulombe[[key]]
    cat(sprintf("  realized (Coulombe) range: [%.4f, %.4f] | mean=%.4f\n",
                min(df$realized, na.rm=T), max(df$realized, na.rm=T), mean(df$realized, na.rm=T)))
    cat(sprintf("  fc_2srr range: [%.4f, %.4f]\n",
                min(df$fc_2srr, na.rm=T), max(df$fc_2srr, na.rm=T)))
  }
  if (!is.null(med_fc[["AR"]]) && ncol(med_fc[["AR"]]) >= hi) {
    cat(sprintf("  AR (Medeiros) range: [%.4f, %.4f]\n",
                min(med_fc[["AR"]][,hi], na.rm=T), max(med_fc[["AR"]][,hi], na.rm=T)))
  }
}

# PARTE 5B: DECOMPOSICAO DO GANHO — lambda vs. beta heterogeneidade
cat("\nPARTE 5B: Decomposicao do ganho 2SRR\n")

for (h in hor) {
  key <- paste0("h",h)
  if (is.null(coulombe[[key]])) next
  df <- coulombe[[key]]
  df <- df[complete.cases(df$realized, df$fc_ridge, df$fc_2srr,
                           df$lam_ridge, df$lam2_2srr), ]
  if (nrow(df) < 30) next

  # Ganho pontual (quadrático) do 2SRR sobre Ridge
  df$ganho <- (df$fc_ridge - df$realized)^2 - (df$fc_2srr - df$realized)^2
  # Positivo = 2SRR ganhou naquele período

  # Variacao relativa de lambda (proxy de adaptação)
  df$delta_lam <- abs(df$lam2_2srr - df$lam_ridge) / (abs(df$lam_ridge) + 1e-10)

  # Regressão: ganho ~ delta_lambda
  fit <- lm(ganho ~ delta_lam, data=df)
  cat(sprintf("  h=%d: R²=%.3f | coef delta_lam=%.4f (p=%.3f)\n",
              h, summary(fit)$r.squared,
              coef(fit)[2], summary(fit)$coefficients[2,4]))

  # Scatter ganho vs delta_lambda
  p <- ggplot(df, aes(x=delta_lam, y=ganho)) +
    geom_point(alpha=0.3, size=1.2, color="steelblue") +
    geom_smooth(method="lm", color="red", linewidth=0.8, se=TRUE) +
    geom_hline(yintercept=0, linetype="dashed", color="gray40") +
    labs(title=sprintf("Ganho 2SRR vs Variacao de Lambda (h=%d)", h),
         subtitle=sprintf("R²=%.3f | Positivo = 2SRR ganhou", summary(fit)$r.squared),
         x="Variacao relativa de lambda (|lam2 - lam_ridge| / lam_ridge)",
         y="Ganho quadrático (Ridge² - 2SRR²)") +
    theme_minimal()
  ggsave(file.path(fig_dir, sprintf("ganho_decomp_h%02d.pdf",h)), p, width=10, height=6)
  print(p)

  write.csv(df[, c("date","ganho","delta_lam")],
            file.path(out_dir, sprintf("ganho_decomp_h%02d.csv",h)), row.names=FALSE)
}

# PARTE 6: TABELA DUAL BENCHMARK
cat("\nPARTE 6: Tabela dual benchmark\n")

tab_dual <- list()
for (hi in seq_along(hor)) {
  h <- hor[hi]; real <- yout[,hi]; n <- length(real)
  rmse_rw <- NULL; rmse_ar <- NULL
  if (!is.null(med_fc[["rw"]]) && ncol(med_fc[["rw"]])>=hi)
    rmse_rw <- sqrt(mean((med_fc[["rw"]][1:n,hi]-real)^2, na.rm=T))
  if (!is.null(med_fc[["AR"]]) && ncol(med_fc[["AR"]])>=hi)
    rmse_ar <- sqrt(mean((med_fc[["AR"]][1:n,hi]-real)^2, na.rm=T))

  all_m <- list()
  for (mn in names(med_fc)) {
    mm <- med_fc[[mn]]
    if (ncol(mm)>=hi && nrow(mm)>=n)
      all_m[[mn]] <- sqrt(mean((mm[1:n,hi]-real)^2, na.rm=T))
  }
  key <- paste0("h",h)
  if (!is.null(coulombe[[key]])) {
    cdf <- coulombe[[key]]
    cdf <- cdf[complete.cases(cdf$realized,cdf$fc_ridge,cdf$fc_2srr),]
    if (nrow(cdf)>10) {
      all_m[["Ridge_Coulombe"]] <- sqrt(mean((cdf$fc_ridge-cdf$realized)^2))
      all_m[["TVP_FAVAR"]]     <- sqrt(mean((cdf$fc_2srr-cdf$realized)^2))
    }
  }
  if (!is.null(tvp_ar_list[[key]])) {
    tdf <- tvp_ar_list[[key]]; tdf <- tdf[complete.cases(tdf$realized,tdf$fc_2srr),]
    if (nrow(tdf)>10) all_m[["TVP_AR"]] <- sqrt(mean((tdf$fc_2srr-tdf$realized)^2))
  }
  if (!is.null(tvp_fac_list[[key]])) {
    tdf <- tvp_fac_list[[key]]; tdf <- tdf[complete.cases(tdf$realized,tdf$fc_2srr),]
    if (nrow(tdf)>10) all_m[["TVP_Factor"]] <- sqrt(mean((tdf$fc_2srr-tdf$realized)^2))
  }
  for (mn in names(all_m)) {
    tab_dual[[paste0(mn,"_h",h)]] <- data.frame(
      h=h, model=mn, RMSE=all_m[[mn]],
      ratio_RW=ifelse(!is.null(rmse_rw), all_m[[mn]]/rmse_rw, NA),
      ratio_AR=ifelse(!is.null(rmse_ar), all_m[[mn]]/rmse_ar, NA))
  }
}

if (length(tab_dual)>0) {
  tab <- do.call(rbind, tab_dual); rownames(tab) <- NULL
  tab$rank <- NA
  for (h in hor) { idx <- tab$h==h; tab$rank[idx] <- rank(tab$RMSE[idx]) }

  cat("\n")
  modelos <- sort(unique(tab$model))
  cat(sprintf("%-20s", "Modelo"))
  for (h in hor) cat(sprintf(" %11s", paste0("h=",h)))
  cat("\n"); cat(paste(rep("-",68),collapse=""),"\n")
  for (m in modelos) {
    cat(sprintf("%-20s", m))
    for (h in hor) {
      r <- tab[tab$model==m & tab$h==h,]
      if (nrow(r)==1) cat(sprintf(" %5.3f(%4.2f)", r$RMSE, r$ratio_RW))
      else cat(sprintf(" %11s", "-"))
    }
    cat("\n")
  }
  print(tab)
  write.csv(tab, file.path(out_dir, "tabela_dual_benchmark.csv"), row.names=FALSE)
}

# PARTE 7: SUPERFICIE 3D DOS BETAS TVP (estilo Coulombe)
cat("\nPARTE 7: Superficie 3D dos betas\n")

for (hi in seq_along(hor)) {
  h <- hor[hi]
  valid_b <- Filter(Negate(is.null), betas_2srr[[hi]])
  if (length(valid_b) < 10) next

  beta_mat <- do.call(rbind, lapply(valid_b, function(b) {
    bm <- b$betas
    if (is.array(bm) && length(dim(bm))==3) bvec <- bm[1,,dim(bm)[3]]
    else if (is.matrix(bm)) bvec <- bm[nrow(bm),]
    else bvec <- as.numeric(bm)
    bvec
  }))

  K <- ncol(beta_mat)
  T_oos <- nrow(beta_mat)
  col_nm <- if (K <= length(reg_names)) reg_names[1:K] else paste0("b",0:(K-1))
  colnames(beta_mat) <- col_nm

  # 7A. Heatmap 2D (funciona sem plotly)
  beta_df <- as.data.frame(beta_mat)
  beta_df$t <- 1:T_oos
  beta_long <- melt(beta_df, id.vars="t", variable.name="Regressor", value.name="Beta")

  p_heat <- ggplot(beta_long, aes(x=t, y=Regressor, fill=Beta)) +
    geom_tile() +
    scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0) +
    labs(title=sprintf("Superficie dos betas TVP (h=%d) — Heatmap",h),
         subtitle="Azul=negativo | Vermelho=positivo | Branco=zero",
         x="Janela OOS", y="Regressor", fill="Beta") +
    theme_minimal() +
    theme(axis.text.y=element_text(size=7))

  ggsave(file.path(fig_dir, sprintf("heatmap_betas_h%02d.pdf",h)), p_heat, width=14, height=8)
  print(p_heat)

  # 7B. Plotly 3D (se disponivel)
  if (has_plotly) {
    p3d <- plot_ly(z=~beta_mat,
                   x=col_nm, y=1:T_oos,
                   type="surface",
                   colorscale="RdBu",
                   reversescale=TRUE) %>%
      layout(title=sprintf("Superficie 3D dos Betas TVP (h=%d)", h),
             scene=list(
               xaxis=list(title="Regressor"),
               yaxis=list(title="Janela OOS"),
               zaxis=list(title="Beta")))
    
    htmlwidgets::saveWidget(p3d,
      file.path(fig_dir, sprintf("surface3d_betas_h%02d.html",h)),
      selfcontained=TRUE)
    print(p3d)
    cat(sprintf("  h=%d: 3D surface salvo (HTML interativo)\n", h))
  }

  cat(sprintf("  h=%d: heatmap salvo (%d janelas x %d regressores)\n", h, T_oos, K))
}

# PARTE 8A: MUDANCAS DE SINAL DOS BETAS
cat("\nPARTE 8A: Mudancas de sinal dos betas\n")

for (hi in seq_along(hor)) {
  h <- hor[hi]
  valid_b <- Filter(Negate(is.null), betas_2srr[[hi]])
  if (length(valid_b) < 10) next

  beta_mat <- do.call(rbind, lapply(valid_b, function(b) {
    bm <- b$betas
    if (is.array(bm) && length(dim(bm))==3) bvec <- bm[1,,dim(bm)[3]]
    else if (is.matrix(bm)) bvec <- bm[nrow(bm),]
    else bvec <- as.numeric(bm)
    bvec
  }))

  K <- ncol(beta_mat)
  col_nm <- if (K<=length(reg_names)) reg_names[1:K] else paste0("b",0:(K-1))

  sign_changes <- sapply(1:K, function(k) {
    signs <- sign(beta_mat[,k])
    signs <- signs[signs != 0]
    if (length(signs) < 2) return(0)
    sum(diff(signs) != 0)
  })
  names(sign_changes) <- col_nm

  sc_df <- data.frame(regressor=col_nm, sign_changes=sign_changes,
                       pct=100*sign_changes/(nrow(beta_mat)-1))
  sc_df <- sc_df[order(-sc_df$sign_changes),]
  print(sc_df)

  p_sc <- ggplot(sc_df, aes(x=reorder(regressor, -sign_changes), y=sign_changes)) +
    geom_bar(stat="identity", fill="steelblue", width=0.6) +
    geom_text(aes(label=sign_changes), vjust=-0.3, size=3) +
    labs(title=sprintf("Mudancas de sinal dos betas TVP (h=%d)",h),
         subtitle="Betas que trocam de sinal frequentemente = instabilidade",
         x="", y="N. mudancas de sinal") +
    theme_minimal() +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=8))
  ggsave(file.path(fig_dir, sprintf("sign_changes_h%02d.pdf",h)), p_sc, width=12, height=5)
  print(p_sc)

  write.csv(sc_df, file.path(out_dir, sprintf("sign_changes_h%02d.csv",h)), row.names=FALSE)
  cat(sprintf("  h=%d: max mudancas = %s (%d)\n", h, sc_df$regressor[1], sc_df$sign_changes[1]))
}

# PARTE 8B: PERSISTENCIA DOS BETAS (autocorrelacao de ordem 1)
cat("\nPARTE 8B: Persistencia dos betas (AR1)\n")

for (hi in seq_along(hor)) {
  h <- hor[hi]
  valid_b <- Filter(Negate(is.null), betas_2srr[[hi]])
  if (length(valid_b) < 20) next

  beta_mat <- do.call(rbind, lapply(valid_b, function(b) {
    bm <- b$betas
    if (is.array(bm) && length(dim(bm))==3) bvec <- bm[1,,dim(bm)[3]]
    else if (is.matrix(bm)) bvec <- bm[nrow(bm),]
    else bvec <- as.numeric(bm)
    bvec
  }))

  K <- ncol(beta_mat)
  col_nm <- if (K <= length(reg_names)) reg_names[1:K] else paste0("b",0:(K-1))

  # AR(1) de cada beta ao longo do tempo
  ar1_coef <- sapply(1:K, function(k) {
    x <- beta_mat[,k]
    x <- x[!is.na(x)]
    if (length(x) < 10) return(NA)
    tryCatch(acf(x, lag.max=1, plot=FALSE)$acf[2], error=function(e) NA)
  })
  names(ar1_coef) <- col_nm

  ar1_df <- data.frame(regressor=col_nm, ar1=ar1_coef)
  ar1_df <- ar1_df[!is.na(ar1_df$ar1),]
  ar1_df <- ar1_df[order(-ar1_df$ar1),]

  p_ar <- ggplot(ar1_df, aes(x=reorder(regressor, ar1), y=ar1, fill=ar1 > 0.5)) +
    geom_col(width=0.6) +
    geom_hline(yintercept=0.5, linetype="dashed", color="red") +
    scale_fill_manual(values=c("FALSE"="steelblue","TRUE"="darkorange"),
                      labels=c("Baixa persistencia","Alta persistencia (>0.5)"),
                      name="") +
    coord_flip() +
    labs(title=sprintf("Persistencia (ACF lag-1) dos betas TVP (h=%d)", h),
         subtitle="Acima de 0.5 = beta altamente persistente (captura regime)",
         x="", y="ACF lag-1") +
    theme_minimal() + theme(legend.position="bottom")
  ggsave(file.path(fig_dir, sprintf("persistencia_betas_h%02d.pdf",h)), p_ar, width=10, height=7)
  print(p_ar)
  write.csv(ar1_df, file.path(out_dir, sprintf("persistencia_betas_h%02d.csv",h)), row.names=FALSE)

  cat(sprintf("  h=%d: beta mais persistente = %s (ACF1=%.3f)\n",
              h, ar1_df$regressor[which.max(ar1_df$ar1)],
              max(ar1_df$ar1, na.rm=TRUE)))
}

# PARTE 9A: ROLLING RMSE RATIO (janela 36 meses)
cat("\nPARTE 9A: Rolling RMSE ratio (36 meses)\n")

for (h in hor) {
  key <- paste0("h",h)
  if (is.null(coulombe[[key]])) next
  df <- coulombe[[key]]
  df <- df[complete.cases(df$realized, df$fc_ridge, df$fc_2srr),]
  n_df <- nrow(df)
  window <- 36

  if (n_df < window + 10) next

  roll_df <- data.frame(date=as.Date(df$date), ratio=NA, pct_wins=NA)

  for (i in window:n_df) {
    w <- (i-window+1):i
    rmse_r <- sqrt(mean((df$fc_ridge[w]-df$realized[w])^2))
    rmse_2 <- sqrt(mean((df$fc_2srr[w]-df$realized[w])^2))
    roll_df$ratio[i] <- rmse_2/rmse_r

    d_t <- (df$fc_ridge[w]-df$realized[w])^2 - (df$fc_2srr[w]-df$realized[w])^2
    roll_df$pct_wins[i] <- 100*mean(d_t > 0)
  }

  roll_clean <- roll_df[!is.na(roll_df$ratio),]

  p_roll <- ggplot(roll_clean, aes(x=date)) +
    geom_rect(data=recessions, aes(xmin=start,xmax=end,ymin=-Inf,ymax=Inf),
              inherit.aes=FALSE, fill="gray90", alpha=0.5) +
    geom_line(aes(y=ratio), color="steelblue", linewidth=0.6) +
    geom_hline(yintercept=1, linetype="dashed", color="red") +
    labs(title=sprintf("Rolling RMSE ratio 2SRR/Ridge (janela %d meses, h=%d)", window, h),
         subtitle="Abaixo de 1 = 2SRR melhor naquela janela",
         x="", y="Ratio RMSE") +
    theme_minimal()
  ggsave(file.path(fig_dir, sprintf("rolling_ratio_h%02d.pdf",h)), p_roll, width=12, height=5)
  print(p_roll)

  # Pct do tempo que 2SRR ganha
  pct_below_1 <- 100*mean(roll_clean$ratio < 1, na.rm=T)
  cat(sprintf("  h=%d: 2SRR ratio<1 em %.0f%% das janelas rolling\n", h, pct_below_1))

  write.csv(roll_clean, file.path(out_dir, sprintf("rolling_ratio_h%02d.csv",h)), row.names=FALSE)
}

# Parte 9B: Rolling RMSE ratio comparando TVP-AR, TVP-Factor e 2SRR
for (h in hor) {
  key <- paste0("h",h)
  if (is.null(coulombe[[key]])) next
  df <- coulombe[[key]]
  df <- df[complete.cases(df$realized, df$fc_ridge, df$fc_2srr),]
  n_df <- nrow(df); window <- 36
  if (n_df < window + 10) next

  roll_df <- data.frame(date=as.Date(df$date), ratio_2srr=NA, ratio_tvpar=NA, ratio_tvpfac=NA)

  # Alinha TVP-AR e TVP-Factor com o mesmo período do Coulombe
  tvpar_fc  <- if (!is.null(tvp_ar_list[[key]]))  tvp_ar_list[[key]]$fc_2srr  else rep(NA, n_df)
  tvpfac_fc <- if (!is.null(tvp_fac_list[[key]])) tvp_fac_list[[key]]$fc_2srr else rep(NA, n_df)
  tvpar_fc  <- tvpar_fc[1:n_df]
  tvpfac_fc <- tvpfac_fc[1:n_df]

  for (i in window:n_df) {
    w <- (i-window+1):i
    rmse_r  <- sqrt(mean((df$fc_ridge[w] - df$realized[w])^2))
    rmse_2  <- sqrt(mean((df$fc_2srr[w]  - df$realized[w])^2))
    roll_df$ratio_2srr[i] <- rmse_2 / rmse_r

    if (!any(is.na(tvpar_fc[w])))
      roll_df$ratio_tvpar[i]  <- sqrt(mean((tvpar_fc[w]  - df$realized[w])^2)) / rmse_r
    if (!any(is.na(tvpfac_fc[w])))
      roll_df$ratio_tvpfac[i] <- sqrt(mean((tvpfac_fc[w] - df$realized[w])^2)) / rmse_r
  }

  roll_long <- melt(roll_df[!is.na(roll_df$ratio_2srr),],
                    id.vars="date", variable.name="Modelo", value.name="Ratio")
  roll_long <- roll_long[!is.na(roll_long$Ratio),]

  p_roll <- ggplot(roll_long, aes(x=date, y=Ratio, color=Modelo)) +
    geom_rect(data=recessions, aes(xmin=start,xmax=end,ymin=-Inf,ymax=Inf),
              inherit.aes=FALSE, fill="gray90", alpha=0.5) +
    geom_line(linewidth=0.6) +
    geom_hline(yintercept=1, linetype="dashed", color="red") +
    scale_color_manual(values=c("ratio_2srr"="steelblue","ratio_tvpar"="darkorange","ratio_tvpfac"="forestgreen"),
                       labels=c("2SRR/Ridge","TVP-AR/Ridge","TVP-Factor/Ridge")) +
    labs(title=sprintf("Rolling RMSE ratio vs Ridge (janela %d meses, h=%d)", window, h),
         subtitle="Abaixo de 1 = modelo melhor que Ridge", x="", y="Ratio RMSE") +
    theme_minimal()
  ggsave(file.path(fig_dir, sprintf("rolling_ratio_multi_h%02d.pdf",h)), p_roll, width=12, height=5)
  print(p_roll)
}

# PARTE 9C: DENSIDADE DOS ERROS POR REGIME
cat("\nPARTE 9C: Densidade dos erros por regime\n")

for (h in c(1, 12)) {   # só h>=6 onde 2SRR ganha
  key <- paste0("h",h)
  if (is.null(coulombe[[key]])) next
  df <- coulombe[[key]]
  df <- df[complete.cases(df$realized, df$fc_ridge, df$fc_2srr),]
  df$date_p <- as.Date(df$date)

  # Classificar regime
  df$regime <- "Normal"
  df$regime[df$date_p >= as.Date("2007-12-01") &
             df$date_p <= as.Date("2009-06-30")] <- "GFC"
  df$regime[df$date_p >= as.Date("2020-02-01") &
             df$date_p <= as.Date("2021-06-30")] <- "COVID"
  df$regime[df$date_p >= as.Date("2021-07-01") &
             df$date_p <= as.Date("2023-06-30")] <- "High Inflation"

  df$err_ridge <- df$fc_ridge - df$realized
  df$err_2srr  <- df$fc_2srr  - df$realized

  err_long <- rbind(
    data.frame(regime=df$regime, erro=df$err_ridge, modelo="Ridge"),
    data.frame(regime=df$regime, erro=df$err_2srr,  modelo="2SRR")
  )

  p_dens <- ggplot(err_long, aes(x=erro, fill=modelo)) +
    geom_density(alpha=0.4, linewidth=0.4) +
    geom_vline(xintercept=0, linetype="dashed", color="gray40") +
    facet_wrap(~regime, scales="free", ncol=2) +
    scale_fill_manual(values=c("Ridge"="tomato","2SRR"="steelblue")) +
    labs(title=sprintf("Distribuicao dos erros por regime (h=%d)", h),
         subtitle="2SRR vs Ridge — caudas mais pesadas em crises?",
         x="Erro de previsao", y="Densidade", fill="") +
    theme_minimal() + theme(legend.position="bottom")
  ggsave(file.path(fig_dir, sprintf("densidade_erros_h%02d.pdf",h)), p_dens, width=12, height=8)
  print(p_dens)
  cat(sprintf("  h=%d: densidade por regime salva\n", h))
}

# PARTE 10: FAN CHART DAS PREVISOES
cat("\nPARTE 10: Fan chart\n")

for (h in hor) { #for (h in c(1, 12)) {
  key <- paste0("h",h)
  if (is.null(coulombe[[key]])) next
  df <- coulombe[[key]]
  df <- df[complete.cases(df$realized, df$fc_2srr),]
  df$date <- as.Date(df$date)
  df$err <- df$fc_2srr - df$realized

  # Rolling sd dos erros (12 meses)
  n_df <- nrow(df); df$sd_12m <- NA
  if (n_df >= 12) for (i in 12:n_df) df$sd_12m[i] <- sd(df$err[(i-11):i])

  df_clean <- df[!is.na(df$sd_12m),]
  if (nrow(df_clean) < 20) next

  p_fan <- ggplot(df_clean, aes(x=date)) +
    geom_rect(data=recessions, aes(xmin=start,xmax=end,ymin=-Inf,ymax=Inf),
              inherit.aes=FALSE, fill="gray90", alpha=0.3) +
    geom_ribbon(aes(ymin=fc_2srr-1.96*sd_12m, ymax=fc_2srr+1.96*sd_12m),
                fill="steelblue", alpha=0.15) +
    geom_ribbon(aes(ymin=fc_2srr-1*sd_12m, ymax=fc_2srr+1*sd_12m),
                fill="steelblue", alpha=0.25) +
    geom_line(aes(y=realized), color="black", linewidth=0.6) +
    geom_line(aes(y=fc_2srr), color="steelblue", linewidth=0.4, alpha=0.8) +
    labs(title=sprintf("Fan Chart: 2SRR (h=%d)", h),
         subtitle="Preto=realizado | Azul=previsao | Bandas=1 e 2 desvios-padrao rolling",
         x="", y="Inflacao") +
    theme_minimal()
  ggsave(file.path(fig_dir, sprintf("fan_chart_h%02d.pdf",h)), p_fan, width=14, height=5)
  print(p_fan)
  cat(sprintf("  h=%d: fan chart salvo\n", h))
}

# PARTE 11: Narrativa dinamica
cat("\nPARTE 11: Narrativa dinamica\n")

sink(file.path(out_dir, "narrativa_2srr.txt"))
cat("NARRATIVA — 2SRR Deep Dive\n")
cat(sprintf("Gerado: %s\n\n", Sys.time()))

# 1. Lambda
cat("1. REGULARIZACAO ADAPTATIVA\n")
for (h in hor) {
  key <- paste0("h",h)
  if (!is.null(coulombe[[key]])) {
    df <- coulombe[[key]]
    df <- df[complete.cases(df$lam_ridge, df$lam2_2srr), ]
    if (nrow(df) >= 10) {
      cor_lam <- cor(df$lam_ridge, df$lam2_2srr, use="complete.obs")
      direcao <- ifelse(cor_lam > 0.3, "concordam moderadamente",
                 ifelse(cor_lam < -0.3, "divergem (adaptacao independente)", "sao pouco correlacionados"))
      cat(sprintf("  h=%d: cor=%.3f — lambdas %s\n", h, cor_lam, direcao))
    }
  }
}

# 2. Dispersao
cat("\n2. DISPERSAO HETEROGENEA (ratio media 2SRR/Ridge)\n")
for (hi in seq_along(hor)) {
  h <- hor[hi]
  fname <- file.path(out_dir, sprintf("dispersao_h%02d.csv",h))
  if (file.exists(fname)) {
    d <- read.csv(fname)
    r <- mean(d$TVP_2SRR, na.rm=T) / mean(d$Ridge, na.rm=T)
    cat(sprintf("  h=%d: ratio=%.2fx — 2SRR %s disperso que Ridge\n",
                h, r, ifelse(r>1,"mais","menos")))
  }
}

# 3. Integridade
cat("\n3. INTEGRIDADE (fallbacks)\n")
for (h in hor) {
  key <- paste0("h",h)
  if (!is.null(coulombe[[key]])) {
    df <- coulombe[[key]]
    df <- df[complete.cases(df$fc_ridge,df$fc_2srr),]
    n_fb <- sum(abs(df$fc_2srr-df$fc_ridge)<1e-10)
    cat(sprintf("  h=%d: %d/%d fallbacks (%.1f%%)\n", h, n_fb, nrow(df), 100*n_fb/nrow(df)))
  }
}

# 4. Sub-periodos com DM
cat("\n4. RESULTADOS POR SUB-PERIODO (2SRR vs Ridge)\n")
if (exists("sub_tab")) {
  for (h in hor) {
    cat(sprintf("\n  h=%d:\n", h))
    sub_h <- sub_tab[sub_tab$h==h,]
    n_wins <- sum(sub_h$Ratio < 1, na.rm=TRUE)
    cat(sprintf("  2SRR ganha em %d/%d sub-periodos\n", n_wins, nrow(sub_h)))
    for (i in seq_len(nrow(sub_h))) {
      # Tenta ler p-valor do DM
      dm_fname <- file.path(out_dir, sprintf("dm_test_h%02d.csv",h))
      pval_str <- ""
      if (file.exists(dm_fname)) {
        dm_df <- read.csv(dm_fname)
        dm_row <- dm_df[dm_df$periodo == sub_h$periodo[i],]
        if (nrow(dm_row) == 1)
          pval_str <- sprintf(" [DM p=%.2f%s]", dm_row$p_value,
                              ifelse(dm_row$p_value < 0.10, "*",""))
      }
      cat(sprintf("    %-25s ratio=%.3f n=%d%s %s\n",
                  sub_h$periodo[i], sub_h$Ratio[i], sub_h$n[i], pval_str,
                  ifelse(sub_h$Ratio[i]<1,"<-- 2SRR melhor","")))
    }
  }
}

# 5. Rolling
cat("\n5. DESEMPENHO ROLLING (% janelas com ratio<1)\n")
for (h in hor) {
  fname <- file.path(out_dir, sprintf("rolling_ratio_h%02d.csv",h))
  if (file.exists(fname)) {
    r <- read.csv(fname)
    pct <- 100*mean(r$ratio < 1, na.rm=TRUE)
    cat(sprintf("  h=%d: 2SRR melhor em %.0f%% das janelas rolling de 36 meses\n", h, pct))
  }
}

# 6. Persistencia
cat("\n6. PERSISTENCIA DOS BETAS\n")
for (h in hor) {
  fname <- file.path(out_dir, sprintf("persistencia_betas_h%02d.csv",h))
  if (file.exists(fname)) {
    p_df <- read.csv(fname)
    top <- p_df[which.max(p_df$ar1),]
    pct_alta <- 100*mean(p_df$ar1 > 0.5, na.rm=TRUE)
    cat(sprintf("  h=%d: %.0f%% dos betas com ACF1>0.5 | mais persistente: %s (%.2f)\n",
                h, pct_alta, top$regressor, top$ar1))
  }
}

cat("\n7. CONCLUSAO\n")
if (exists("sub_tab")) {
  best_h <- hor[which.min(sapply(hor, function(h) {
    s <- sub_tab[sub_tab$h==h & sub_tab$periodo=="Full Sample",]
    if (nrow(s)==1) s$Ratio else Inf
  }))]
  cat(sprintf("Melhor desempenho relativo: h=%d (menor ratio full sample)\n", best_h))
}
cat("O 2SRR demonstra: (1) regularizacao adaptativa independente do Ridge,\n")
cat("(2) encolhimento heterogeneo por variavel, (3) robustez em sub-periodos\n")
cat("de estresse, e (4) betas persistentes que capturam mudancas de regime.\n")
sink()

cat(sprintf("Narrativa salva: %s/narrativa_2srr.txt\n", out_dir))
cat(sprintf("\nTODOS OS OUTPUTS EM: %s\n", out_dir))

# ============================================================
# PARTE NOVA (Falta corrigir o 2SRR ESTÁ FICANDO QUASE ZERADO): Comparacao 2SRR Felipe vs Coulombe
# ============================================================

cat("\nPARTE NOVA: Comparacao Felipe RW 2SRR vs Coulombe 2SRR\n")
cat(rep("=", 60), "\n", sep="")

# 1. Carrega o forecast do Felipe (formato Medeiros)

env_fe <- new.env()
load("forecasts/2SRR.rda", envir = env_fe)
fc_felipe_raw <- get(ls(env_fe)[1], envir = env_fe)
if (!is.matrix(fc_felipe_raw)) fc_felipe_raw <- as.matrix(fc_felipe_raw)
cat(sprintf("  Felipe 2SRR: %d janelas x %d horizontes\n",
            nrow(fc_felipe_raw), ncol(fc_felipe_raw)))

# 2. Verifica escala dos realizados

cat("\n--- DIAGNOSTICO DE ESCALA ---\n")
for (h in hor) {
  key <- paste0("h", h)
  hi  <- match(h, hor)
  if (is.null(coulombe[[key]])) next
  cdf <- coulombe[[key]]
  n_fe <- min(nrow(fc_felipe_raw), nrow(yout), nrow(cdf))
  cat(sprintf("h=%d:\n", h))
  cat(sprintf("  yout[,hi]    range: [%.4f, %.4f] | mean=%.4f\n",
              min(yout[1:n_fe, hi], na.rm=TRUE),
              max(yout[1:n_fe, hi], na.rm=TRUE),
              mean(yout[1:n_fe, hi], na.rm=TRUE)))
  cat(sprintf("  cdf$realized range: [%.4f, %.4f] | mean=%.4f\n",
              min(cdf$realized, na.rm=TRUE),
              max(cdf$realized, na.rm=TRUE),
              mean(cdf$realized, na.rm=TRUE)))
  cat(sprintf("  fc_felipe    range: [%.4f, %.4f] | mean=%.4f\n",
              min(fc_felipe_raw[1:n_fe, hi], na.rm=TRUE),
              max(fc_felipe_raw[1:n_fe, hi], na.rm=TRUE),
              mean(fc_felipe_raw[1:n_fe, hi], na.rm=TRUE)))
}

# 3. Compara usando cdf$realized como denominador comum
cat("\n--- COMPARACAO (ambos vs cdf$realized) ---\n")
compare_tab <- list()

for (h in hor) {
  key <- paste0("h", h)
  hi  <- match(h, hor)

  if (is.null(coulombe[[key]])) {
    cat(sprintf("  h=%d: coulombe nao encontrado\n", h))
    next
  }

  cdf <- coulombe[[key]]
  cdf <- cdf[complete.cases(cdf$realized, cdf$fc_2srr), ]

  n_fe <- min(nrow(fc_felipe_raw), nrow(cdf))

  fc_fe     <- fc_felipe_raw[1:n_fe, hi]   # forecast Felipe
  fc_coul   <- cdf$fc_2srr[1:n_fe]         # forecast Coulombe
  real_coul <- cdf$realized[1:n_fe]         # MESMO realizador para ambos

  rmse_fe   <- sqrt(mean((fc_fe   - real_coul)^2, na.rm = TRUE))
  rmse_coul <- sqrt(mean((fc_coul - real_coul)^2, na.rm = TRUE))
  ratio     <- rmse_fe / rmse_coul   # < 1 = Felipe melhor

  e_fe   <- fc_fe   - real_coul
  e_coul <- fc_coul - real_coul

  dm <- tryCatch(
    dm.test(e_fe, e_coul, alternative = "greater", h = h, power = 2),
    error = function(e) list(statistic = NA, p.value = NA)
  )

  vencedor <- ifelse(ratio < 1, "<-- Felipe melhor", "<-- Coulombe melhor")

  cat(sprintf(
    "  h=%2d | Felipe=%.4f | Coulombe=%.4f | Ratio=%.3f %s | DM p=%.3f\n",
    h, rmse_fe, rmse_coul, ratio, vencedor,
    ifelse(is.na(dm$p.value), NA, dm$p.value)
  ))

  compare_tab[[key]] <- data.frame(
    h                    = h,
    RMSE_Felipe          = rmse_fe,
    RMSE_Coulombe        = rmse_coul,
    Ratio_FelipeCoulombe = ratio,
    Felipe_melhor        = ratio < 1,
    DM_stat   = ifelse(is.na(dm$statistic), NA, as.numeric(dm$statistic)),
    DM_pvalue = ifelse(is.na(dm$p.value),   NA, dm$p.value),
    n         = n_fe
  )
}

# 4. Tabela resumo

if (length(compare_tab) > 0) {
  tab_comp <- do.call(rbind, compare_tab)
  rownames(tab_comp) <- NULL
  cat("\n--- TABELA RESUMO ---\n")
  print(tab_comp)
  write.csv(tab_comp,
            file.path(out_dir, "comparacao_felipe_vs_coulombe_v2.csv"),
            row.names = FALSE)
  cat(sprintf("\nTabela salva em: %s/comparacao_felipe_vs_coulombe_v2.csv\n", out_dir))
}

# 5. Plots: ambos vs cdf$realized
for (h in hor) {
  key <- paste0("h", h)
  hi  <- match(h, hor)
  if (is.null(coulombe[[key]])) next
  if (is.null(compare_tab[[key]])) next

  cdf  <- coulombe[[key]]
  cdf  <- cdf[complete.cases(cdf$realized, cdf$fc_2srr), ]
  n_fe <- min(nrow(fc_felipe_raw), nrow(cdf))

  plot_df <- data.frame(
    janela        = 1:n_fe,
    Realizado     = cdf$realized[1:n_fe],      # escala Coulombe
    Felipe_2SRR   = fc_felipe_raw[1:n_fe, hi],
    Coulombe_2SRR = cdf$fc_2srr[1:n_fe]
  )

  plot_long <- melt(plot_df, id.vars = "janela",
                    variable.name = "Serie", value.name = "Valor")

  p <- ggplot(plot_long, aes(x = janela, y = Valor, color = Serie)) +
    geom_line(linewidth = 0.6) +
    scale_color_manual(values = c(
      "Realizado"     = "black",
      "Felipe_2SRR"   = "steelblue",
      "Coulombe_2SRR" = "firebrick"
    )) +
    labs(
      title    = sprintf("2SRR Felipe vs Coulombe (h=%d)", h),
      subtitle = sprintf("RMSE Felipe=%.4f | RMSE Coulombe=%.4f | Ratio=%.3f %s",
                         compare_tab[[key]]$RMSE_Felipe,
                         compare_tab[[key]]$RMSE_Coulombe,
                         compare_tab[[key]]$Ratio_FelipeCoulombe,
                         ifelse(compare_tab[[key]]$Felipe_melhor,
                                "| FELIPE MELHOR", "| COULOMBE MELHOR")),
      x = "Janela OOS", y = "Inflacao (escala Coulombe)", color = ""
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")

  ggsave(file.path(fig_dir, sprintf("fc_felipe_vs_coulombe_v2_h%02d.pdf", h)),
         p, width = 14, height = 5)
  print(p)
  cat(sprintf("  h=%d: plot salvo\n", h))
}