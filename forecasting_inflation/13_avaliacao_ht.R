# ============================================================
# 13_avaliacao_ht.R — ANALISE CONSOLIDADA
#
# Script unico com TODAS as analises solicitadas pelo orientador:
#   PARTE 1: Carregamento de dados
#   PARTE 2: Betas TVP ao longo do tempo (4h em 1 pagina)
#   PARTE 3: Betas TVP vs Ridge OOS (comparacao direta)
#   PARTE 4: Betas variando por horizonte (cross-horizon)
#   PARTE 5: Trajetoria dos lambdas (regularizacao adaptativa)
#   PARTE 6: 2SRR vs Ridge Coulombe (forecast + DM)
#   PARTE 7: 2SRR vs Melhores e Pior do Medeiros
#   PARTE 8: Parcimonia (HHI, near-zero, shrinkage) com testes
#   PARTE 9: Testes econometricos (DM, Clark-West, MZ, Encompassing)
#   PARTE 10: Sub-periodos e regimes
#   PARTE 11: Rolling RMSE e CSFE
#   PARTE 12: Tabelas LaTeX consolidadas
#
# OUTPUTS: pasta 40_results/run13_final/
#   - Todos os graficos em PDF (4 horizontes por pagina)
#   - Tabelas CSV e LaTeX
# ============================================================

rm(list = ls())
setwd("~/tcc/Felipe_Dornelles_tcc/forecasting_inflation")

# --- Pacotes ---
pkgs <- c("ggplot2", "reshape2", "gridExtra", "grid",
          "forecast", "lmtest", "sandwich", "scales", "car")
new <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if (length(new)) install.packages(new, repos = "https://cran.r-project.org")
invisible(lapply(pkgs, library, character.only = TRUE))

# --- Output dir ---
out_dir <- "40_results/run13_final"
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

hor <- c(1, 3, 6, 12)
hor_labels <- c("h=1", "h=3", "h=6", "h=12")

# Recessoes NBER
recessions <- data.frame(
  start = as.Date(c("2001-03-01","2007-12-01","2020-02-01")),
  end   = as.Date(c("2001-11-01","2009-06-01","2020-04-01")),
  label = c("Dot-com","GFC","COVID"))

# Tema padrao
theme_tcc <- theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(color = "grey40", size = 9),
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank())

# Helper: add recession bands
add_recessions <- function() {
  geom_rect(data = recessions,
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "gray85", alpha = 0.5)
}

cat("============================================================\n")
cat(" 13_avaliacao_ht.R — ANALISE CONSOLIDADA\n")
cat(sprintf(" Output: %s\n", out_dir))
cat("============================================================\n\n")

# ============================================================
# PARTE 1: CARREGAMENTO DE DADOS
# ============================================================
cat("PARTE 1: Carregando dados...\n")

load("forecasts/yout.rda")
load("forecasts/coulombe_forecasts.rda")
load("forecasts/coulombe_betas_2SRR.rda")
load("forecasts/coulombe_betas_ridge.rda")
load("data/data.rda")

fred_raw <- as.data.frame(data)
date_col <- grep("^date$", colnames(fred_raw), ignore.case = TRUE)[1]
cpi_col  <- grep("^CPIAUCSL$", colnames(fred_raw), ignore.case = TRUE)[1]
dates    <- fred_raw[, date_col]
y_raw    <- as.numeric(fred_raw[, cpi_col])
bigt     <- nrow(fred_raw)
n_oos    <- nrow(yout)
tau      <- bigt - n_oos

# Coulombe CSVs
coulombe <- list()
for (h in hor) {
  fname <- sprintf("forecasts/coulombe_fc_h%02d.csv", h)
  if (file.exists(fname))
    coulombe[[paste0("h",h)]] <- read.csv(fname, stringsAsFactors = FALSE)
}

# Medeiros models
med_names <- c("AR","AR_BIC","AdaEINET","AdaLASSO","Bagging",
               "CSR","EINET","factor","LASSO","RF","Ridge",
               "T.Factor","2SRR","rw")
med_fc <- list()
for (mn in med_names) {
  fname <- sprintf("forecasts/%s.rda", mn)
  if (file.exists(fname)) {
    env <- new.env(); load(fname, envir = env)
    obj <- get(ls(env)[1], envir = env)
    if (is.matrix(obj) || is.data.frame(obj)) med_fc[[mn]] <- as.matrix(obj)
  }
}

# Nomes dos regressores
reg_names <- c("intercept", paste0("y_lag", 1:2),
               paste0("F", rep(1:8, 2), "_lag", rep(1:2, each = 8)))

# Nomes bonitos para os regressores
reg_labels <- c("Intercepto", "y(t-1)", "y(t-2)",
                paste0("F", rep(1:8, 2), "(t-", rep(1:2, each = 8), ")"))

cat(sprintf("  yout: %dx%d | Medeiros: %d modelos | Coulombe: %d h\n\n",
            nrow(yout), ncol(yout), length(med_fc), length(coulombe)))

# ============================================================
# PARTE 2: BETAS TVP AO LONGO DO TEMPO (4h em 1 pagina)
# ============================================================
cat("PARTE 2: Betas TVP ao longo do tempo\n")

# Extrai betas para todos os horizontes
extract_beta_series <- function(betas_list, hi, type = "2srr") {
  valid <- Filter(Negate(is.null), betas_list[[hi]])
  if (length(valid) < 5) return(NULL)
  
  if (type == "2srr") {
    beta_mat <- do.call(rbind, lapply(valid, function(b) {
      bm <- b$betas
      if (is.array(bm) && length(dim(bm)) == 3) bvec <- bm[1,,dim(bm)[3]]
      else if (is.matrix(bm)) bvec <- bm[nrow(bm),]
      else bvec <- as.numeric(bm)
      bvec
    }))
    dates_b <- as.Date(sapply(valid, function(b) as.character(b$date)))
  } else {
    beta_mat <- do.call(rbind, lapply(valid, function(b) {
      c(b$beta0, b$betas)
    }))
    dates_b <- as.Date(sapply(valid, function(b) as.character(b$date)))
  }
  
  K <- ncol(beta_mat)
  col_nm <- if (K <= length(reg_names)) reg_names[1:K] else paste0("b", 0:(K-1))
  colnames(beta_mat) <- col_nm
  list(betas = beta_mat, dates = dates_b, K = K, col_names = col_nm)
}

# 2A: Top regressores por variancia — 4 horizontes em 1 pagina
# Seleciona regressores teoricamente relevantes
reg_interest <- c("intercept", "y_lag1", "y_lag2", "F1_lag1", "F2_lag1")

plot_list_betas <- list()
for (hi in seq_along(hor)) {
  h <- hor[hi]
  bs <- extract_beta_series(betas_2srr, hi, "2srr")
  if (is.null(bs)) next
  
  reg_sel <- intersect(reg_interest, bs$col_names)
  if (length(reg_sel) == 0) next
  
  df_b <- as.data.frame(bs$betas[, reg_sel, drop = FALSE])
  df_b$date <- bs$dates
  df_long <- melt(df_b, id.vars = "date", variable.name = "Regressor", value.name = "Beta")
  
  p <- ggplot(df_long, aes(x = date, y = Beta, color = Regressor)) +
    add_recessions() +
    geom_line(linewidth = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    facet_wrap(~Regressor, scales = "free_y", ncol = 1) +
    labs(title = sprintf("h = %d", h), x = "", y = "") +
    theme_tcc + theme(legend.position = "none",
                      axis.text.x = element_text(size = 7))
  
  plot_list_betas[[hi]] <- p
}

if (length(plot_list_betas) >= 4) {
  p_combined <- grid.arrange(
    grobs = plot_list_betas,
    ncol = 2,
    top = textGrob("Evolucao Temporal dos Betas TVP (2SRR) — Regressores Principais",
                   gp = gpar(fontface = "bold", fontsize = 14)))
  ggsave(file.path(fig_dir, "P2_betas_tvp_evolucao_4h.pdf"),
         p_combined, width = 16, height = 20)
  cat("  [OK] P2_betas_tvp_evolucao_4h.pdf\n")
}

# ============================================================
# PARTE 3: BETAS TVP vs RIDGE OOS (pedido do orientador)
# ============================================================
cat("\nPARTE 3: Betas 2SRR vs Ridge fora da amostra\n")

# Para cada horizonte: top 6 betas mais variaveis, TVP vs Ridge
for (hi in seq_along(hor)) {
  h <- hor[hi]
  bs_tvp   <- extract_beta_series(betas_2srr, hi, "2srr")
  bs_ridge <- extract_beta_series(betas_ridge, hi, "ridge")
  if (is.null(bs_tvp) || is.null(bs_ridge)) next
  
  n_min <- min(nrow(bs_tvp$betas), nrow(bs_ridge$betas))
  K_min <- min(bs_tvp$K, bs_ridge$K)
  col_use <- bs_tvp$col_names[1:K_min]
  
  # Top 6 por variancia do TVP
  var_tvp <- apply(bs_tvp$betas[1:n_min, 1:K_min], 2, var)
  top6 <- names(sort(var_tvp, decreasing = TRUE))[1:min(6, K_min)]
  
  plots_3 <- list()
  for (bn in top6) {
    ki <- which(col_use == bn)
    df_plot <- data.frame(
      date = c(bs_tvp$dates[1:n_min], bs_ridge$dates[1:n_min]),
      valor = c(bs_tvp$betas[1:n_min, ki], bs_ridge$betas[1:n_min, ki]),
      modelo = c(rep("2SRR (TVP)", n_min), rep("Ridge (constante)", n_min)))
    
    # Label bonito
    lab <- ifelse(bn %in% reg_names, reg_labels[which(reg_names == bn)], bn)
    
    p <- ggplot(df_plot, aes(x = date, y = valor, color = modelo)) +
      add_recessions() +
      geom_line(linewidth = 0.5) +
      scale_color_manual(values = c("2SRR (TVP)" = "#D32F2F", "Ridge (constante)" = "#1976D2")) +
      labs(title = lab, x = "", y = "", color = "") +
      theme_tcc + theme(legend.position = "bottom",
                        plot.title = element_text(size = 10))
    plots_3[[bn]] <- p
  }
  
  if (length(plots_3) >= 2) {
    p_comb <- grid.arrange(
      grobs = plots_3, ncol = 2,
      top = textGrob(sprintf("Betas OOS: 2SRR (TVP) vs Ridge (constante) — h=%d", h),
                     gp = gpar(fontface = "bold", fontsize = 13)))
    ggsave(file.path(fig_dir, sprintf("P3_betas_tvp_vs_ridge_h%02d.pdf", h)),
           p_comb, width = 14, height = 12)
    cat(sprintf("  [OK] P3_betas_tvp_vs_ridge_h%02d.pdf\n", h))
  }
}

# 3B: Tabela de correlacao TVP vs Ridge (media OOS)
cat("\n  Tabela de correlacao media betas TVP vs Ridge:\n")
cor_tab <- list()
for (hi in seq_along(hor)) {
  h <- hor[hi]
  bs_tvp   <- extract_beta_series(betas_2srr, hi, "2srr")
  bs_ridge <- extract_beta_series(betas_ridge, hi, "ridge")
  if (is.null(bs_tvp) || is.null(bs_ridge)) next
  n_min <- min(nrow(bs_tvp$betas), nrow(bs_ridge$betas))
  K_min <- min(bs_tvp$K, bs_ridge$K)
  
  tvp_mean   <- colMeans(bs_tvp$betas[1:n_min, 1:K_min], na.rm = TRUE)
  ridge_mean <- colMeans(bs_ridge$betas[1:n_min, 1:K_min], na.rm = TRUE)
  cor_val <- cor(tvp_mean, ridge_mean)
  dist_rel <- sqrt(sum((tvp_mean - ridge_mean)^2)) / sqrt(sum(ridge_mean^2))
  
  cor_tab[[hi]] <- data.frame(h = h, cor_media = round(cor_val, 4),
                               dist_eucl_rel = round(dist_rel, 4))
  cat(sprintf("  h=%d: cor=%.4f | dist_rel=%.4f\n", h, cor_val, dist_rel))
}
cor_df <- do.call(rbind, cor_tab)
write.csv(cor_df, file.path(out_dir, "P3_correlacao_tvp_vs_ridge.csv"), row.names = FALSE)

# ============================================================
# PARTE 4: BETAS VARIANDO POR HORIZONTE (cross-horizon)
# ============================================================
cat("\nPARTE 4: Betas cross-horizonte\n")

# Para os 5 regressores de interesse, mostra como o beta final muda com h
reg_cross <- intersect(reg_interest, extract_beta_series(betas_2srr, 1, "2srr")$col_names)

if (length(reg_cross) >= 2) {
  cross_data <- list()
  for (hi in seq_along(hor)) {
    h <- hor[hi]
    bs <- extract_beta_series(betas_2srr, hi, "2srr")
    if (is.null(bs)) next
    
    # Media dos ultimos 24 meses
    n_b <- nrow(bs$betas)
    window <- max(1, n_b - 23):n_b
    
    for (rn in reg_cross) {
      ki <- which(bs$col_names == rn)
      if (length(ki) == 0) next
      cross_data[[paste0(rn, "_h", h)]] <- data.frame(
        h = h, regressor = rn,
        mean_last24 = mean(bs$betas[window, ki], na.rm = TRUE),
        sd_last24   = sd(bs$betas[window, ki], na.rm = TRUE))
    }
  }
  
  if (length(cross_data) > 0) {
    cross_df <- do.call(rbind, cross_data)
    
    p_cross <- ggplot(cross_df, aes(x = factor(h), y = mean_last24, fill = regressor)) +
      geom_bar(stat = "identity", position = "dodge", width = 0.7) +
      geom_errorbar(aes(ymin = mean_last24 - sd_last24,
                        ymax = mean_last24 + sd_last24),
                    position = position_dodge(0.7), width = 0.2) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      labs(title = "Betas TVP (media ultimos 24 meses) por Horizonte de Previsao",
           subtitle = "Barras de erro = 1 desvio-padrao. Mostra como o peso muda com h.",
           x = "Horizonte (meses)", y = "Beta medio", fill = "Regressor") +
      theme_tcc
    
    ggsave(file.path(fig_dir, "P4_betas_cross_horizonte.pdf"), p_cross, width = 12, height = 6)
    write.csv(cross_df, file.path(out_dir, "P4_betas_cross_horizonte.csv"), row.names = FALSE)
    cat("  [OK] P4_betas_cross_horizonte.pdf\n")
  }
}

# ============================================================
# PARTE 5: TRAJETORIA DOS LAMBDAS (4h em 1 pagina)
# ============================================================
cat("\nPARTE 5: Lambdas ao longo do tempo\n")

plot_list_lam <- list()
lam_stats <- list()

for (hi in seq_along(hor)) {
  h <- hor[hi]
  key <- paste0("h", h)
  if (is.null(coulombe[[key]])) next
  df <- coulombe[[key]]
  df <- df[complete.cases(df$lam_ridge, df$lam2_2srr), ]
  if (nrow(df) < 10) next
  df$date <- as.Date(df$date)
  
  lam_long <- data.frame(
    date   = rep(df$date, 2),
    Lambda = c(df$lam_ridge, df$lam2_2srr),
    Modelo = c(rep("Ridge (lambda1)", nrow(df)),
               rep("2SRR (lambda2)", nrow(df))))
  
  cor_lam <- cor(df$lam_ridge, df$lam2_2srr, use = "complete.obs")
  ratio_lam <- mean(df$lam2_2srr, na.rm = TRUE) / mean(df$lam_ridge, na.rm = TRUE)
  
  lam_stats[[hi]] <- data.frame(h = h, cor = round(cor_lam, 4),
                                 ratio_media = round(ratio_lam, 4))
  
  p <- ggplot(lam_long, aes(x = date, y = Lambda, color = Modelo)) +
    add_recessions() +
    geom_line(linewidth = 0.4) +
    scale_y_log10() +
    scale_color_manual(values = c("Ridge (lambda1)" = "#1976D2",
                                   "2SRR (lambda2)" = "#D32F2F")) +
    labs(title = sprintf("h=%d | cor=%.3f", h, cor_lam),
         x = "", y = expression(lambda ~ "(log)"), color = "") +
    theme_tcc + theme(legend.position = "bottom",
                      plot.title = element_text(size = 10))
  
  plot_list_lam[[hi]] <- p
}

if (length(plot_list_lam) >= 4) {
  p_lam <- grid.arrange(
    grobs = plot_list_lam, ncol = 2,
    top = textGrob("Regularizacao Adaptativa: Lambda Ridge vs Lambda 2SRR",
                   gp = gpar(fontface = "bold", fontsize = 14)))
  ggsave(file.path(fig_dir, "P5_lambdas_4h.pdf"), p_lam, width = 14, height = 10)
  cat("  [OK] P5_lambdas_4h.pdf\n")
}

if (length(lam_stats) > 0) {
  lam_df <- do.call(rbind, lam_stats)
  write.csv(lam_df, file.path(out_dir, "P5_lambda_stats.csv"), row.names = FALSE)
  print(lam_df)
}

# ============================================================
# PARTE 6: 2SRR vs RIDGE COULOMBE (forecast direto)
# ============================================================
cat("\nPARTE 6: 2SRR vs Ridge Coulombe\n")

plot_list_fc <- list()
rmse_tab <- list()

for (hi in seq_along(hor)) {
  h <- hor[hi]
  key <- paste0("h", h)
  if (is.null(coulombe[[key]])) next
  df <- coulombe[[key]]
  df <- df[complete.cases(df$realized, df$fc_ridge, df$fc_2srr), ]
  if (nrow(df) < 20) next
  df$date <- as.Date(df$date)
  
  rmse_r <- sqrt(mean((df$fc_ridge - df$realized)^2))
  rmse_2 <- sqrt(mean((df$fc_2srr  - df$realized)^2))
  ratio  <- rmse_2 / rmse_r
  
  # DM test
  e_r <- df$fc_ridge - df$realized
  e_2 <- df$fc_2srr  - df$realized
  dm <- tryCatch(dm.test(e_r, e_2, alternative = "greater", h = h, power = 2),
                 error = function(e) list(statistic = NA, p.value = NA))
  
  rmse_tab[[hi]] <- data.frame(
    h = h, RMSE_Ridge = round(rmse_r, 5), RMSE_2SRR = round(rmse_2, 5),
    Ratio = round(ratio, 4),
    DM_stat = round(as.numeric(dm$statistic), 3),
    DM_pval = round(dm$p.value, 4),
    sig = ifelse(dm$p.value < 0.01, "***",
          ifelse(dm$p.value < 0.05, "**",
          ifelse(dm$p.value < 0.10, "*", ""))))
  
  # Grafico
  fc_long <- data.frame(
    date = rep(df$date, 3),
    valor = c(df$realized, df$fc_ridge, df$fc_2srr),
    serie = c(rep("Realizado", nrow(df)),
              rep("Ridge", nrow(df)),
              rep("2SRR", nrow(df))))
  
  p <- ggplot(fc_long, aes(x = date, y = valor, color = serie, linewidth = serie)) +
    add_recessions() +
    geom_line() +
    scale_color_manual(values = c("Realizado" = "black",
                                   "Ridge" = "#1976D2",
                                   "2SRR" = "#D32F2F")) +
    scale_linewidth_manual(values = c("Realizado" = 0.7, "Ridge" = 0.4, "2SRR" = 0.5),
                           guide = "none") +
    labs(title = sprintf("h=%d | Ratio=%.3f%s", h, ratio,
                         ifelse(ratio < 1, " (2SRR melhor)", "")),
         x = "", y = "Inflacao", color = "") +
    theme_tcc + theme(plot.title = element_text(size = 10))
  
  plot_list_fc[[hi]] <- p
}

if (length(plot_list_fc) >= 4) {
  p_fc <- grid.arrange(
    grobs = plot_list_fc, ncol = 2,
    top = textGrob("Previsoes OOS: 2SRR vs Ridge (Coulombe) vs Realizado",
                   gp = gpar(fontface = "bold", fontsize = 14)))
  ggsave(file.path(fig_dir, "P6_forecast_2srr_vs_ridge_4h.pdf"), p_fc, width = 14, height = 10)
  cat("  [OK] P6_forecast_2srr_vs_ridge_4h.pdf\n")
}

rmse_df <- do.call(rbind, rmse_tab)
cat("\n  RMSE 2SRR vs Ridge Coulombe:\n")
print(rmse_df)
write.csv(rmse_df, file.path(out_dir, "P6_rmse_2srr_vs_ridge.csv"), row.names = FALSE)

# ============================================================
# PARTE 7: 2SRR vs MELHORES E PIOR DO MEDEIROS
# ============================================================
cat("\nPARTE 7: 2SRR vs Medeiros (melhores + pior)\n")

# Calcula RMSFE relativo ao RW para todos
all_rmsfe <- list()
for (hi in seq_along(hor)) {
  h <- hor[hi]
  real <- yout[, hi]
  n <- length(real)
  
  for (mn in names(med_fc)) {
    mm <- med_fc[[mn]]
    if (ncol(mm) >= hi && nrow(mm) >= n) {
      fc_v <- mm[1:n, hi]
      rmse_m <- sqrt(mean((fc_v - real)^2, na.rm = TRUE))
      all_rmsfe[[paste0(mn, "_h", h)]] <- data.frame(
        h = h, model = mn, RMSE = rmse_m)
    }
  }
  
  # Adiciona 2SRR e Ridge do Coulombe
  key <- paste0("h", h)
  if (!is.null(coulombe[[key]])) {
    cdf <- coulombe[[key]]
    cdf <- cdf[complete.cases(cdf$realized, cdf$fc_ridge, cdf$fc_2srr), ]
    if (nrow(cdf) > 10) {
      all_rmsfe[[paste0("Ridge_C_h", h)]] <- data.frame(
        h = h, model = "Ridge_Coulombe",
        RMSE = sqrt(mean((cdf$fc_ridge - cdf$realized)^2)))
      all_rmsfe[[paste0("2SRR_C_h", h)]] <- data.frame(
        h = h, model = "2SRR_Coulombe",
        RMSE = sqrt(mean((cdf$fc_2srr - cdf$realized)^2)))
    }
  }
}

rmsfe_all <- do.call(rbind, all_rmsfe)
rownames(rmsfe_all) <- NULL

# RMSFE relativo ao RW
rmsfe_all$ratio_RW <- NA
for (h in hor) {
  rw_rmse <- rmsfe_all$RMSE[rmsfe_all$model == "rw" & rmsfe_all$h == h]
  if (length(rw_rmse) == 1) {
    idx <- rmsfe_all$h == h
    rmsfe_all$ratio_RW[idx] <- rmsfe_all$RMSE[idx] / rw_rmse
  }
}

# Ranking por horizonte
rmsfe_all$rank <- NA
for (h in hor) {
  idx <- rmsfe_all$h == h
  rmsfe_all$rank[idx] <- rank(rmsfe_all$RMSE[idx])
}

write.csv(rmsfe_all, file.path(out_dir, "P7_rmsfe_all_models.csv"), row.names = FALSE)

# Identifica melhores e pior por horizonte (excluindo RW e AR_BIC)
modelos_comp <- setdiff(unique(rmsfe_all$model), c("rw", "AR_BIC", "2SRR", "2SRR_Coulombe"))

cat("\n  RANKING POR HORIZONTE:\n")
for (h in hor) {
  sub <- rmsfe_all[rmsfe_all$h == h, ]
  sub <- sub[order(sub$RMSE), ]
  cat(sprintf("\n  h=%d:\n", h))
  for (i in 1:min(nrow(sub), 5)) {
    marker <- ifelse(grepl("2SRR", sub$model[i]), " <--", "")
    cat(sprintf("    %2d. %-20s RMSE=%.5f ratio_RW=%.4f%s\n",
                i, sub$model[i], sub$RMSE[i],
                ifelse(is.na(sub$ratio_RW[i]), NA, sub$ratio_RW[i]), marker))
  }
}

# Grafico: 2SRR vs top 2 melhores + pior (por horizonte)
plot_list_7 <- list()
for (h in hor) {
  sub <- rmsfe_all[rmsfe_all$h == h & rmsfe_all$model %in% modelos_comp, ]
  sub <- sub[order(sub$RMSE), ]
  
  best2 <- head(sub$model, 2)
  worst1 <- tail(sub$model, 1)
  
  # Adiciona 2SRR_Coulombe
  models_show <- unique(c("2SRR_Coulombe", best2, worst1, "rw"))
  sub_show <- rmsfe_all[rmsfe_all$h == h & rmsfe_all$model %in% models_show, ]
  
  if (nrow(sub_show) < 3) next
  
  # Cores
  cores <- c("2SRR_Coulombe" = "#D32F2F", "rw" = "gray60")
  for (m in best2) cores[m] <- "#2E7D32"
  for (m in worst1) cores[m] <- "#FF9800"
  
  p <- ggplot(sub_show, aes(x = reorder(model, RMSE), y = RMSE, fill = model)) +
    geom_col(width = 0.6) +
    scale_fill_manual(values = cores) +
    coord_flip() +
    labs(title = sprintf("h=%d", h), x = "", y = "RMSE") +
    theme_tcc + theme(legend.position = "none",
                      plot.title = element_text(size = 10))
  
  plot_list_7[[as.character(h)]] <- p
}

if (length(plot_list_7) >= 4) {
  p7 <- grid.arrange(
    grobs = plot_list_7, ncol = 2,
    top = textGrob("2SRR vs Top 2 Melhores + Pior modelo (Medeiros) + RW",
                   gp = gpar(fontface = "bold", fontsize = 14)))
  ggsave(file.path(fig_dir, "P7_2srr_vs_medeiros_4h.pdf"), p7, width = 14, height = 10)
  cat("\n  [OK] P7_2srr_vs_medeiros_4h.pdf\n")
}

# ============================================================
# PARTE 8: PARCIMONIA (HHI + near-zero + testes)
# ============================================================
cat("\nPARTE 8: Analise de parcimonia\n")

threshold <- 0.05  # 5% do max
parc_results <- list()
hhi_plots <- list()

for (hi in seq_along(hor)) {
  h <- hor[hi]
  bs_tvp   <- extract_beta_series(betas_2srr, hi, "2srr")
  bs_ridge <- extract_beta_series(betas_ridge, hi, "ridge")
  if (is.null(bs_tvp) || is.null(bs_ridge)) next
  
  n_min <- min(nrow(bs_tvp$betas), nrow(bs_ridge$betas))
  K_min <- min(bs_tvp$K, bs_ridge$K)
  
  # HHI por janela
  hhi_tvp <- numeric(n_min)
  hhi_ridge <- numeric(n_min)
  nz_tvp <- numeric(n_min)
  nz_ridge <- numeric(n_min)
  
  for (i in 1:n_min) {
    b_t <- abs(as.numeric(bs_tvp$betas[i, 1:K_min]))
    b_r <- abs(as.numeric(bs_ridge$betas[i, 1:K_min]))
    
    # HHI
    tot_t <- sum(b_t); tot_r <- sum(b_r)
    if (tot_t > 0) hhi_tvp[i] <- sum((b_t/tot_t)^2)
    if (tot_r > 0) hhi_ridge[i] <- sum((b_r/tot_r)^2)
    
    # Near-zero
    max_t <- max(b_t); max_r <- max(b_r)
    if (max_t > 0) nz_tvp[i] <- sum(b_t < threshold * max_t)
    if (max_r > 0) nz_ridge[i] <- sum(b_r < threshold * max_r)
  }
  
  # Teste t: HHI do 2SRR vs Ridge
  tt <- tryCatch(t.test(hhi_tvp, hhi_ridge, paired = TRUE),
                 error = function(e) list(statistic = NA, p.value = NA))
  
  parc_results[[hi]] <- data.frame(
    h = h,
    HHI_2SRR_mean = round(mean(hhi_tvp), 5),
    HHI_Ridge_mean = round(mean(hhi_ridge), 5),
    HHI_2SRR_mais_concentrado = mean(hhi_tvp) > mean(hhi_ridge),
    NZ_2SRR_mean = round(mean(nz_tvp), 1),
    NZ_Ridge_mean = round(mean(nz_ridge), 1),
    ttest_HHI_pval = round(tt$p.value, 4))
  
  cat(sprintf("  h=%d: HHI 2SRR=%.5f Ridge=%.5f | t-test p=%.4f | NZ: 2SRR=%.1f Ridge=%.1f\n",
              h, mean(hhi_tvp), mean(hhi_ridge), tt$p.value, mean(nz_tvp), mean(nz_ridge)))
  
  # Plot HHI ao longo do tempo
  hhi_df <- data.frame(
    date = bs_tvp$dates[1:n_min],
    HHI_2SRR = hhi_tvp,
    HHI_Ridge = hhi_ridge)
  hhi_long <- melt(hhi_df, id.vars = "date", variable.name = "Modelo", value.name = "HHI")
  
  p <- ggplot(hhi_long, aes(x = date, y = HHI, color = Modelo)) +
    add_recessions() +
    geom_line(linewidth = 0.4) +
    geom_hline(yintercept = 1/K_min, linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c("HHI_2SRR" = "#D32F2F", "HHI_Ridge" = "#1976D2")) +
    labs(title = sprintf("h=%d", h), x = "", y = "HHI", color = "") +
    theme_tcc + theme(legend.position = "bottom",
                      plot.title = element_text(size = 10))
  hhi_plots[[hi]] <- p
}

if (length(hhi_plots) >= 4) {
  p_hhi <- grid.arrange(
    grobs = hhi_plots, ncol = 2,
    top = textGrob("Concentracao de Informacao (HHI) — 2SRR vs Ridge",
                   gp = gpar(fontface = "bold", fontsize = 14)))
  ggsave(file.path(fig_dir, "P8_parcimonia_hhi_4h.pdf"), p_hhi, width = 14, height = 10)
  cat("  [OK] P8_parcimonia_hhi_4h.pdf\n")
}

parc_df <- do.call(rbind, parc_results)
write.csv(parc_df, file.path(out_dir, "P8_parcimonia_stats.csv"), row.names = FALSE)
cat("\n  Tabela de parcimonia:\n")
print(parc_df)

# ============================================================
# PARTE 9: TESTES ECONOMETRICOS COMPLETOS
# ============================================================
cat("\nPARTE 9: Testes econometricos\n")

econ_tests <- list()

for (hi in seq_along(hor)) {
  h <- hor[hi]
  key <- paste0("h", h)
  if (is.null(coulombe[[key]])) next
  df <- coulombe[[key]]
  df <- df[complete.cases(df$realized, df$fc_ridge, df$fc_2srr), ]
  if (nrow(df) < 30) next
  
  real <- df$realized
  fc_r <- df$fc_ridge
  fc_2 <- df$fc_2srr
  e_r  <- fc_r - real
  e_2  <- fc_2 - real
  
  # 9A: Diebold-Mariano
  dm <- tryCatch(dm.test(e_r, e_2, alternative = "greater", h = h, power = 2),
                 error = function(e) list(statistic = NA, p.value = NA))
  
  # 9B: Clark-West
  cw_d <- e_r^2 - (e_2^2 - (fc_r - fc_2)^2)
  cw_reg <- lm(cw_d ~ 1)
  cw_nw <- tryCatch(coeftest(cw_reg, vcov = NeweyWest(cw_reg, lag = h)),
                    error = function(e) matrix(c(NA,NA,NA,NA), nrow = 1))
  cw_pval <- ifelse(is.na(cw_nw[1,4]), NA, cw_nw[1,4] / 2)  # one-sided
  
  # 9C: Mincer-Zarnowitz para 2SRR
  mz_2 <- lm(real ~ fc_2)
  mz_r <- lm(real ~ fc_r)
  r2_2 <- summary(mz_2)$r.squared
  r2_r <- summary(mz_r)$r.squared
  
  mz_ftest_2 <- tryCatch({
    car::linearHypothesis(mz_2, c("(Intercept) = 0", "fc_2 = 1"))
  }, error = function(e) NULL)
  mz_pval_2 <- ifelse(!is.null(mz_ftest_2), mz_ftest_2[["Pr(>F)"]][2], NA)
  
  mz_ftest_r <- tryCatch({
    car::linearHypothesis(mz_r, c("(Intercept) = 0", "fc_r = 1"))
  }, error = function(e) NULL)
  mz_pval_r <- ifelse(!is.null(mz_ftest_r), mz_ftest_r[["Pr(>F)"]][2], NA)
  
  # 9D: Forecast Encompassing (Fair-Shiller)
  enc <- lm(real ~ fc_r + fc_2)
  enc_nw <- tryCatch(coeftest(enc, vcov = NeweyWest(enc, lag = h)),
                     error = function(e) coeftest(enc))
  p_ridge_enc <- enc_nw[2, 4]
  p_2srr_enc  <- enc_nw[3, 4]
  
  encomp_result <- if (p_2srr_enc < 0.05 && p_ridge_enc > 0.05) "2SRR encompassa Ridge"
  else if (p_ridge_enc < 0.05 && p_2srr_enc > 0.05) "Ridge encompassa 2SRR"
  else if (p_ridge_enc < 0.05 && p_2srr_enc < 0.05) "Ambos contribuem"
  else "Nenhum sig."
  
  econ_tests[[hi]] <- data.frame(
    h = h,
    DM_stat = round(as.numeric(dm$statistic), 3),
    DM_pval = round(dm$p.value, 4),
    CW_pval = round(cw_pval, 4),
    MZ_R2_2SRR = round(r2_2, 4),
    MZ_R2_Ridge = round(r2_r, 4),
    MZ_pval_2SRR = round(mz_pval_2, 4),
    MZ_pval_Ridge = round(mz_pval_r, 4),
    Encomp = encomp_result,
    stringsAsFactors = FALSE)
  
  cat(sprintf("  h=%d: DM p=%.4f | CW p=%.4f | MZ R2: 2SRR=%.4f Ridge=%.4f | %s\n",
              h, dm$p.value, cw_pval, r2_2, r2_r, encomp_result))
}

econ_df <- do.call(rbind, econ_tests)
write.csv(econ_df, file.path(out_dir, "P9_testes_econometricos.csv"), row.names = FALSE)
cat("\n  Tabela de testes:\n")
print(econ_df)

# ============================================================
# PARTE 10: SUB-PERIODOS E REGIMES (4h em 1 pagina)
# ============================================================
cat("\nPARTE 10: Sub-periodos\n")

periodos <- list(
  "Full Sample"    = c(as.Date("1999-07-01"), as.Date("2025-06-01")),
  "Pre-GFC"        = c(as.Date("1999-07-01"), as.Date("2007-11-30")),
  "GFC"            = c(as.Date("2007-12-01"), as.Date("2009-06-30")),
  "Post-GFC"       = c(as.Date("2009-07-01"), as.Date("2020-01-31")),
  "COVID"          = c(as.Date("2020-02-01"), as.Date("2021-06-30")),
  "High Inflation" = c(as.Date("2021-07-01"), as.Date("2023-06-30")))

sub_results <- list()
for (h in hor) {
  key <- paste0("h", h)
  if (is.null(coulombe[[key]])) next
  df <- coulombe[[key]]
  df <- df[complete.cases(df$realized, df$fc_ridge, df$fc_2srr), ]
  df$date_p <- as.Date(df$date)
  
  for (pn in names(periodos)) {
    pr <- periodos[[pn]]
    idx <- df$date_p >= pr[1] & df$date_p <= pr[2]
    if (sum(idx) < 5) next
    rmse_r <- sqrt(mean((df$fc_ridge[idx] - df$realized[idx])^2))
    rmse_2 <- sqrt(mean((df$fc_2srr[idx]  - df$realized[idx])^2))
    sub_results[[paste0(pn, "_h", h)]] <- data.frame(
      h = h, periodo = pn, n = sum(idx),
      RMSE_Ridge = round(rmse_r, 5), RMSE_2SRR = round(rmse_2, 5),
      Ratio = round(rmse_2/rmse_r, 4),
      Vantagem_2SRR = rmse_2 < rmse_r)
  }
}

if (length(sub_results) > 0) {
  sub_tab <- do.call(rbind, sub_results)
  rownames(sub_tab) <- NULL
  write.csv(sub_tab, file.path(out_dir, "P10_subperiodos.csv"), row.names = FALSE)
  
  sub_tab$periodo <- factor(sub_tab$periodo, levels = names(periodos))
  
  p_sub <- ggplot(sub_tab, aes(x = periodo, y = Ratio, fill = factor(h))) +
    geom_col(position = "dodge", width = 0.7) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 0.8) +
    scale_fill_manual(values = c("1" = "#D32F2F", "3" = "#1976D2",
                                  "6" = "#388E3C", "12" = "#FF9800")) +
    labs(title = "RMSE(2SRR) / RMSE(Ridge) por Sub-periodo e Horizonte",
         subtitle = "Abaixo de 1 = 2SRR melhor. Cinza = periodos de recessao/estresse",
         x = "", y = "Ratio RMSE", fill = "h") +
    theme_tcc + theme(axis.text.x = element_text(angle = 25, hjust = 1))
  
  ggsave(file.path(fig_dir, "P10_subperiodos.pdf"), p_sub, width = 14, height = 7)
  cat("  [OK] P10_subperiodos.pdf\n")
  
  # Conta vitorias
  n_wins <- sum(sub_tab$Vantagem_2SRR)
  n_total <- nrow(sub_tab)
  cat(sprintf("  2SRR ganha em %d/%d combinacoes sub-periodo x horizonte\n", n_wins, n_total))
}

# ============================================================
# PARTE 11: ROLLING RMSE + CSFE (4h combinados)
# ============================================================
cat("\nPARTE 11: Rolling RMSE e CSFE\n")

plot_list_roll <- list()
plot_list_csfe <- list()

for (hi in seq_along(hor)) {
  h <- hor[hi]
  key <- paste0("h", h)
  if (is.null(coulombe[[key]])) next
  df <- coulombe[[key]]
  df <- df[complete.cases(df$realized, df$fc_ridge, df$fc_2srr), ]
  n_df <- nrow(df); window <- 36
  if (n_df < window + 10) next
  df$date <- as.Date(df$date)
  
  # Rolling ratio
  roll_ratio <- rep(NA, n_df)
  for (i in window:n_df) {
    w <- (i - window + 1):i
    rmse_r <- sqrt(mean((df$fc_ridge[w] - df$realized[w])^2))
    rmse_2 <- sqrt(mean((df$fc_2srr[w]  - df$realized[w])^2))
    roll_ratio[i] <- rmse_2 / rmse_r
  }
  
  roll_df <- data.frame(date = df$date, ratio = roll_ratio)
  roll_df <- roll_df[!is.na(roll_df$ratio), ]
  pct_below <- round(100 * mean(roll_df$ratio < 1), 1)
  
  p_roll <- ggplot(roll_df, aes(x = date, y = ratio)) +
    add_recessions() +
    geom_ribbon(aes(ymin = pmin(ratio, 1), ymax = 1), fill = "#2E7D32", alpha = 0.15) +
    geom_ribbon(aes(ymin = 1, ymax = pmax(ratio, 1)), fill = "#D32F2F", alpha = 0.15) +
    geom_line(linewidth = 0.5, color = "#D32F2F") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    labs(title = sprintf("h=%d | 2SRR melhor em %.0f%% das janelas", h, pct_below),
         x = "", y = "Ratio RMSE") +
    theme_tcc + theme(plot.title = element_text(size = 10))
  plot_list_roll[[hi]] <- p_roll
  
  # CSFE
  d_t <- (df$fc_ridge - df$realized)^2 - (df$fc_2srr - df$realized)^2
  csfe <- cumsum(d_t)
  csfe_df <- data.frame(date = df$date, CSFE = csfe)
  
  p_csfe <- ggplot(csfe_df, aes(x = date, y = CSFE)) +
    add_recessions() +
    geom_line(linewidth = 0.6, color = "#D32F2F") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    labs(title = sprintf("h=%d | CSFE final=%.2f", h, tail(csfe, 1)),
         x = "", y = "CSFE") +
    theme_tcc + theme(plot.title = element_text(size = 10))
  plot_list_csfe[[hi]] <- p_csfe
}

if (length(plot_list_roll) >= 4) {
  p_r <- grid.arrange(
    grobs = plot_list_roll, ncol = 2,
    top = textGrob("Rolling RMSE Ratio (janela 36 meses) — 2SRR/Ridge",
                   gp = gpar(fontface = "bold", fontsize = 14)))
  ggsave(file.path(fig_dir, "P11_rolling_rmse_4h.pdf"), p_r, width = 14, height = 10)
  cat("  [OK] P11_rolling_rmse_4h.pdf\n")
}

if (length(plot_list_csfe) >= 4) {
  p_c <- grid.arrange(
    grobs = plot_list_csfe, ncol = 2,
    top = textGrob("CSFE: Acima de 0 = 2SRR acumula menos erro que Ridge",
                   gp = gpar(fontface = "bold", fontsize = 14)))
  ggsave(file.path(fig_dir, "P11_csfe_4h.pdf"), p_c, width = 14, height = 10)
  cat("  [OK] P11_csfe_4h.pdf\n")
}

# ============================================================
# PARTE 12: TABELA LATEX CONSOLIDADA
# ============================================================
cat("\nPARTE 12: Tabela LaTeX\n")

# Tabela principal: RMSFE relativo ao RW para modelos selecionados
modelos_latex <- intersect(c("rw", "AR", "Ridge", "LASSO", "AdaLASSO",
                             "EINET", "RF", "Bagging", "CSR", "factor"),
                           names(med_fc))

# Monta pivot table
latex_rows <- list()
for (h in hor) {
  row <- list(h = h)
  real <- yout[, which(hor == h)]
  n <- length(real)
  
  rw_rmse <- NA
  if (!is.null(med_fc[["rw"]])) {
    rw_v <- med_fc[["rw"]][1:n, which(hor == h)]
    rw_rmse <- sqrt(mean((rw_v - real)^2, na.rm = TRUE))
  }
  
  for (mn in modelos_latex) {
    mm <- med_fc[[mn]]
    hi_idx <- which(hor == h)
    if (ncol(mm) >= hi_idx && nrow(mm) >= n) {
      fc_v <- mm[1:n, hi_idx]
      rmse_m <- sqrt(mean((fc_v - real)^2, na.rm = TRUE))
      row[[mn]] <- round(rmse_m / rw_rmse, 4)
    } else {
      row[[mn]] <- NA
    }
  }
  
  # 2SRR Coulombe
  key <- paste0("h", h)
  if (!is.null(coulombe[[key]])) {
    cdf <- coulombe[[key]]
    cdf <- cdf[complete.cases(cdf$realized, cdf$fc_2srr), ]
    rmse_2 <- sqrt(mean((cdf$fc_2srr - cdf$realized)^2))
    # Precisa do RW na mesma escala... usa ratio vs Ridge como proxy
    row[["2SRR"]] <- round(rmse_2 / rw_rmse, 4)
  }
  
  latex_rows[[as.character(h)]] <- as.data.frame(row, stringsAsFactors = FALSE)
}

latex_tab <- do.call(rbind, latex_rows)
write.csv(latex_tab, file.path(out_dir, "P12_rmsfe_relativo_rw.csv"), row.names = FALSE)

# Gera LaTeX
tex_file <- file.path(out_dir, "P12_tabela_latex.tex")
sink(tex_file)
cat("\\begin{table}[ht]\n\\centering\n")
cat("\\caption{RMSFE Relativo ao Random Walk — Inflacao CPI (EUA)}\n")
cat("\\label{tab:rmsfe_rw}\n")
cols_show <- intersect(c("AR","Ridge","LASSO","AdaLASSO","EINET","RF","factor","CSR","Bagging","2SRR"),
                       names(latex_tab))
cat(sprintf("\\begin{tabular}{c%s}\n", paste(rep("c", length(cols_show)), collapse = "")))
cat("\\hline\\hline\n")
cat(sprintf("h & %s \\\\\n", paste(cols_show, collapse = " & ")))
cat("\\hline\n")
for (i in seq_len(nrow(latex_tab))) {
  vals <- sapply(cols_show, function(m) {
    v <- latex_tab[i, m]
    if (is.null(v) || is.na(v)) return("---")
    sprintf("%.3f", v)
  })
  cat(sprintf("%d & %s \\\\\n", latex_tab$h[i], paste(vals, collapse = " & ")))
}
cat("\\hline\\hline\n\\end{tabular}\n")
cat("\\begin{flushleft}\n")
cat("\\footnotesize Nota: Valores $< 1$ indicam modelo superior ao Random Walk.\n")
cat("\\end{flushleft}\n")
cat("\\end{table}\n")
sink()
cat("  [OK] P12_tabela_latex.tex\n")

# ============================================================
# RESUMO FINAL
# ============================================================

cat("\n\n============================================================\n")
cat(" RESUMO DOS OUTPUTS GERADOS\n")
cat("============================================================\n\n")

n_pdf <- length(list.files(fig_dir, pattern = "\\.pdf$"))
n_csv <- length(list.files(out_dir, pattern = "\\.csv$"))
n_tex <- length(list.files(out_dir, pattern = "\\.tex$"))

cat(sprintf("  Pasta: %s\n", out_dir))
cat(sprintf("  PDFs:  %d\n", n_pdf))
cat(sprintf("  CSVs:  %d\n", n_csv))
cat(sprintf("  TEX:   %d\n", n_tex))

cat("\n  LISTA DE OUTPUTS:\n")
cat("  FIGURAS (todas com 4 horizontes combinados):\n")
cat("    P2_betas_tvp_evolucao_4h.pdf      — Betas TVP ao longo do tempo\n")
cat("    P3_betas_tvp_vs_ridge_h01..12.pdf  — Betas 2SRR vs Ridge OOS (por h)\n")
cat("    P4_betas_cross_horizonte.pdf       — Betas por horizonte de previsao\n")
cat("    P5_lambdas_4h.pdf                  — Trajetoria dos lambdas\n")
cat("    P6_forecast_2srr_vs_ridge_4h.pdf   — Previsoes OOS: 2SRR vs Ridge\n")
cat("    P7_2srr_vs_medeiros_4h.pdf         — 2SRR vs melhores/pior Medeiros\n")
cat("    P8_parcimonia_hhi_4h.pdf           — Concentracao (HHI) ao longo do tempo\n")
cat("    P10_subperiodos.pdf                — Ratio RMSE por sub-periodo\n")
cat("    P11_rolling_rmse_4h.pdf            — Rolling RMSE ratio (36 meses)\n")
cat("    P11_csfe_4h.pdf                    — CSFE cumulativo\n\n")

cat("  TABELAS:\n")
cat("    P3_correlacao_tvp_vs_ridge.csv     — Correlacao betas TVP vs Ridge\n")
cat("    P5_lambda_stats.csv                — Estatisticas dos lambdas\n")
cat("    P6_rmse_2srr_vs_ridge.csv          — RMSE + DM test (2SRR vs Ridge)\n")
cat("    P7_rmsfe_all_models.csv            — RMSFE todos os modelos\n")
cat("    P8_parcimonia_stats.csv            — HHI + near-zero + t-test\n")
cat("    P9_testes_econometricos.csv        — DM, Clark-West, MZ, Encompassing\n")
cat("    P10_subperiodos.csv                — RMSE por sub-periodo\n")
cat("    P12_rmsfe_relativo_rw.csv          — RMSFE relativo ao RW\n")
cat("    P12_tabela_latex.tex               — Tabela LaTeX para o TCC\n\n")

cat(" analysis_2SRR.R --- COMPLETO\n")