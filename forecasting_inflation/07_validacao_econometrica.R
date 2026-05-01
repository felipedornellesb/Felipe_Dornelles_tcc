# ============================================================
# 07_validacao_econometrica.R  v3.0
#
# PRE-REQUISITO: 06_coulombe_2SRR_pipeline.R v8.0 ja rodou
#   com fc_2srr preenchido (nao pode ser NA)
#
# CONTEUDO:
#   PARTE 1: Carrega todos os forecasts (Coulombe + Medeiros)
#   PARTE 2: Metricas basicas (RMSE, MAE, MAPE) — tabela
#   PARTE 3: Diebold-Mariano pairwise (2SRR vs TODOS)
#   PARTE 4: Clark-West (nested model test)
#   PARTE 5: Mincer-Zarnowitz (eficiencia)
#   PARTE 6: Forecast Encompassing (Fair-Shiller)
#   PARTE 7: Giacomini-White (conditional predictive ability)
#   PARTE 8: CSSED + Fluctuation Test (Giacomini-Rossi)
#   PARTE 9: Model Confidence Set (MCS)
#   PARTE 10: Analise TVP dos betas (todos horizontes)
#   PARTE 11: Analise por sub-amostras (pre/pos COVID)
#   PARTE 12: Graficos de qualidade publicacao
#   PARTE 13: Tabela LaTeX consolidada
# ============================================================

rm(list = ls())
setwd("~/tcc/Felipe_Dornelles_tcc/forecasting_inflation")

# --- Pacotes ---
pkgs <- c("forecast", "lmtest", "sandwich", "ggplot2", "reshape2",
          "gridExtra", "car", "xtable", "scales")
new <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if (length(new)) install.packages(new, repos = "https://cran.r-project.org")
invisible(lapply(pkgs, library, character.only = TRUE))

mcs_ok <- tryCatch({ library(MCS); TRUE }, error = function(e) {
  tryCatch({
    install.packages("MCS", repos = "https://cran.r-project.org")
    library(MCS); TRUE
  }, error = function(e2) FALSE)
})

hor <- c(1, 3, 6, 12)
dir.create("results", showWarnings = FALSE)
dir.create("results/figures", showWarnings = FALSE)

cat("  VALIDACAO ECONOMETRICA COMPLETA v3.0\n")

# ============================================================
# PARTE 1: CARREGA TODOS OS FORECASTS
# ============================================================
cat("=== PARTE 1: Carregando forecasts ===\n\n")

# 1A. Forecasts Coulombe (2SRR + Ridge)
coulombe <- list()
for (h in hor) {
  fname <- sprintf("forecasts/coulombe_fc_h%02d.csv", h)
  if (file.exists(fname)) {
    df <- read.csv(fname, stringsAsFactors = FALSE)
    df <- df[complete.cases(df$realized, df$fc_ridge, df$fc_2srr), ]
    coulombe[[paste0("h", h)]] <- df
    cat(sprintf("  Coulombe h=%2d: %d obs validas\n", h, nrow(df)))
  } else {
    cat(sprintf("  [AVISO] %s nao encontrado\n", fname))
  }
}

# Verifica se 2SRR tem dados
has_2srr <- any(sapply(coulombe, function(df) {
  !all(is.na(df$fc_2srr))
}))
if (!has_2srr) {
  stop(paste0(
    "\n*** ERRO CRITICO ***\n",
    "fc_2srr esta 100% NA em todos os horizontes.\n",
    "Voce precisa re-rodar o pipeline v8.0 (com fGarch instalado)\n",
    "ANTES de rodar este script de validacao.\n",
    "Execute: install.packages('fGarch')\n",
    "E depois re-rode o 06_coulombe_2SRR_pipeline.R v8.0\n"
  ))
}

# 1B. Forecasts Medeiros
cat("\n  Modelos Medeiros:\n")
med_forecasts <- list()
ffiles <- list.files("forecasts", pattern = "^for_.*\\.rda$", full.names = TRUE)
load("forecasts/yout.rda")  # realized values

for (ff in ffiles) {
  mn <- gsub("^for_|\\.rda$", "", basename(ff))
  tryCatch({
    env <- new.env()
    load(ff, envir = env)
    obj <- get(ls(env)[1], envir = env)
    if (is.matrix(obj) || is.data.frame(obj)) {
      med_forecasts[[mn]] <- as.matrix(obj)
      cat(sprintf("    %-20s: %d x %d\n", mn, nrow(obj), ncol(obj)))
    }
  }, error = function(e) NULL)
}
cat(sprintf("\n  Total: %d modelos Medeiros + 2 Coulombe (Ridge, 2SRR)\n\n",
            length(med_forecasts)))

# ============================================================
# PARTE 2: METRICAS BASICAS
# ============================================================
cat("=== PARTE 2: Metricas por horizonte ===\n\n")

all_metrics <- list()

for (h in hor) {
  hi <- which(hor == h)
  key <- paste0("h", h)
  if (is.null(coulombe[[key]])) next

  df <- coulombe[[key]]
  n_oos <- nrow(df)
  real <- df$realized

  # Todos os modelos num data.frame
  models <- data.frame(
    Ridge_Coulombe = df$fc_ridge,
    SRR2 = df$fc_2srr
  )

  # Adiciona Medeiros
  for (mn in names(med_forecasts)) {
    mm <- med_forecasts[[mn]]
    if (ncol(mm) >= hi && nrow(mm) >= n_oos) {
      fc <- mm[1:n_oos, hi]
      if (sum(!is.na(fc)) > n_oos * 0.5) {
        models[[mn]] <- fc
      }
    }
  }

  # Calcula metricas
  metrics <- data.frame(
    model = names(models),
    RMSE  = sapply(models, function(fc) sqrt(mean((fc - real)^2, na.rm = TRUE))),
    MAE   = sapply(models, function(fc) mean(abs(fc - real), na.rm = TRUE)),
    MAPE  = sapply(models, function(fc) mean(abs((fc - real) / real), na.rm = TRUE) * 100),
    n_valid = sapply(models, function(fc) sum(!is.na(fc) & !is.na(real))),
    stringsAsFactors = FALSE
  )
  metrics <- metrics[order(metrics$RMSE), ]
  metrics$rank <- seq_len(nrow(metrics))
  metrics$h <- h

  # Ratio relativo ao 2SRR
  rmse_2srr <- metrics$RMSE[metrics$model == "SRR2"]
  metrics$ratio_vs_2SRR <- metrics$RMSE / rmse_2srr

  all_metrics[[key]] <- metrics

  cat(sprintf("  h = %d (n = %d):\n", h, n_oos))
  cat(sprintf("  %-25s %8s %8s %8s %6s\n", "Modelo", "RMSE", "MAE", "MAPE%", "Rank"))
  cat(paste(rep("-", 60), collapse = ""), "\n")
  for (i in seq_len(nrow(metrics))) {
    star <- ifelse(metrics$model[i] == "SRR2", " <--", "")
    cat(sprintf("  %-25s %8.4f %8.4f %8.2f %6d%s\n",
                metrics$model[i], metrics$RMSE[i], metrics$MAE[i],
                metrics$MAPE[i], metrics$rank[i], star))
  }
  cat("\n")
}

# ============================================================
# PARTE 3: DIEBOLD-MARIANO PAIRWISE (2SRR vs TODOS)
# ============================================================
cat("=== PARTE 3: Diebold-Mariano (2SRR vs cada modelo) ===\n\n")

dm_results <- list()

for (h in hor) {
  hi <- which(hor == h)
  key <- paste0("h", h)
  if (is.null(coulombe[[key]])) next

  df <- coulombe[[key]]
  real <- df$realized
  e_2srr <- df$fc_2srr - real
  n_oos <- nrow(df)

  cat(sprintf("  h = %d:\n", h))
  cat(sprintf("  %-25s %8s %8s %12s\n", "Modelo", "DM_stat", "p-valor", "Resultado"))
  cat(paste(rep("-", 60), collapse = ""), "\n")

  # vs Ridge Coulombe
  all_models <- list(Ridge_Coulombe = df$fc_ridge - real)

  # vs Medeiros
  for (mn in names(med_forecasts)) {
    mm <- med_forecasts[[mn]]
    if (ncol(mm) >= hi && nrow(mm) >= n_oos) {
      e_m <- mm[1:n_oos, hi] - real
      if (sum(!is.na(e_m)) > n_oos * 0.5) {
        all_models[[mn]] <- e_m
      }
    }
  }

  for (mn in names(all_models)) {
    e_other <- all_models[[mn]]
    valid <- !is.na(e_2srr) & !is.na(e_other)

    dm <- tryCatch(
      dm.test(e1 = e_other[valid], e2 = e_2srr[valid],
              h = h, alternative = "greater", power = 2),
      error = function(e) NULL)

    if (!is.null(dm)) {
      sig <- ifelse(dm$p.value < 0.01, "***",
             ifelse(dm$p.value < 0.05, "**",
             ifelse(dm$p.value < 0.10, "*", "")))
      result <- ifelse(dm$p.value < 0.05, "2SRR GANHA",
                ifelse(dm$p.value > 0.95, "2SRR PERDE", "Empate"))
      cat(sprintf("  %-25s %8.3f %8.4f %12s %s\n",
                  mn, dm$statistic, dm$p.value, result, sig))

      dm_results[[paste0(key, "_", mn)]] <- data.frame(
        h = h, model = mn, DM_stat = dm$statistic,
        p_value = dm$p.value, result = result)
    }
  }
  cat("\n")
}

# ============================================================
# PARTE 4: CLARK-WEST TEST
# ============================================================
cat("=== PARTE 4: Clark-West (Ridge restrito vs 2SRR irrestrito) ===\n\n")

for (h in hor) {
  key <- paste0("h", h)
  if (is.null(coulombe[[key]])) next
  df <- coulombe[[key]]
  real <- df$realized
  fc1  <- df$fc_ridge  # modelo restrito (constante)
  fc2  <- df$fc_2srr   # modelo irrestrito (TVP)

  e1 <- fc1 - real
  e2 <- fc2 - real

  # Clark-West: d_t = e1^2 - [e2^2 - (fc1 - fc2)^2]
  cw_d <- e1^2 - (e2^2 - (fc1 - fc2)^2)
  valid <- !is.na(cw_d)
  cw_d <- cw_d[valid]

  if (length(cw_d) > 20) {
    # Regressao de cw_d contra constante com HAC
    cw_reg <- lm(cw_d ~ 1)
    cw_nw <- tryCatch(
      coeftest(cw_reg, vcov = NeweyWest(cw_reg, lag = h)),
      error = function(e) coeftest(cw_reg))

    cw_stat <- cw_nw[1, 3]  # t-statistic
    cw_pval <- cw_nw[1, 4] / 2  # one-sided

    sig <- ifelse(cw_pval < 0.01, "***",
           ifelse(cw_pval < 0.05, "**",
           ifelse(cw_pval < 0.10, "*", "")))

    cat(sprintf("  h=%2d: CW_stat=%.3f p=%.4f %s %s\n",
                h, cw_stat, cw_pval, sig,
                ifelse(cw_pval < 0.05,
                       "[2SRR sig. melhor que Ridge constante]",
                       "[Nao rejeita H0]")))
  }
}
cat("\n")

# ============================================================
# PARTE 5: MINCER-ZARNOWITZ
# ============================================================
cat("=== PARTE 5: Mincer-Zarnowitz ===\n\n")

for (h in hor) {
  key <- paste0("h", h)
  if (is.null(coulombe[[key]])) next
  df <- coulombe[[key]]
  real <- df$realized

  cat(sprintf("  h = %d:\n", h))
  for (mn in c("Ridge", "2SRR")) {
    fc_use <- if (mn == "Ridge") df$fc_ridge else df$fc_2srr
    mz <- lm(real ~ fc_use)
    a <- coef(mz)[1]; b <- coef(mz)[2]
    r2 <- summary(mz)$r.squared

    ftest <- tryCatch(
      car::linearHypothesis(mz, c("(Intercept) = 0", "fc_use = 1")),
      error = function(e) NULL)

    nw <- tryCatch(
      coeftest(mz, vcov = NeweyWest(mz, lag = h)),
      error = function(e) NULL)

    cat(sprintf("    %s: a=%.4f b=%.4f R2=%.4f", mn, a, b, r2))
    if (!is.null(ftest)) {
      fp <- ftest[["Pr(>F)"]][2]
      cat(sprintf(" | F(a=0,b=1): p=%.4f %s",
                  fp, ifelse(fp > 0.05, "[Eficiente]", "[Vies]")))
    }
    cat("\n")
  }
  cat("\n")
}

# ============================================================
# PARTE 6: FORECAST ENCOMPASSING
# ============================================================
cat("=== PARTE 6: Forecast Encompassing (Fair-Shiller) ===\n\n")

for (h in hor) {
  key <- paste0("h", h)
  if (is.null(coulombe[[key]])) next
  df <- coulombe[[key]]
  real <- df$realized

  enc <- lm(real ~ df$fc_ridge + df$fc_2srr)
  enc_nw <- tryCatch(
    coeftest(enc, vcov = NeweyWest(enc, lag = h)),
    error = function(e) coeftest(enc))

  p_ridge <- enc_nw[2, 4]
  p_2srr  <- enc_nw[3, 4]

  interp <- if (p_2srr < 0.05 && p_ridge > 0.05) "2SRR ENCOMPASSA Ridge"
  else if (p_ridge < 0.05 && p_2srr > 0.05) "Ridge encompassa 2SRR"
  else if (p_ridge < 0.05 && p_2srr < 0.05) "Ambos contribuem"
  else "Nenhum significativo"

  cat(sprintf("  h=%2d: b_Ridge=%.3f(p=%.3f) b_2SRR=%.3f(p=%.3f) -> %s\n",
              h, enc_nw[2,1], p_ridge, enc_nw[3,1], p_2srr, interp))
}
cat("\n")

# ============================================================
# PARTE 7: GIACOMINI-WHITE TEST
# ============================================================
cat("=== PARTE 7: Giacomini-White (Conditional Predictive Ability) ===\n\n")

for (h in hor) {
  key <- paste0("h", h)
  if (is.null(coulombe[[key]])) next
  df <- coulombe[[key]]
  real <- df$realized

  e_r <- (df$fc_ridge - real)^2
  e_2 <- (df$fc_2srr  - real)^2
  d_t <- e_r - e_2  # positivo = 2SRR melhor

  n <- length(d_t)
  if (n < h + 10) next

  # Instrumentos: constante + lag de d_t
  # GW regressa d_t em h_t-1 (instrumentos disponiveis em t-1)
  d_lag <- c(rep(NA, h), d_t[1:(n - h)])
  valid <- !is.na(d_lag) & !is.na(d_t)

  if (sum(valid) < 20) next

  gw_reg <- lm(d_t[valid] ~ d_lag[valid])
  gw_nw  <- tryCatch(
    coeftest(gw_reg, vcov = NeweyWest(gw_reg, lag = h)),
    error = function(e) coeftest(gw_reg))

  # Teste F conjunto: todos os coeficientes = 0
  gw_f <- tryCatch(
    car::linearHypothesis(gw_reg,
      c("(Intercept) = 0", "d_lag[valid] = 0"),
      vcov = NeweyWest(gw_reg, lag = h)),
    error = function(e) NULL)

  if (!is.null(gw_f)) {
    f_stat <- gw_f[["F"]][2]
    f_pval <- gw_f[["Pr(>F)"]][2]
    cat(sprintf("  h=%2d: GW F=%.3f p=%.4f %s\n",
                h, f_stat, f_pval,
                ifelse(f_pval < 0.05,
                       "[Rejeita H0: diferenca e previsivel]",
                       "[Nao rejeita: diferenca constante no tempo]")))
  }
}
cat("\n")

# ============================================================
# PARTE 8: CSSED + FLUCTUATION TEST
# ============================================================
cat("=== PARTE 8: CSSED e Fluctuation Test ===\n\n")

for (h in hor) {
  key <- paste0("h", h)
  if (is.null(coulombe[[key]])) next
  df <- coulombe[[key]]
  real <- df$realized

  d_t   <- (df$fc_ridge - real)^2 - (df$fc_2srr - real)^2
  cssed <- cumsum(d_t)
  n     <- length(d_t)

  pct_wins <- mean(d_t > 0) * 100

  # Fluctuation test: DM rolling com janela m
  m <- floor(n * 0.3)  # 30% da amostra
  if (m > 30) {
    dm_roll <- rep(NA, n)
    for (i in m:n) {
      window <- (i - m + 1):i
      e1_w <- (df$fc_ridge[window] - real[window])
      e2_w <- (df$fc_2srr[window]  - real[window])
      dm_w <- tryCatch({
        tt <- dm.test(e1_w, e2_w, h = h, power = 2)
        tt$statistic
      }, error = function(e) NA)
      dm_roll[i] <- dm_w
    }

    # Sup statistic
    sup_dm <- max(abs(dm_roll), na.rm = TRUE)
    # Sob H0, cv ~3.39 para 10% e ~3.76 para 5% (Giacomini-Rossi Table 1)
    cat(sprintf("  h=%2d: CSSED_final=%.2f | 2SRR ganha %.0f%% | SupDM=%.2f %s\n",
                h, tail(cssed, 1), pct_wins, sup_dm,
                ifelse(sup_dm > 3.39, "[Instabilidade detectada]",
                       "[Ranking estavel]")))
  } else {
    cat(sprintf("  h=%2d: CSSED_final=%.2f | 2SRR ganha %.0f%%\n",
                h, tail(cssed, 1), pct_wins))
  }

  # Salva CSSED
  df_cssed <- data.frame(date = df$date, cssed = cssed, d_t = d_t)
  write.csv(df_cssed, sprintf("results/cssed_h%02d.csv", h), row.names = FALSE)
}
cat("\n")

# ============================================================
# PARTE 9: MODEL CONFIDENCE SET
# ============================================================
if (mcs_ok) {
  cat("=== PARTE 9: Model Confidence Set ===\n\n")

  for (h in hor) {
    hi <- which(hor == h)
    key <- paste0("h", h)
    if (is.null(coulombe[[key]])) next
    df <- coulombe[[key]]
    real <- df$realized
    n_oos <- nrow(df)

    loss_list <- list(
      Ridge_C = (df$fc_ridge - real)^2,
      SRR2    = (df$fc_2srr  - real)^2)

    for (mn in names(med_forecasts)) {
      mm <- med_forecasts[[mn]]
      if (ncol(mm) >= hi && nrow(mm) >= n_oos) {
        fc <- mm[1:n_oos, hi]
        if (sum(!is.na(fc)) > n_oos * 0.5)
          loss_list[[mn]] <- (fc - real)^2
      }
    }

    loss_mat <- do.call(cbind, loss_list)
    loss_mat <- loss_mat[complete.cases(loss_mat), , drop = FALSE]

    if (ncol(loss_mat) >= 2 && nrow(loss_mat) >= 30) {
      cat(sprintf("  MCS h=%d: %d modelos, %d obs\n",
                  h, ncol(loss_mat), nrow(loss_mat)))
      mcs_r <- tryCatch(
        MCSprocedure(Loss = loss_mat, alpha = 0.10, B = 5000,
                     statistic = "Tmax"),
        error = function(e) {
          tryCatch(
            MCSprocedure(Loss = loss_mat, alpha = 0.10, B = 5000,
                         statistic = "TR"),
            error = function(e2) NULL)
        })
      if (!is.null(mcs_r)) { print(mcs_r); cat("\n") }
    }
  }
} else {
  cat("=== PARTE 9: MCS nao disponivel (pacote MCS nao instalou) ===\n\n")
}

# ============================================================
# PARTE 10: ANALISE TVP DOS BETAS (todos horizontes)
# ============================================================
cat("=== PARTE 10: Analise TVP dos Betas ===\n\n")

if (file.exists("forecasts/coulombe_betas_2SRR.rda") &&
    file.exists("forecasts/coulombe_betas_ridge.rda")) {

  load("forecasts/coulombe_betas_2SRR.rda")
  load("forecasts/coulombe_betas_ridge.rda")

  for (hi in seq_along(hor)) {
    h <- hor[hi]
    valid_b <- Filter(Negate(is.null), betas_2srr[[hi]])
    if (length(valid_b) < 5) {
      cat(sprintf("  h=%d: %d janelas (insuficiente)\n", h, length(valid_b)))
      next
    }

    cat(sprintf("\n  --- h=%d | %d janelas ---\n", h, length(valid_b)))

    beta_series <- do.call(rbind, lapply(valid_b, function(b) {
      bm <- b$betas
      if (is.array(bm) && length(dim(bm)) == 3) bvec <- bm[1,,dim(bm)[3]]
      else if (is.matrix(bm)) bvec <- bm[nrow(bm),]
      else bvec <- as.numeric(bm)
      c(t = b$t, bvec)
    }))

    beta_df <- as.data.frame(beta_series)
    colnames(beta_df) <- c("t", paste0("beta", 0:(ncol(beta_df) - 2)))
    for (j in seq_len(ncol(beta_df)))
      beta_df[, j] <- as.numeric(as.character(beta_df[, j]))

    # Variancia e CV
    bcols <- grep("^beta", colnames(beta_df), value = TRUE)
    var_b <- sapply(bcols, function(c) var(beta_df[[c]], na.rm = TRUE))
    cv_b  <- sapply(bcols, function(c) {
      m <- mean(beta_df[[c]], na.rm = TRUE)
      s <- sd(beta_df[[c]], na.rm = TRUE)
      ifelse(abs(m) > 1e-10, abs(s/m), NA)
    })

    var_df <- data.frame(param = bcols, var = var_b, cv = cv_b)
    var_df <- var_df[order(-var_df$var), ]

    n_tvp <- sum(var_df$cv > 0.5, na.rm = TRUE)
    cat(sprintf("  TVP genuinos (CV>0.5): %d/%d\n", n_tvp, nrow(var_df)))

    # Top 5
    for (i in seq_len(min(5, nrow(var_df))))
      cat(sprintf("    %-10s var=%.6f cv=%.3f\n",
                  var_df$param[i], var_df$var[i], var_df$cv[i]))

    write.csv(beta_df, sprintf("results/beta_traj_h%02d.csv", h), row.names = FALSE)
    write.csv(var_df, sprintf("results/beta_var_h%02d.csv", h), row.names = FALSE)
  }
}
cat("\n")

# ============================================================
# PARTE 11: SUB-AMOSTRAS (pre/pos COVID)
# ============================================================
cat("=== PARTE 11: Analise por sub-amostra ===\n\n")

covid_date <- as.Date("2020-03-01")

for (h in hor) {
  key <- paste0("h", h)
  if (is.null(coulombe[[key]])) next
  df <- coulombe[[key]]
  df$date_parsed <- as.Date(df$date)
  real <- df$realized

  pre  <- df$date_parsed < covid_date
  post <- df$date_parsed >= covid_date

  if (sum(pre) < 20 || sum(post) < 10) next

  for (period in c("Pre-COVID", "Pos-COVID")) {
    idx <- if (period == "Pre-COVID") pre else post
    r <- real[idx]
    rmse_r <- sqrt(mean((df$fc_ridge[idx] - r)^2, na.rm = TRUE))
    rmse_2 <- sqrt(mean((df$fc_2srr[idx] - r)^2, na.rm = TRUE))

    dm_sub <- tryCatch(
      dm.test(e1 = df$fc_ridge[idx] - r, e2 = df$fc_2srr[idx] - r,
              h = h, alternative = "greater", power = 2),
      error = function(e) NULL)

    cat(sprintf("  h=%2d %s (n=%d): Ridge=%.4f 2SRR=%.4f ratio=%.3f",
                h, period, sum(idx), rmse_r, rmse_2, rmse_2/rmse_r))
    if (!is.null(dm_sub))
      cat(sprintf(" DM p=%.3f", dm_sub$p.value))
    cat("\n")
  }
}
cat("\n")

# ============================================================
# PARTE 12: GRAFICOS
# ============================================================
cat("=== PARTE 12: Graficos ===\n\n")

for (h in hor) {
  key <- paste0("h", h)
  if (is.null(coulombe[[key]])) next
  df <- coulombe[[key]]
  df$date_parsed <- as.Date(df$date)

  # 12A: Forecasts vs Realizado
  p1 <- ggplot(df, aes(x = date_parsed)) +
    geom_line(aes(y = realized, color = "Realizado"), linewidth = 0.5) +
    geom_line(aes(y = fc_ridge, color = "Ridge"), linewidth = 0.4, alpha = 0.7) +
    geom_line(aes(y = fc_2srr, color = "2SRR"), linewidth = 0.4, alpha = 0.7) +
    labs(title = sprintf("Forecasts h=%d", h),
         x = "", y = "Inflacao acumulada", color = "") +
    theme_minimal() + theme(legend.position = "bottom")
  ggsave(sprintf("results/figures/forecast_h%02d.pdf", h), p1,
         width = 10, height = 5)

print(p1)

  # 12B: CSSED
  d_t <- (df$fc_ridge - df$realized)^2 - (df$fc_2srr - df$realized)^2
  df$cssed <- cumsum(d_t)
  p2 <- ggplot(df, aes(x = date_parsed, y = cssed)) +
    geom_line(linewidth = 0.6, color = "steelblue") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    labs(title = sprintf("CSSED h=%d (acima de 0 = 2SRR melhor)", h),
         x = "", y = "CSSED") +
    theme_minimal()
  ggsave(sprintf("results/figures/cssed_h%02d.pdf", h), p2,
         width = 10, height = 4)

print(p2)

  cat(sprintf("  h=%d: graficos salvos\n", h))
}

# 12C: Betas TVP (h=1)
if (file.exists("results/beta_traj_h01.csv")) {
  bt <- read.csv("results/beta_traj_h01.csv")
  if (ncol(bt) > 3) {
    # Seleciona top 4 betas mais variaveis
    bcols <- grep("^beta", colnames(bt), value = TRUE)
    vars <- sapply(bcols, function(c) var(bt[[c]], na.rm = TRUE))
    top4 <- names(sort(vars, decreasing = TRUE))[1:min(4, length(vars))]

    bt_long <- reshape2::melt(bt[, c("t", top4)], id.vars = "t",
                              variable.name = "beta", value.name = "valor")

    p3 <- ggplot(bt_long, aes(x = t, y = valor, color = beta)) +
      geom_line(linewidth = 0.5) +
      facet_wrap(~ beta, scales = "free_y", ncol = 2) +
      labs(title = "Betas TVP (h=1) — Top 4 mais variaveis",
           x = "Observacao", y = "Valor do beta") +
      theme_minimal() + theme(legend.position = "none")
    ggsave("results/figures/betas_tvp_h01.pdf", p3, width = 10, height = 8)
    cat("  Betas TVP h=1: salvo\n")
  }
}
cat("\n")

print(p3)

# ============================================================
# PARTE 13: TABELA LATEX
# ============================================================
cat("=== PARTE 13: Tabela LaTeX ===\n\n")

if (length(all_metrics) > 0) {
  tab_all <- do.call(rbind, all_metrics)
  tab_all <- tab_all[, c("h", "model", "RMSE", "MAE", "rank", "ratio_vs_2SRR")]

  write.csv(tab_all, "results/tabela_completa.csv", row.names = FALSE)

  # LaTeX
  tex <- xtable(tab_all,
                caption = "Forecast Comparison: RMSE and MAE by horizon",
                label = "tab:forecast_comparison",
                digits = c(0, 0, 0, 4, 4, 0, 3))
  print(tex, file = "results/tabela_forecast.tex",
        include.rownames = FALSE,
        booktabs = TRUE,
        caption.placement = "top")
  cat("  Tabela LaTeX: results/tabela_forecast.tex\n")

  # Tabela resumo DM
  if (length(dm_results) > 0) {
    dm_tab <- do.call(rbind, dm_results)
    write.csv(dm_tab, "results/dm_results.csv", row.names = FALSE)

    tex_dm <- xtable(dm_tab,
                     caption = "Diebold-Mariano Test: 2SRR vs All Models",
                     label = "tab:dm_test",
                     digits = c(0, 0, 0, 3, 4, 0))
    print(tex_dm, file = "results/tabela_dm.tex",
          include.rownames = FALSE, booktabs = TRUE,
          caption.placement = "top")
    cat("  Tabela DM LaTeX: results/tabela_dm.tex\n")
  }
}

cat("\n============================================================\n")
cat("  07_validacao_econometrica.R v3.0 --- COMPLETO\n")
cat("\nOutputs gerados:\n")
cat("  results/tabela_completa.csv\n")
cat("  results/tabela_forecast.tex\n")
cat("  results/tabela_dm.tex (se DM rodou)\n")
cat("  results/dm_results.csv\n")
cat("  results/cssed_hXX.csv (por horizonte)\n")
cat("  results/beta_traj_hXX.csv (por horizonte)\n")
cat("  results/beta_var_hXX.csv (por horizonte)\n")
cat("  results/figures/forecast_hXX.pdf\n")
cat("  results/figures/cssed_hXX.pdf\n")
cat("  results/figures/betas_tvp_h01.pdf\n")
