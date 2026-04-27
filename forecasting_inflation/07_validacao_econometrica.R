# ============================================================
# 07_validacao_econometrica.R
#
# Testes estatisticos profissionais para validar o 2SRR:
#   1. Diebold-Mariano (DM) — pairwise
#   2. Model Confidence Set (MCS) — Hansen et al. (2011)
#   3. Mincer-Zarnowitz — eficiencia informacional
#   4. Giacomini-White (2006) — conditional predictive ability
#   5. Analise TVP dos betas
#   6. Fluctuation Test (Giacomini-Rossi, 2010)
# ============================================================

rm(list = ls())
setwd("~/tcc/Felipe_Dornelles_tcc/forecasting_inflation")

# Pacotes
pkgs <- c("forecast", "MCS", "lmtest", "sandwich", "ggplot2",
          "reshape2", "gridExtra", "car")
new  <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if (length(new)) install.packages(new, repos = "https://cran.r-project.org")
invisible(lapply(pkgs, library, character.only = TRUE))

# ============================================================
# 1. CARREGA PREVISOES DO 2SRR/RIDGE
# ============================================================
hor <- c(1, 3, 6, 12)

cat("============================================================\n")
cat("  VALIDACAO ECONOMETRICA DO 2SRR\n")
cat("============================================================\n\n")

results_summary <- list()

for (h in hor) {
  fname <- sprintf("forecasts/coulombe_fc_h%02d.csv", h)
  if (!file.exists(fname)) {
    cat(sprintf("[AVISO] %s nao encontrado, pulando h=%d\n", fname, h))
    next
  }
  
  df <- read.csv(fname, stringsAsFactors = FALSE)
  
  # Remove NAs
  valid <- complete.cases(df$realized, df$fc_ridge, df$fc_2srr)
  df <- df[valid, ]
  
  if (nrow(df) < 20) {
    cat(sprintf("[AVISO] h=%d: apenas %d obs validas, insuficiente\n", h, nrow(df)))
    next
  }
  
  realized <- df$realized
  fc_ridge <- df$fc_ridge
  fc_2srr  <- df$fc_2srr
  err_ridge <- df$err_ridge
  err_2srr  <- df$err_2srr
  n_obs <- length(realized)
  
  cat(sprintf("\n########################################\n"))
  cat(sprintf("  HORIZONTE h = %d | n = %d observacoes\n", h, n_obs))
  cat(sprintf("########################################\n\n"))
  
  # ----------------------------------------------------------
  # METRICAS BASICAS
  # ----------------------------------------------------------
  rmse_r <- sqrt(mean(err_ridge^2))
  rmse_2 <- sqrt(mean(err_2srr^2))
  mae_r  <- mean(abs(err_ridge))
  mae_2  <- mean(abs(err_2srr))
  mape_r <- mean(abs(err_ridge / realized)) * 100
  mape_2 <- mean(abs(err_2srr / realized)) * 100
  
  cat("--- Metricas de Erro ---\n")
  cat(sprintf("  %-8s  RMSE    MAE     MAPE\n", ""))
  cat(sprintf("  %-8s  %.4f  %.4f  %.2f%%\n", "Ridge", rmse_r, mae_r, mape_r))
  cat(sprintf("  %-8s  %.4f  %.4f  %.2f%%\n", "2SRR", rmse_2, mae_2, mape_2))
  cat(sprintf("  Ratio RMSE(2SRR/Ridge) = %.4f\n", rmse_2 / rmse_r))
  cat(sprintf("  Ratio MAE(2SRR/Ridge)  = %.4f\n\n", mae_2 / mae_r))
  
  # ----------------------------------------------------------
  # 1.A TESTE DE DIEBOLD-MARIANO
  # H0: igual precisao preditiva
  # H1 (one-sided): 2SRR mais preciso que Ridge
  # ----------------------------------------------------------
  cat("--- Teste Diebold-Mariano ---\n")
  
  # Usando MSE (power=2)
  dm_mse <- tryCatch(
    dm.test(e1 = err_ridge, e2 = err_2srr, h = h,
            alternative = "greater", power = 2),
    error = function(e) NULL
  )
  
  # Usando MAE (power=1)
  dm_mae <- tryCatch(
    dm.test(e1 = err_ridge, e2 = err_2srr, h = h,
            alternative = "greater", power = 1),
    error = function(e) NULL
  )
  
  if (!is.null(dm_mse)) {
    cat(sprintf("  DM (MSE): stat=%.3f, p-valor=%.4f %s\n",
                dm_mse$statistic, dm_mse$p.value,
                ifelse(dm_mse$p.value < 0.05, "*** SIGNIFICATIVO",
                       ifelse(dm_mse$p.value < 0.10, "* marginalmente sig.", ""))))
  }
  if (!is.null(dm_mae)) {
    cat(sprintf("  DM (MAE): stat=%.3f, p-valor=%.4f %s\n",
                dm_mae$statistic, dm_mae$p.value,
                ifelse(dm_mae$p.value < 0.05, "*** SIGNIFICATIVO",
                       ifelse(dm_mae$p.value < 0.10, "* marginalmente sig.", ""))))
  }
  cat("\n")
  
  # ----------------------------------------------------------
  # 1.B REGRESSAO DE MINCER-ZARNOWITZ
  # realizado_t = alpha + beta * previsao_t + eps_t
  # H0: alpha=0, beta=1 (previsao eficiente)
  # ----------------------------------------------------------
  cat("--- Teste Mincer-Zarnowitz ---\n")
  
  for (model_name in c("Ridge", "2SRR")) {
    fc_use <- if (model_name == "Ridge") fc_ridge else fc_2srr
    
    mz_reg <- lm(realized ~ fc_use)
    mz_sum <- summary(mz_reg)
    
    alpha_hat <- coef(mz_reg)[1]
    beta_hat  <- coef(mz_reg)[2]
    r2        <- mz_sum$r.squared
    
    # Teste F conjunto: H0: alpha=0 E beta=1
    # Usando linearHypothesis do pacote car
    ftest <- tryCatch(
      linearHypothesis(mz_reg, c("(Intercept) = 0", "fc_use = 1")),
      error = function(e) NULL
    )
    
    cat(sprintf("  %s: alpha=%.4f, beta=%.4f, R2=%.4f\n",
                model_name, alpha_hat, beta_hat, r2))
    
    if (!is.null(ftest)) {
      f_pval <- ftest$`Pr(>F)`[2]
      cat(sprintf("    H0(alpha=0,beta=1): F=%.3f, p=%.4f %s\n",
                  ftest$F[2], f_pval,
                  ifelse(f_pval > 0.05, "[NAO rejeita -> EFICIENTE]",
                         "[Rejeita -> viés detectado]")))
    }
    
    # Teste com Newey-West HAC para autocorrelacao nos erros
    nw_se <- tryCatch({
      coeftest(mz_reg, vcov = NeweyWest(mz_reg, lag = h))
    }, error = function(e) NULL)
    
    if (!is.null(nw_se)) {
      cat(sprintf("    Newey-West: alpha p=%.4f | beta p=%.4f (H0:coef=0)\n",
                  nw_se[1, 4], nw_se[2, 4]))
    }
  }
  cat("\n")
  
  # ----------------------------------------------------------
  # 1.C FORECAST ENCOMPASSING TEST (Fair-Shiller, 1990)
  # realizado_t = alpha + beta1*fc_ridge_t + beta2*fc_2srr_t + eps
  # Se beta2 sig. e beta1 nao sig. -> 2SRR encompassa Ridge
  # ----------------------------------------------------------
  cat("--- Forecast Encompassing (Fair-Shiller) ---\n")
  
  enc_reg <- lm(realized ~ fc_ridge + fc_2srr)
  enc_sum <- summary(enc_reg)
  
  # HAC standard errors
  enc_nw <- tryCatch(
    coeftest(enc_reg, vcov = NeweyWest(enc_reg, lag = h)),
    error = function(e) coeftest(enc_reg)
  )
  
  cat(sprintf("  realized = %.4f + %.4f*Ridge + %.4f*2SRR\n",
              coef(enc_reg)[1], coef(enc_reg)[2], coef(enc_reg)[3]))
  cat(sprintf("  HAC p-valores: Ridge=%.4f | 2SRR=%.4f\n",
              enc_nw[2, 4], enc_nw[3, 4]))
  
  if (enc_nw[3, 4] < 0.05 && enc_nw[2, 4] > 0.05) {
    cat("  >>> 2SRR ENCOMPASSA Ridge (2SRR contem info nao captada pelo Ridge)\n")
  } else if (enc_nw[2, 4] < 0.05 && enc_nw[3, 4] > 0.05) {
    cat("  >>> Ridge encompassa 2SRR\n")
  } else if (enc_nw[2, 4] < 0.05 && enc_nw[3, 4] < 0.05) {
    cat("  >>> Ambos contribuem — combinacao e ideal\n")
  } else {
    cat("  >>> Nenhum contribui significativamente (amostra pequena?)\n")
  }
  cat("\n")
  
  # ----------------------------------------------------------
  # 1.D CUMULATIVE SUM OF SQUARED ERROR DIFFERENCES (CSSED)
  # Giacomini-Rossi style: visualiza quando 2SRR ganha/perde
  # ----------------------------------------------------------
  cat("--- CSSED (Relative Performance ao longo do tempo) ---\n")
  
  d_t <- err_ridge^2 - err_2srr^2   # positivo = 2SRR melhor
  cssed <- cumsum(d_t)
  
  # Salva para grafico
  df_cssed <- data.frame(
    date  = df$date,
    cssed = cssed,
    d_t   = d_t
  )
  write.csv(df_cssed, sprintf("results/cssed_h%02d.csv", h), row.names = FALSE)
  
  # Estatisticas do CSSED
  pct_2srr_wins <- mean(d_t > 0) * 100
  cat(sprintf("  CSSED final = %.4f (%s)\n",
              tail(cssed, 1),
              ifelse(tail(cssed, 1) > 0, "2SRR acumula vantagem",
                     "Ridge acumula vantagem")))
  cat(sprintf("  2SRR ganha em %.1f%% dos periodos\n\n", pct_2srr_wins))
  
  # ----------------------------------------------------------
  # Armazena resultados para tabela final
  # ----------------------------------------------------------
  results_summary[[paste0("h", h)]] <- data.frame(
    h = h,
    n_obs = n_obs,
    RMSE_Ridge = rmse_r,
    RMSE_2SRR  = rmse_2,
    Ratio_RMSE = rmse_2 / rmse_r,
    MAE_Ridge  = mae_r,
    MAE_2SRR   = mae_2,
    Ratio_MAE  = mae_2 / mae_r,
    DM_MSE_stat = ifelse(!is.null(dm_mse), dm_mse$statistic, NA),
    DM_MSE_pval = ifelse(!is.null(dm_mse), dm_mse$p.value, NA),
    DM_MAE_stat = ifelse(!is.null(dm_mae), dm_mae$statistic, NA),
    DM_MAE_pval = ifelse(!is.null(dm_mae), dm_mae$p.value, NA),
    pct_2SRR_wins = pct_2srr_wins
  )
}

# ============================================================
# 2. MODEL CONFIDENCE SET (MCS) — Incluindo modelos Medeiros
# ============================================================
cat("\n============================================================\n")
cat("  MODEL CONFIDENCE SET (Hansen et al., 2011)\n")
cat("============================================================\n\n")

# Tenta carregar previsoes dos modelos Medeiros
medeiros_models <- list()
forecast_files <- list.files("forecasts", pattern = "^for_.*\\.rda$", full.names = TRUE)

if (length(forecast_files) > 0) {
  for (ff in forecast_files) {
    model_name <- gsub("^for_|\\.rda$", "", basename(ff))
    tryCatch({
      env <- new.env()
      load(ff, envir = env)
      obj_names <- ls(env)
      if (length(obj_names) > 0) {
        medeiros_models[[model_name]] <- get(obj_names[1], envir = env)
      }
    }, error = function(e) {
      cat(sprintf("  [AVISO] Falha ao carregar %s: %s\n", ff, e$message))
    })
  }
  cat(sprintf("  Modelos Medeiros carregados: %d\n", length(medeiros_models)))
  cat(sprintf("  Nomes: %s\n\n", paste(names(medeiros_models), collapse = ", ")))
} else {
  cat("  Nenhum forecast Medeiros encontrado em forecasts/\n")
  cat("  MCS sera rodado apenas com Ridge vs 2SRR\n\n")
}

# MCS por horizonte
for (h in hor) {
  hi <- which(hor == h)
  fname <- sprintf("forecasts/coulombe_fc_h%02d.csv", h)
  if (!file.exists(fname)) next
  
  df <- read.csv(fname, stringsAsFactors = FALSE)
  valid <- complete.cases(df$realized, df$fc_ridge, df$fc_2srr)
  df <- df[valid, ]
  if (nrow(df) < 30) next
  
  # Monta matriz de losses
  loss_list <- list(
    Ridge = (df$fc_ridge - df$realized)^2,
    SRR2  = (df$fc_2srr  - df$realized)^2
  )
  
  # Tenta adicionar modelos Medeiros
  # (Precisa alinhar temporalmente — assume mesma janela OOS)
  for (mn in names(medeiros_models)) {
    mm <- medeiros_models[[mn]]
    if (is.matrix(mm) || is.data.frame(mm)) {
      if (ncol(mm) >= hi && nrow(mm) >= nrow(df)) {
        fc_med <- mm[1:nrow(df), hi]
        if (!all(is.na(fc_med))) {
          loss_list[[mn]] <- (fc_med - df$realized)^2
        }
      }
    }
  }
  
  # Converte para matriz
  loss_mat <- do.call(cbind, loss_list)
  loss_mat <- loss_mat[complete.cases(loss_mat), , drop = FALSE]
  
  if (ncol(loss_mat) >= 2 && nrow(loss_mat) >= 30) {
    cat(sprintf("  MCS h=%d: %d modelos, %d obs\n", h, ncol(loss_mat), nrow(loss_mat)))
    
    mcs_result <- tryCatch({
      MCSprocedure(
        Loss = loss_mat,
        alpha = 0.10,
        B = 5000,
        statistic = "Tmax"
      )
    }, error = function(e) {
      cat(sprintf("    MCS falhou: %s\n", e$message))
      # Tenta com Tr (menos restritivo)
      tryCatch(
        MCSprocedure(Loss = loss_mat, alpha = 0.10, B = 5000, statistic = "TR"),
        error = function(e2) NULL
      )
    })
    
    if (!is.null(mcs_result)) {
      cat("  Resultado MCS (alpha=0.10):\n")
      print(mcs_result)
      cat("\n")
    }
  }
}

# ============================================================
# 3. TABELA RESUMO FINAL
# ============================================================
cat("\n============================================================\n")
cat("  TABELA RESUMO FINAL\n")
cat("============================================================\n\n")

if (length(results_summary) > 0) {
  tab_final <- do.call(rbind, results_summary)
  print(tab_final, digits = 4)
  write.csv(tab_final, "results/tabela_resumo_final.csv", row.names = FALSE)
  cat("\nTabela salva em results/tabela_resumo_final.csv\n")
}

# ============================================================
# 4. ANALISE DOS PARAMETROS VARIANTES NO TEMPO (TVP)
# ============================================================
cat("\n============================================================\n")
cat("  ANALISE TVP DOS BETAS\n")
cat("============================================================\n\n")

if (file.exists("forecasts/coulombe_betas_2SRR.rda") &&
    file.exists("forecasts/coulombe_betas_ridge.rda")) {
  
  load("forecasts/coulombe_betas_2SRR.rda")
  load("forecasts/coulombe_betas_ridge.rda")
  
  hi1 <- 1  # h=1
  
  # Extrai betas 2SRR validos
  valid_b <- Filter(Negate(is.null), betas_2srr[[hi1]])
  
  if (length(valid_b) > 5) {
    cat(sprintf("  %d janelas de betas 2SRR disponveis para h=1\n\n", length(valid_b)))
    
    # Extrai o ultimo vetor de betas de cada janela
    beta_series <- do.call(rbind, lapply(valid_b, function(b) {
      bm <- b$betas
      if (is.matrix(bm)) {
        bvec <- bm[nrow(bm), ]
      } else {
        bvec <- as.numeric(bm)
      }
      c(t = b$t, bvec)
    }))
    
    beta_df <- as.data.frame(beta_series)
    colnames(beta_df) <- c("t", paste0("beta", 0:(ncol(beta_df) - 2)))
    
    # Converte para numerico
    for (j in seq_len(ncol(beta_df))) {
      beta_df[, j] <- as.numeric(as.character(beta_df[, j]))
    }
    
    # Calcula variancia de cada beta ao longo do tempo
    cat("  Variancia temporal dos betas (TVP diagnostic):\n")
    cat("  Betas com alta variancia = genuinamente variantes no tempo\n")
    cat("  Betas com var ~0 = efetivamente constantes (Ridge bastava)\n\n")
    
    beta_cols <- grep("^beta", colnames(beta_df), value = TRUE)
    var_betas <- sapply(beta_cols, function(col) var(beta_df[[col]], na.rm = TRUE))
    
    var_df <- data.frame(
      parametro = names(var_betas),
      variancia = var_betas,
      cv = sapply(beta_cols, function(col) {
        m <- mean(beta_df[[col]], na.rm = TRUE)
        s <- sd(beta_df[[col]], na.rm = TRUE)
        ifelse(abs(m) > 1e-10, abs(s / m), NA)
      })
    )
    var_df <- var_df[order(-var_df$variancia), ]
    
    cat(sprintf("  %-12s  %12s  %12s\n", "Parametro", "Variancia", "CV"))
    cat(paste(rep("-", 40), collapse = ""), "\n")
    for (i in seq_len(min(nrow(var_df), 15))) {
      cat(sprintf("  %-12s  %12.6f  %12.4f\n",
                  var_df$parametro[i], var_df$variancia[i], var_df$cv[i]))
    }
    
    # Quantos betas sao "genuinamente TVP"?
    # Criterio: CV > 0.5 (coeficiente de variacao substancial)
    n_tvp <- sum(var_df$cv > 0.5, na.rm = TRUE)
    n_const <- sum(var_df$cv <= 0.5, na.rm = TRUE)
    cat(sprintf("\n  Betas genuinamente TVP (CV > 0.5): %d/%d\n", n_tvp, nrow(var_df)))
    cat(sprintf("  Betas efetivamente constantes:     %d/%d\n\n", n_const, nrow(var_df)))
    
    # Salva para graficos externos
    write.csv(beta_df, "results/beta_trajectories_h1.csv", row.names = FALSE)
    write.csv(var_df, "results/beta
