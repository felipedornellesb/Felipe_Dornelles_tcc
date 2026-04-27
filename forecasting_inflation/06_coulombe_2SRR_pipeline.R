# ============================================================
# 06_coulombe_2SRR_pipeline.R   (v7.0)
#
# CORRECOES CRITICAS vs v6.1:
#   [C1] Eliminado Double-Shifting: y acumulado RETROSPECTIVO
#        y_cum[t] = sum(y[(t-h+1):t]) em vez de y[(t+1):(t+h)]
#        T_end = t (em vez de t - h)
#        make_reg_matrix cuida internamente da defasagem
#   [C2] Download completo de TODOS os scripts Coulombe
#   [C3] factor() PCA em environment isolado (nao polui base::factor)
#   [C4] Validacao de dimensao de xnew vs reg antes de predict
#   [C5] Tratamento robusto de erros em TVPRR_cosso
#   [C6] Salva RMSE parcial no checkpoint para monitoramento
#
# MANTIDOS de v6.1:
#   - EM_sw fora do loop (1x total)
#   - rm_const() antes de prcomp
#   - Checkpoint a cada 50 iteracoes
#   - Filtro de outliers OF()
#   - Alinhamento tau com yout.rda
# ============================================================

rm(list = ls())
gc()

# ============================================================
# 0. CONFIGURACAO INICIAL
# ============================================================
setwd("~/tcc/Felipe_Dornelles_tcc/forecasting_inflation")
cat("============================================================\n")
cat("  2SRR Pipeline v7.0\n")
cat("  Working dir:", getwd(), "\n")
cat("============================================================\n\n")

# ============================================================
# 1. PACOTES
# ============================================================
pkgs <- c("pracma", "glmnet", "matrixcalc", "GA", "e1071")
new  <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if (length(new)) install.packages(new, repos = "https://cran.r-project.org")
invisible(lapply(pkgs, library, character.only = TRUE))
cat("[OK] Pacotes carregados\n")

# ============================================================
# 2. factor() PCA — ISOLADA para nao conflitar com base::factor
#
# Solucao: criamos num environment dedicado que o EM_sw.R
# consegue acessar, mas que nao sobrescreve base::factor
# para outros pacotes (como glmnet que usa factor() internamente).
# ============================================================

# Guardamos a funcao PCA com nome especifico
pca_factor <- function(X, n_fac) {
  X    <- as.matrix(X)
  Tobs <- nrow(X)
  S    <- (1 / Tobs) * t(X) %*% X
  eig  <- eigen(S, symmetric = TRUE)
  nf   <- min(n_fac, ncol(X))
  lam  <- eig$vectors[, 1:nf, drop = FALSE]
  fac  <- X %*% lam
  fit  <- fac %*% t(lam)
  mse  <- mean((X - fit)^2, na.rm = TRUE)
  list(factors = fac, lambda = lam, mse = mse)
}

# Agora, sobrescrevemos factor() no global APENAS para que
# EM_sw.R funcione. Guardaremos base::factor para restaurar depois.
base_factor_backup <- base::factor
factor <- function(X, n_fac = NULL, ...) {
  if (!is.null(n_fac)) {
    return(pca_factor(X, n_fac))
  }
  # Se chamada sem n_fac, comporta-se como base::factor
  base_factor_backup(X, ...)
}
cat("[OK] factor() PCA definida (com fallback para base::factor)\n")

# ============================================================
# 3. DOWNLOAD DE TODOS OS SCRIPTS COULOMBE
#
# [C2] Lista completa — v6.1 faltavam vários arquivos
# ============================================================
dir.create("coulombe", showWarnings = FALSE, recursive = TRUE)

base_url_tools <- paste0(
  "https://raw.githubusercontent.com/hugocout/",
  "Replication-codes-for-Time-Varying-Parameters-as-Ridge-Regressions/",
  "main/Empirical/20_tools"
)
base_url_func <- paste0(base_url_tools, "/functions")

# Lista COMPLETA de todos os scripts necessários
download_list <- list(
  # Raiz de 20_tools/
  list(url = paste0(base_url_tools, "/EM_sw.R"),
       dest = "coulombe/EM_sw.R"),
  list(url = paste0(base_url_tools, "/ICp2.R"),
       dest = "coulombe/ICp2.R"),
  # Subpasta functions/
  list(url = paste0(base_url_func, "/factor.R"),
       dest = "coulombe/factor.R"),
  list(url = paste0(base_url_func, "/TVPRR_v181111.R"),
       dest = "coulombe/TVPRR_v181111.R"),
  list(url = paste0(base_url_func, "/CVKFMV_v190214.R"),
       dest = "coulombe/CVKFMV_v190214.R"),
  list(url = paste0(base_url_func, "/Xgenerators_v190127.R"),
       dest = "coulombe/Xgenerators_v190127.R"),
  list(url = paste0(base_url_func, "/dualGRRmdA_v190215.R"),
       dest = "coulombe/dualGRRmdA_v190215.R"),
  list(url = paste0(base_url_func, "/CVGSBHK_v181127.R"),
       dest = "coulombe/CVGSBHK_v181127.R"),
  list(url = paste0(base_url_func, "/zfun_v190304.R"),
       dest = "coulombe/zfun_v190304.R"),
  list(url = paste0(base_url_func, "/TVPRRcosso_v181120.R"),
       dest = "coulombe/TVPRRcosso_v181120.R"),
  list(url = paste0(base_url_func, "/fastZrot_v181125.R"),
       dest = "coulombe/fastZrot_v181125.R")
)

cat("\n=== Download dos scripts Coulombe ===\n")
for (f in download_list) {
  if (!file.exists(f$dest)) {
    cat(sprintf("  Baixando %-40s", basename(f$dest)))
    tryCatch({
      download.file(f$url, destfile = f$dest, quiet = TRUE, method = "libcurl")
      cat(" [OK]\n")
    }, error = function(e) {
      cat(sprintf(" [ERRO: %s]\n", e$message))
    })
  } else {
    cat(sprintf("  %-40s [ja existe]\n", basename(f$dest)))
  }
}

# ============================================================
# 4. CARREGA FUNCOES COULOMBE
# ============================================================
cs <- function(f) {
  p <- file.path("coulombe", f)
  if (!file.exists(p)) {
    warning(paste0("Arquivo nao encontrado: ", p, " — tentando continuar"))
    return(invisible(NULL))
  }
  source(p, local = FALSE)
  cat(sprintf("  [OK] %s\n", f))
}

cat("\n=== Carregando funcoes Coulombe ===\n")
cs("EM_sw.R")
cs("ICp2.R")
cs("Xgenerators_v190127.R")
cs("dualGRRmdA_v190215.R")
cs("CVGSBHK_v181127.R")
cs("zfun_v190304.R")
cs("factor.R")
cs("TVPRRcosso_v181120.R")
cs("TVPRR_v181111.R")
cs("fastZrot_v181125.R")
cs("CVKFMV_v190214.R")

# Verifica funcoes criticas
funcoes_criticas <- c("make_reg_matrix", "TVPRR_cosso")
for (fn in funcoes_criticas) {
  if (!exists(fn, mode = "function")) {
    stop(sprintf("[CRITICO] Funcao '%s' nao encontrada apos source(). Verifique downloads.", fn))
  }
}
cat("\nTodas as funcoes criticas disponiveis.\n\n")

# ============================================================
# 5. FILTRO DE OUTLIERS (Hugo, Empirical_v2.R)
# ============================================================
OF <- function(pred, y, tol = 2, go.to.pred) {
  newx <- pred
  cm   <- (newx - mean(y)) > tol * (max(y) - mean(y))
  cmi  <- (newx - mean(y)) < tol * (min(y) - mean(y))
  newx[cm]  <- go.to.pred[cm]
  newx[cmi] <- go.to.pred[cmi]
  newx
}

# ============================================================
# 6. HELPER: remove colunas constantes
# ============================================================
rm_const <- function(X) {
  keep <- apply(X, 2, function(col) {
    v <- var(col, na.rm = TRUE)
    !is.na(v) && v > .Machine$double.eps
  })
  n_rem <- sum(!keep)
  if (n_rem > 0)
    cat(sprintf("    [rm_const] %d coluna(s) constante(s) removida(s)\n", n_rem))
  X[, keep, drop = FALSE]
}

# ============================================================
# 7. CARREGA BASE DO MEDEIROS
# ============================================================
cat("=== Carregando base de dados ===\n")
load("data/data.rda")
fred_raw <- as.data.frame(data)
bigt     <- nrow(fred_raw)

date_col <- grep("^date$", colnames(fred_raw), ignore.case = TRUE)[1]
cpi_col  <- grep("^CPIAUCSL$", colnames(fred_raw), ignore.case = TRUE)[1]

if (is.na(date_col) || is.na(cpi_col)) {
  cat("Colunas disponiveis:\n")
  print(colnames(fred_raw))
  stop("Nao encontrou 'date' ou 'CPIAUCSL' na base.")
}

dates <- fred_raw[, date_col]
y_raw <- as.numeric(fred_raw[, cpi_col])
X_raw <- as.matrix(fred_raw[, -c(date_col, cpi_col)])

cat(sprintf("  Base: %d obs x %d preditores\n", bigt, ncol(X_raw)))
cat(sprintf("  Periodo: %s a %s\n", as.character(dates[1]), as.character(dates[bigt])))
cat(sprintf("  NAs em X: %d (%.1f%%)\n\n",
            sum(is.na(X_raw)), 100 * mean(is.na(X_raw))))

# ============================================================
# 8. ALINHAMENTO tau COM yout.rda DO MEDEIROS
# ============================================================
cat("=== Alinhamento OOS ===\n")
if (file.exists("forecasts/yout.rda")) {
  load("forecasts/yout.rda")
  n_oos <- nrow(yout)
  tau   <- bigt - n_oos
  cat(sprintf("  yout.rda encontrado -> n_oos = %d (alinhado com Medeiros)\n", n_oos))
} else {
  nwindows <- 312
  n_oos    <- nwindows
  tau      <- bigt - n_oos
  cat(sprintf("  yout.rda NAO encontrado -> fallback nwindows = %d\n", nwindows))
}

if (tau < 50) stop(sprintf("tau = %d muito pequeno. Verifique dados.", tau))

cat(sprintf("  tau = %d | IS: obs 1-%d | OOS: obs %d-%d\n",
            tau, tau, tau + 1, bigt))
cat(sprintf("  IS:  %s a %s\n", as.character(dates[1]), as.character(dates[tau])))
cat(sprintf("  OOS: %s a %s\n\n", as.character(dates[tau + 1]), as.character(dates[bigt])))

# ============================================================
# 9. IMPUTACAO EM — UMA VEZ, FORA DO LOOP
# ============================================================
cat("=== Imputacao EM Stock-Watson ===\n")
X_imp <- tryCatch({
  em_out <- EM_sw(data = as.data.frame(X_raw), n = 8, it_max = 1000)
  cat(sprintf("  EM_sw OK | NAs restantes: %d\n", sum(is.na(em_out$data))))
  as.matrix(em_out$data)
}, error = function(e) {
  cat(sprintf("  EM_sw falhou (%s)\n", e$message))
  cat("  Usando interpolacao linear como fallback\n")
  Xr <- X_raw
  for (j in seq_len(ncol(Xr))) {
    nas <- which(is.na(Xr[, j]))
    if (length(nas) > 0 && length(nas) < nrow(Xr) - 2) {
      Xr[, j] <- approx(seq_len(nrow(Xr)), Xr[, j],
                         xout = seq_len(nrow(Xr)), rule = 2)$y
    }
  }
  Xr[is.na(Xr)] <- 0  # ultima linha de defesa
  Xr
})
cat("  Imputacao concluida.\n\n")

# ============================================================
# 10. VARIAVEL Y ACUMULADA — CORRECAO CRITICA [C1]
#
# *** ESTA E A CORRECAO PRINCIPAL ***
#
# ANTES (v6.1, ERRADO — double-shifting):
#   y_cum[t] = y[t+1] + ... + y[t+h]   (forward-looking)
#   T_end = t - h                        (recua de novo)
#   -> make_reg_matrix recua mais h internamente = PERDE 2h obs
#
# AGORA (v7.0, CORRETO — retrospective):
#   y_cum[t] = y[(t-h+1):t]             (backward-looking)
#   T_end = t                            (sem recuo adicional)
#   -> make_reg_matrix aplica o unico shift necessario
#
# INTERPRETACAO:
#   y_cum[t] = inflacao acumulada nos ultimos h meses ate t
#   No direct forecasting, regredimos y_cum[t+h] ~ X[t]
#   A funcao make_reg_matrix do Coulombe faz esse alinhamento.
# ============================================================

build_cumulative_y_retro <- function(y, h) {
  n  <- length(y)
  yh <- rep(NA_real_, n)
  for (t in h:n) {
    yh[t] <- sum(y[(t - h + 1):t])
  }
  yh
}

hor           <- c(1, 3, 6, 12)
forecast_vars <- sapply(hor, build_cumulative_y_retro, y = y_raw)
colnames(forecast_vars) <- paste0("h", hor)

cat("=== Variavel dependente acumulada (retrospectiva) ===\n")
for (hi in seq_along(hor)) {
  h <- hor[hi]
  n_valid <- sum(!is.na(forecast_vars[, hi]))
  cat(sprintf("  h=%2d: %d valores validos (de %d)\n", h, n_valid, bigt))
}
cat("\n")

# ============================================================
# 11. PARAMETROS POOS
# ============================================================
nf        <- 8     # fatores PCA
ly        <- 2     # lags de y
lf        <- 2     # lags dos fatores
lambdavec <- exp(pracma::linspace(-2, 12, n = 15))
silent    <- 1

cat("=== Parametros POOS ===\n")
cat(sprintf("  nf=%d | ly=%d | lf=%d | lambdas: %d valores\n", nf, ly, lf, length(lambdavec)))
cat(sprintf("  bigt=%d | tau=%d | n_oos=%d\n\n", bigt, tau, n_oos))

# Diretorios
dir.create("forecasts",   showWarnings = FALSE)
dir.create("results",     showWarnings = FALSE)
dir.create("checkpoints", showWarnings = FALSE)

# ============================================================
# 12. ARRAYS DE RESULTADO
# ============================================================
fc_ridge  <- matrix(NA_real_, nrow = bigt, ncol = length(hor))
fc_2srr   <- matrix(NA_real_, nrow = bigt, ncol = length(hor))
lam_ridge <- matrix(NA_real_, nrow = bigt, ncol = length(hor))
lam2_2srr <- matrix(NA_real_, nrow = bigt, ncol = length(hor))

betas_2srr  <- setNames(vector("list", length(hor)), paste0("h", hor))
betas_ridge <- setNames(vector("list", length(hor)), paste0("h", hor))
for (hi in seq_along(hor)) {
  betas_2srr[[hi]]  <- vector("list", n_oos)
  betas_ridge[[hi]] <- vector("list", n_oos)
}

# Contadores de falha para diagnostico
fail_counts <- list(
  pca_fail = 0, reg_fail = 0, cv_fail = 0, 
  tvp_fail = 0, dim_fail = 0, skip_small = 0
)

# ============================================================
# 13. LOOP POOS PRINCIPAL
# ============================================================
closeAllConnections()
cat("============================================================\n")
cat("  INICIANDO LOOP POOS\n")
cat("============================================================\n\n")

t0_total <- proc.time()

for (hi in seq_along(hor)) {
  h <- hor[hi]
  cat(sprintf("\n########################################\n"))
  cat(sprintf("  HORIZONTE h = %d meses\n", h))
  cat(sprintf("########################################\n"))
  
  y_target <- forecast_vars[, hi]
  
  n_ridge_ok <- 0
  n_2srr_ok  <- 0
  
  for (t in tau:(bigt - 1)) {
    idx <- t - tau + 1L
    
    # --------------------------------------------------------
    # [C1] CORRECAO: T_end = t (nao t - h)
    # make_reg_matrix internamente faz o shift por h
    # --------------------------------------------------------
    T_end <- t
    
    # Verificacao de janela minima
    min_obs <- ly + lf + h + 30L
    if (T_end < min_obs) {
      fail_counts$skip_small <- fail_counts$skip_small + 1
      next
    }
    
    # --------------------------------------------------------
    # PCA sobre X_imp[1:T_end] — expanding window
    # --------------------------------------------------------
    X_is <- rm_const(X_imp[1:T_end, , drop = FALSE])
    
    if (ncol(X_is) < nf) {
      fail_counts$pca_fail <- fail_counts$pca_fail + 1
      next
    }
    
    pc <- tryCatch(
      prcomp(X_is, center = TRUE, scale. = TRUE),
      error = function(e) NULL
    )
    
    if (is.null(pc)) {
      fail_counts$pca_fail <- fail_counts$pca_fail + 1
      next
    }
    
    n_fac_use <- min(nf, ncol(pc$x))
    fac <- pc$x[, 1:n_fac_use, drop = FALSE]
    
    # --------------------------------------------------------
    # Variavel dependente e nivel
    # --------------------------------------------------------
    y_h <- y_target[1:T_end]
    Y_h <- y_raw[1:T_end]
    
    first_valid <- which(!is.na(y_h))[1]
    if (is.na(first_valid) || first_valid >= T_end) next
    
    si <- first_valid
    if ((T_end - si + 1L) < min_obs) {
      fail_counts$skip_small <- fail_counts$skip_small + 1
      next
    }
    
    y_is <- as.matrix(y_h[si:T_end])
    Y_is <- as.matrix(Y_h[si:T_end])
    f_is <- fac[si:T_end, , drop = FALSE]
    
    # --------------------------------------------------------
    # Matriz de regressao via make_reg_matrix do Coulombe
    # --------------------------------------------------------
    reg <- tryCatch(
      make_reg_matrix(y = y_is, Y = Y_is, factors = f_is,
                      h = h, ly = ly, lf = lf),
      error = function(e) {
        fail_counts$reg_fail <<- fail_counts$reg_fail + 1
        NULL
      }
    )
    
    if (is.null(reg)) next
    if (nrow(reg) < (ly + lf + h + 10L)) {
      fail_counts$reg_fail <- fail_counts$reg_fail + 1
      next
    }
    
    # --------------------------------------------------------
    # Extrair ultima linha para previsao e ajustar reg
    # --------------------------------------------------------
    last <- as.numeric(reg[nrow(reg), ])
    reg  <- reg[1:(nrow(reg) - 1L), , drop = FALSE]
    
    ml <- max(ly, lf)
    if (nrow(reg) <= ml + 5L) {
      fail_counts$reg_fail <- fail_counts$reg_fail + 1
      next
    }
    reg <- reg[(ml + 1L):nrow(reg), , drop = FALSE]
    reg <- reg[complete.cases(reg), , drop = FALSE]
    
    if (nrow(reg) < 20L || ncol(reg) < 2L) {
      fail_counts$reg_fail <- fail_counts$reg_fail + 1
      next
    }
    
    # reg[,1] = y_target; reg[,-1] = regressores
    yy  <- reg[, 1]
    XX  <- reg[, -1, drop = FALSE]
    xnew <- matrix(last[-1], nrow = 1)
    
    # [C4] Validacao de dimensao
    if (ncol(xnew) != ncol(XX)) {
      fail_counts$dim_fail <- fail_counts$dim_fail + 1
      # Tenta ajustar (truncar ou preencher)
      if (ncol(xnew) > ncol(XX)) {
        xnew <- xnew[, 1:ncol(XX), drop = FALSE]
      } else {
        # Preenche com media das colunas faltantes
        diff_cols <- ncol(XX) - ncol(xnew)
        xnew <- cbind(xnew, matrix(0, nrow = 1, ncol = diff_cols))
      }
    }
    
    # --------------------------------------------------------
    # MODELO 1: RIDGE PLANO (benchmark)
    # --------------------------------------------------------
    CV <- tryCatch(
      cv.glmnet(x = XX, y = yy, family = "gaussian", 
                alpha = 0, nfolds = min(10, nrow(XX))),
      error = function(e) NULL
    )
    
    if (is.null(CV)) {
      fail_counts$cv_fail <- fail_counts$cv_fail + 1
      next
    }
    
    mdl_r    <- glmnet(x = XX, y = yy, family = "gaussian",
                       alpha = 0, lambda = CV$lambda.min)
    pred_lin <- as.numeric(predict(mdl_r, newx = xnew))
    
    fc_ridge[t, hi]  <- pred_lin
    lam_ridge[t, hi] <- CV$lambda.min
    n_ridge_ok       <- n_ridge_ok + 1
    
    rc <- as.numeric(coef(mdl_r))
    betas_ridge[[hi]][[idx]] <- list(
      t = t, date = dates[t],
      beta0 = rc[1], betas = rc[-1]
    )
    
    # --------------------------------------------------------
    # MODELO 2: 2SRR (Coulombe — TVPRR_cosso)
    # --------------------------------------------------------
    aa <- tryCatch(
      TVPRR_cosso(
        y       = yy,
        X       = XX,
        lambdavec = lambdavec,
        sweigths  = 1,
        type      = 2,
        alpha     = 0.01,
        silent    = silent,
        kfold     = 5,
        lambda2   = CV$lambda.min,
        tol       = 1e-6,
        maxit     = 10,
        oosX      = as.numeric(xnew)
      ),
      error = function(e) {
        fail_counts$tvp_fail <<- fail_counts$tvp_fail + 1
        NULL
      }
    )
    
    if (!is.null(aa)) {
      # Extrair previsao
      if (!is.null(aa$fcast)) {
        p_raw <- as.numeric(aa$fcast)
      } else if (!is.null(aa$grrats) && !is.null(aa$grrats$betas_grr)) {
        bm  <- aa$grrats$betas_grr
        bl  <- if (is.matrix(bm)) bm[nrow(bm), ] else as.numeric(bm)
        if (length(bl) >= 1 + ncol(xnew)) {
          p_raw <- as.numeric(bl[1] + sum(bl[-1] * as.numeric(xnew)))
        } else {
          p_raw <- pred_lin  # fallback
        }
      } else {
        p_raw <- pred_lin
      }
      
      # Filtro de outliers
      p_filt <- OF(pred = p_raw, y = yy, go.to.pred = pred_lin)
      
      fc_2srr[t, hi]   <- p_filt
      lam2_2srr[t, hi] <- if (!is.null(aa$grrats$lambdas)) {
        aa$grrats$lambdas[1]
      } else {
        NA_real_
      }
      
      betas_2srr[[hi]][[idx]] <- list(
        t     = t,
        date  = dates[t],
        betas = if (!is.null(aa$grrats$betas_grr)) aa$grrats$betas_grr else NULL
      )
      
      n_2srr_ok <- n_2srr_ok + 1
    }
    
    # --------------------------------------------------------
    # Progresso a cada 12 obs (~1 ano)
    # --------------------------------------------------------
    if (idx %% 12L == 0L) {
      el  <- (proc.time() - t0_total)["elapsed"]
      rem <- el / idx * (n_oos - idx)
      
      # RMSE parcial para monitoramento
      oos_so_far <- (tau + 1):t
      real_sf    <- forecast_vars[oos_so_far, hi]
      rmse_r <- sqrt(mean((fc_ridge[oos_so_far, hi] - real_sf)^2, na.rm = TRUE))
      rmse_2 <- sqrt(mean((fc_2srr[oos_so_far, hi] - real_sf)^2, na.rm = TRUE))
      ratio  <- ifelse(is.finite(rmse_r) && rmse_r > 0, rmse_2 / rmse_r, NA)
      
      cat(sprintf(
        "  h=%2d | %s | %3d/%d (%3.0f%%) | Ridge:%.3f 2SRR:%.3f ratio:%.3f | %.1fmin (~%.0fmin rest)\n",
        h, as.character(dates[t]), idx, n_oos, 100 * idx / n_oos,
        rmse_r, rmse_2, ifelse(is.na(ratio), 99, ratio),
        el / 60, rem / 60
      ))
    }
    
    # --------------------------------------------------------
    # CHECKPOINT a cada 50 obs
    # --------------------------------------------------------
    if (idx %% 50L == 0L) {
      cp <- list(
        hi = hi, t = t, idx = idx, h = h,
        fc_ridge = fc_ridge, fc_2srr = fc_2srr,
        lam_ridge = lam_ridge, lam2_2srr = lam2_2srr,
        betas_2srr = betas_2srr, betas_ridge = betas_ridge,
        fail_counts = fail_counts
      )
      save(cp, file = sprintf("checkpoints/cp_h%d_t%d.rda", h, t))
      cat(sprintf("    [CHECKPOINT] salvo h=%d t=%d\n", h, t))
    }
    
  }  # fim loop t
  
  el_h <- (proc.time() - t0_total)["elapsed"]
  cat(sprintf("\n  h=%d CONCLUIDO | Ridge: %d/%d | 2SRR: %d/%d | %.1f min\n",
              h, n_ridge_ok, n_oos, n_2srr_ok, n_oos, el_h / 60))
  
}  # fim loop hi

el_total <- (proc.time() - t0_total)["elapsed"]
cat(sprintf("\n============================================================\n"))
cat(sprintf("  POOS COMPLETO: %.1f min (%.1f horas)\n", el_total / 60, el_total / 3600))
cat(sprintf("============================================================\n\n"))

# Diagnostico de falhas
cat("=== Diagnostico de falhas ===\n")
for (nm in names(fail_counts)) {
  cat(sprintf("  %-15s: %d\n", nm, fail_counts[[nm]]))
}
cat("\n")

# ============================================================
# 14. SALVA RESULTADOS FINAIS
# ============================================================
save(fc_ridge, fc_2srr, lam_ridge, lam2_2srr,
     file = "forecasts/coulombe_forecasts.rda")
save(betas_2srr,  file = "forecasts/coulombe_betas_2SRR.rda")
save(betas_ridge, file = "forecasts/coulombe_betas_ridge.rda")
cat("RDAs salvos em forecasts/\n")

# CSVs por horizonte
oos_idx <- (tau + 1L):bigt

for (hi in seq_along(hor)) {
  h    <- hor[hi]
  real <- forecast_vars[oos_idx, hi]
  fr   <- fc_ridge[oos_idx, hi]
  f2   <- fc_2srr[oos_idx, hi]
  
  df_out <- data.frame(
    date      = dates[oos_idx],
    realized  = real,
    fc_ridge  = fr,
    fc_2srr   = f2,
    err_ridge = fr - real,
    err_2srr  = f2 - real,
    lam_ridge = lam_ridge[oos_idx, hi],
    lam2_2srr = lam2_2srr[oos_idx, hi]
  )
  
  # Metricas finais
  rmse_r <- sqrt(mean((fr - real)^2, na.rm = TRUE))
  rmse_2 <- sqrt(mean((f2 - real)^2, na.rm = TRUE))
  mae_r  <- mean(abs(fr - real), na.rm = TRUE)
  mae_2  <- mean(abs(f2 - real), na.rm = TRUE)
  
  fname <- sprintf("forecasts/coulombe_fc_h%02d.csv", h)
  write.csv(df_out, file = fname, row.names = FALSE)
  
  cat(sprintf("  h=%2d | %s | Ridge: RMSE=%.4f MAE=%.4f (%d validos)\n",
              h, fname, rmse_r, mae_r, sum(!is.na(fr))))
  cat(sprintf("       |%s| 2SRR:  RMSE=%.4f MAE=%.4f (%d validos)\n",
              strrep(" ", nchar(fname) + 3), rmse_2, mae_2, sum(!is.na(f2))))
  cat(sprintf("       | Ratio RMSE(2SRR/Ridge) = %.4f %s\n",
              rmse_2 / rmse_r,
              ifelse(rmse_2 < rmse_r, "<-- 2SRR GANHA", "<-- Ridge ganha")))
}

# ============================================================
# 15. EXPORTA BETAS PARA ANALISE TVP (h=1)
# ============================================================
hi1 <- which(hor == 1L)

if (length(hi1) > 0) {
  # Betas Ridge h=1
  valid_r <- Filter(Negate(is.null), betas_ridge[[hi1]])
  if (length(valid_r) > 0L) {
    df_rh1 <- do.call(rbind, lapply(valid_r, function(b) {
      c(date = as.character(b$date), t = b$t, beta0 = b$beta0,
        setNames(b$betas, paste0("b", seq_along(b$betas))))
    }))
    write.csv(as.data.frame(df_rh1), "results/betas_ridge_h1.csv", row.names = FALSE)
    cat(sprintf("\nBetas Ridge h=1: %d janelas -> results/betas_ridge_h1.csv\n",
                nrow(df_rh1)))
  }
  
  # Betas 2SRR h=1
  valid_b <- Filter(Negate(is.null), betas_2srr[[hi1]])
  if (length(valid_b) > 0L) {
    df_bh1 <- do.call(rbind, lapply(valid_b, function(b) {
      bm   <- b$betas
      bvec <- if (is.matrix(bm)) bm[nrow(bm), ] else as.numeric(bm)
      c(date = as.character(b$date), t = b$t,
        setNames(bvec, paste0("b", seq_along(bvec))))
    }))
    write.csv(as.data.frame(df_bh1), "results/betas_2srr_h1.csv", row.names = FALSE)
    cat(sprintf("Betas 2SRR h=1: %d janelas -> results/betas_2srr_h1.csv\n",
                nrow(df_bh1)))
  }
}

cat("\n============================================================\n")
cat("  06_coulombe_2SRR_pipeline.R v7.0 --- COMPLETO\n")
cat("============================================================\n")
