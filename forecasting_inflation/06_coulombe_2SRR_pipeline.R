# ============================================================
# 06_coulombe_2SRR_pipeline.R   (v6 — correcoes criticas)
#
# BUGS CORRIGIDOS vs v5.1:
#
#   [BUG1] EM_sw dentro do loop t — ELIMINADO.
#          Cada chamada EM_sw sobre ~470x116 com it_max=1000
#          demorava 2-5 seg. Multiplicado por 312*4 iteracoes
#          = 40+ min so de EM_sw. Solucao: EM_sw UMA VEZ
#          antes do loop sobre X_raw completo, gerando X_imp.
#          Perda de informacao negligenciavel vs ganho de tempo.
#
#   [BUG2] PCA dentro do loop t (prcomp duplicado) — ELIMINADO.
#          Rodava prcomp duas vezes por iteracao (X_is e X_full).
#          Solucao: PCA incremental — uma PCA por t sobre
#          X_imp[1:T_end], rapida pois X_imp ja esta limpo.
#
#   [BUG3] last[-1] errado como vetor de regressores.
#          make_reg_matrix retorna matriz cujas colunas sao
#          [y_h | y_{t-1} | y_{t-2} | f1_{t-1} | f1_{t-2} | ...].
#          A coluna 1 e a variavel dependente (y_h), nao
#          intercepto. Portanto last[-1] = regressores corretos
#          e last[1] = y_h_{t} (ultimo valor acumulado in-sample).
#          A subtracao fc = predict - last[1] esta correta
#          SOMENTE para h=1. Para h>1, last[1] e y_acumulado,
#          nao y_nivel. Corrigido: salva forecast em nivel
#          absoluto e converte para erro na analise posterior.
#          Padrao Hugo: o CSVs ja tem 'realized' para comparar.
#
#   [BUG4] Checkpoint ausente. Se o R travasse apos 3h, perdia
#          tudo. Adicionado: save incremental a cada 50 obs.
#
# PARAMETROS CANONICOS (Coulombe 2024, Table 16-17):
#   ly=2, lf=2, nf=8 (fatores PCA), lambdavec=exp(linspace(-2,12,15))
#   alpha=0.01, kfold=5, tol=1e-6, maxit=10
#   nwindows=312 (identico ao hugocout para base mensal USA)
# ============================================================

rm(list = ls())
setwd("~/tcc/Felipe_Dornelles_tcc/forecasting_inflation")
cat("Working dir:", getwd(), "\n")

# ============================================================
# 0. PACOTES
# ============================================================
pkgs <- c("pracma", "glmnet", "timeSeries", "matrixcalc", "GA", "e1071")
new  <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if (length(new)) install.packages(new, repos="http://cran.us.r-project.org")
invisible(lapply(pkgs, library, character.only=TRUE))

# ============================================================
# 1. CARREGA FUNCOES DO COULOMBE
# ============================================================
base_raw <- paste0(
  "https://raw.githubusercontent.com/hugocout/",
  "Replication-codes-for-Time-Varying-Parameters-as-Ridge-Regressions/",
  "main/Empirical/20_tools"
)
extras <- list(
  list(url=paste0(base_raw,"/EM_sw.R"),                    dest="coulombe/EM_sw.R"),
  list(url=paste0(base_raw,"/ICp2.R"),                     dest="coulombe/ICp2.R"),
  list(url=paste0(base_raw,"/functions/factor.R"),         dest="coulombe/factor.R"),
  list(url=paste0(base_raw,"/functions/TVPRR_v181111.R"),  dest="coulombe/TVPRR_v181111.R"),
  list(url=paste0(base_raw,"/functions/CVKFMV_v190214.R"), dest="coulombe/CVKFMV_v190214.R")
)
for (f in extras) {
  if (!file.exists(f$dest)) {
    cat(sprintf("  Baixando %s ...", basename(f$dest)))
    tryCatch({
      download.file(f$url, destfile=f$dest, quiet=TRUE, method="libcurl")
      cat(" OK\n")
    }, error=function(e) cat(sprintf(" ERRO: %s\n", e$message)))
  }
}

cs <- function(f) {
  p <- file.path("coulombe", f)
  if (!file.exists(p)) stop(paste0("Nao encontrado: ", p))
  source(p); cat(sprintf("  [OK] %s\n", f))
}
cat("=== Carregando funcoes Coulombe ===\n")
cs("EM_sw.R"); cs("ICp2.R"); cs("Xgenerators_v190127.R")
cs("dualGRRmdA_v190215.R"); cs("CVGSBHK_v181127.R"); cs("zfun_v190304.R")
cs("factor.R"); cs("TVPRRcosso_v181120.R"); cs("TVPRR_v181111.R")
cs("fastZrot_v181125.R"); cs("CVKFMV_v190214.R")
cat("Todas as funcoes carregadas.\n\n")

# ============================================================
# 2. FILTRO DE OUTLIERS (Hugo, Empirical_v2.R)
# ============================================================
OF <- function(pred, y, tol=2, go.to.pred) {
  newx <- pred
  newx[(newx - mean(y)) >  tol*(max(y) - mean(y))] <- go.to.pred[(newx - mean(y)) >  tol*(max(y) - mean(y))]
  newx[(newx - mean(y)) <  tol*(min(y) - mean(y))] <- go.to.pred[(newx - mean(y)) <  tol*(min(y) - mean(y))]
  newx
}

# ============================================================
# 3. CARREGA BASE DO MEDEIROS
# ============================================================
load("data/data.rda")
fred_raw <- as.data.frame(data)
bigt     <- nrow(fred_raw)
date_col <- grep("^date$",     colnames(fred_raw), ignore.case=TRUE)[1]
cpi_col  <- grep("^CPIAUCSL$", colnames(fred_raw), ignore.case=TRUE)[1]
if (is.na(date_col) || is.na(cpi_col)) {
  print(colnames(fred_raw))
  stop("Nao encontrou 'date' ou 'CPIAUCSL'.")
}
dates <- fred_raw[, date_col]
y_raw <- fred_raw[, cpi_col]
X_raw <- as.matrix(fred_raw[, -c(date_col, cpi_col)])
cat(sprintf("Base: %d obs x %d preditores | %s a %s\n",
            bigt, ncol(X_raw), as.character(dates[1]), as.character(dates[bigt])))

# ============================================================
# 4. IMPUTACAO EM — UMA VEZ, FORA DO LOOP  [FIX BUG1]
#    Hugo roda EM_sw uma vez sobre toda a base antes do POOS.
#    nf=8: numero de fatores para imputacao (independente do
#    nf usado na PCA do regressor).
# ============================================================
cat("\nImputacao EM Stock & Watson (nf=8, it_max=1000)...\n")
X_imp <- tryCatch({
  em_out <- EM_sw(data=as.data.frame(X_raw), n=8, it_max=1000)
  cat("  EM_sw OK\n")
  as.matrix(em_out$data)
}, error=function(e) {
  cat(sprintf("  EM_sw falhou (%s) -> interpolacao linear\n", e$message))
  Xr <- X_raw
  for (j in seq_len(ncol(Xr))) {
    nas <- which(is.na(Xr[,j]))
    if (length(nas) > 0 && length(nas) < nrow(Xr)-2)
      Xr[,j] <- approx(seq_len(nrow(Xr)), Xr[,j], xout=seq_len(nrow(Xr)))$y
  }
  Xr
})
cat("Imputacao concluida.\n\n")

# ============================================================
# 5. VARIAVEL Y ACUMULADA h-PASSOS (direct multi-step forecast)
#    y_cum[t,h] = y[t+1] + y[t+2] + ... + y[t+h]
#    Equivalente ao que Hugo faz com newQ_targets.csv
# ============================================================
build_cumulative_y <- function(y, h) {
  n <- length(y); yh <- rep(NA_real_, n)
  for (t in seq_len(n-h)) yh[t] <- sum(y[(t+1):(t+h)])
  yh
}
hor           <- c(1, 3, 6, 12)
forecast_vars <- sapply(hor, build_cumulative_y, y=y_raw)
colnames(forecast_vars) <- paste0("h", hor)

# ============================================================
# 6. PARAMETROS POOS
# ============================================================
nwindows  <- 312        # janelas OOS
tau       <- bigt - nwindows  # inicio da avaliacao OOS
n_oos     <- nwindows
nf        <- 8          # fatores PCA no regressor (= Hugo mod=2 usa 2, mas base
                        #   mensal tem mais variacao; use 8 como FRED-MD standard;
                        #   pode testar nf=2 para replicar Hugo mais fielmente)
ly        <- 2          # lags de y no regressor
lf        <- 2          # lags dos fatores no regressor
lambdavec <- exp(pracma::linspace(-2, 12, n=15))  # grid lambda (Hugo)
silent    <- 1

cat(sprintf("POOS: bigt=%d | tau=%d | n_oos=%d | nf=%d | ly=%d | lf=%d\n",
            bigt, tau, n_oos, nf, ly, lf))
cat(sprintf("Periodo in-sample: %s a %s\n",
            as.character(dates[1]), as.character(dates[tau])))
cat(sprintf("Periodo OOS:       %s a %s\n\n",
            as.character(dates[tau+1]), as.character(dates[bigt])))

dir.create("forecasts",     showWarnings=FALSE)
dir.create("results",       showWarnings=FALSE)
dir.create("checkpoints",   showWarnings=FALSE)

# ============================================================
# 7. ARRAYS DE RESULTADO
# ============================================================
fc_ridge  <- matrix(NA_real_, nrow=bigt, ncol=length(hor))
fc_2srr   <- matrix(NA_real_, nrow=bigt, ncol=length(hor))
lam_ridge <- matrix(NA_real_, nrow=bigt, ncol=length(hor))
lam2_2srr <- matrix(NA_real_, nrow=bigt, ncol=length(hor))

betas_2srr  <- setNames(vector("list", length(hor)), paste0("h", hor))
betas_ridge <- setNames(vector("list", length(hor)), paste0("h", hor))
for (hi in seq_along(hor)) {
  betas_2srr[[hi]]  <- vector("list", n_oos)
  betas_ridge[[hi]] <- vector("list", n_oos)
}

# ============================================================
# 8. LOOP POOS PRINCIPAL
# ============================================================
closeAllConnections()  # fix Windows: limite de 128 conexoes

cat("=== INICIANDO LOOP POOS ===\n")
t0_total <- proc.time()

for (hi in seq_along(hor)) {
  h <- hor[hi]
  cat(sprintf("\n========================================\n")
  cat(sprintf("  Horizonte h = %d meses\n", h))
  cat(sprintf("========================================\n"))
  y_target <- forecast_vars[, hi]  # vetor: y_cum para este h

  for (t in tau:(bigt-1)) {
    idx   <- t - tau + 1L
    T_end <- t - h   # information set termina em t-h (sem look-ahead)
    if (T_end < (ly + lf + 20L)) next

    # ----------------------------------------------------------
    # PCA INCREMENTAL sobre X_imp[1:T_end]  [FIX BUG2]
    # Uma PCA por t. X_imp ja esta imputado (sem NAs).
    # Muito mais rapido que EM_sw dentro do loop.
    # ----------------------------------------------------------
    X_is  <- X_imp[1:T_end, , drop=FALSE]
    pc    <- prcomp(X_is, center=TRUE, scale.=TRUE)
    fac   <- pc$x[, 1:min(nf, ncol(pc$x)), drop=FALSE]

    # Variavel dependente acumulada
    y_h <- y_target[1:T_end]
    Y_h <- y_raw[1:T_end]    # nivel de y para lags

    # Remove NAs iniciais de y_h (primeiros h obs sao NA por construcao)
    first_valid <- which(!is.na(y_h))[1]
    if (is.na(first_valid)) next
    si <- first_valid
    if ((T_end - si + 1L) < (ly + lf + 10L)) next

    y_is <- as.matrix(y_h[si:T_end])
    Y_is <- as.matrix(Y_h[si:T_end])
    f_is <- fac[si:T_end, , drop=FALSE]

    # ----------------------------------------------------------
    # CONSTROI MATRIZ DE REGRESSAO
    # make_reg_matrix retorna: [y_cum | lags_Y | lags_fatores]
    # coluna 1 = variavel dependente (y_cum_{t})
    # colunas 2:K = regressores (lags de y e fatores)
    # ultima linha = observacao t (para previsao)
    # ----------------------------------------------------------
    reg <- tryCatch(
      make_reg_matrix(y=y_is, Y=Y_is, factors=f_is, h=h, ly=ly, lf=lf),
      error=function(e) NULL
    )
    if (is.null(reg) || nrow(reg) < (ly + lf + h + 5L)) next

    # [FIX BUG3] last = ultima linha de reg (observacao t)
    # Contem: [y_cum_t | Y_{t-1} | Y_{t-2} | f1_{t-1} | f1_{t-2} | ...]
    # last[1]  = y_cum acumulado (NAO eh intercepto)
    # last[-1] = vetor de regressores para prever em t+h
    last <- as.numeric(reg[nrow(reg), ])

    # Remove a ultima linha (que seria o target t+h, nao disponivel)
    # e remove lags iniciais incompletos
    reg <- reg[1:(nrow(reg)-1L), , drop=FALSE]
    ml  <- max(ly, lf)
    if (nrow(reg) <= ml) next
    reg <- reg[(ml+1L):nrow(reg), , drop=FALSE]
    reg <- reg[complete.cases(reg), , drop=FALSE]
    if (nrow(reg) < 15L || ncol(reg) < 2L) next

    # Vetor de regressores para previsao (sem coluna dependente)
    xnew <- matrix(last[-1], nrow=1)

    # ----------------------------------------------------------
    # m=1: RIDGE PLANO
    # ----------------------------------------------------------
    CV <- tryCatch(
      cv.glmnet(x=reg[, -1, drop=FALSE], y=reg[, 1],
                family="gaussian", alpha=0),
      error=function(e) NULL
    )
    if (is.null(CV)) next

    mdl_r    <- glmnet(x=reg[,-1,drop=FALSE], y=reg[,1],
                       family="gaussian", alpha=0, lambda=CV$lambda.min)
    pred_lin <- as.numeric(predict(mdl_r, newx=xnew))

    # Forecast em nivel acumulado (sem subtrair last[1])
    # A conversao para variacao eh feita na analise (07_compare.R)
    fc_ridge[t, hi]  <- pred_lin
    lam_ridge[t, hi] <- CV$lambda.min

    # Salva betas Ridge (intercepto + coefs)
    rc <- as.numeric(coef(mdl_r))  # [intercepto, b1, b2, ...]
    betas_ridge[[hi]][[idx]] <- list(
      t     = t,
      date  = dates[t],
      beta0 = rc[1],
      betas = rc[-1]
    )

    # ----------------------------------------------------------
    # m=2: 2SRR (Coulombe — TVPRR_cosso)
    # ----------------------------------------------------------
    aa <- tryCatch(
      TVPRR_cosso(
        y         = reg[, 1],
        X         = reg[, -1, drop=FALSE],
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
      error=function(e) NULL
    )

    if (!is.null(aa)) {
      # Extrai previsao do objeto aa
      if (!is.null(aa$fcast)) {
        p_raw <- as.numeric(aa$fcast)
      } else {
        bm  <- aa$grrats$betas_grr
        bl  <- if (is.matrix(bm)) bm[nrow(bm),] else as.numeric(bm)
        # bl tem dimensao K+1 (intercepto incluso) ou K?
        # TVPRR_cosso inclui intercepto: bl[1]=intercept, bl[-1]=betas
        p_raw <- as.numeric(bl[1] + sum(bl[-1] * as.numeric(xnew)))
      }

      # Filtro de outliers do Hugo
      p_filt <- OF(pred=p_raw, y=reg[,1], go.to.pred=pred_lin)

      fc_2srr[t, hi]   <- p_filt
      lam2_2srr[t, hi] <- if (!is.null(aa$grrats$lambdas)) aa$grrats$lambdas[1] else NA_real_

      betas_2srr[[hi]][[idx]] <- list(
        t     = t,
        date  = dates[t],
        betas = aa$grrats$betas_grr  # matriz T x (K+1): trajetoria TVP
      )
    }

    # Progresso a cada 12 obs (~1 ano)
    if (idx %% 12L == 0L) {
      el <- (proc.time() - t0_total)["elapsed"]
      remaining <- el / idx * (n_oos - idx)
      cat(sprintf("  h=%2d | %s | t=%d/%d (%3.0f%%) | %.1f min | restam ~%.0f min\n",
                  h, as.character(dates[t]), t, bigt-1,
                  100*idx/n_oos, el/60, remaining/60))
    }

    # [FIX BUG4] CHECKPOINT: salva estado a cada 50 iteracoes
    # Se o R travar, reinicia a partir do ultimo checkpoint
    if (idx %% 50L == 0L) {
      cp <- list(hi=hi, t=t, idx=idx,
                 fc_ridge=fc_ridge, fc_2srr=fc_2srr,
                 lam_ridge=lam_ridge, lam2_2srr=lam2_2srr,
                 betas_2srr=betas_2srr, betas_ridge=betas_ridge)
      save(cp, file=sprintf("checkpoints/cp_h%d_t%d.rda", h, t))
    }

  }  # fim loop t
  cat(sprintf("  h=%d concluido. %.1f min acumulados\n",
              h, (proc.time()-t0_total)["elapsed"]/60))
}  # fim loop hi

cat(sprintf("\nPOOS completo: %.1f min\n\n",
            (proc.time()-t0_total)["elapsed"]/60))

# ============================================================
# 9. SALVA RESULTADOS FINAIS
# ============================================================
save(fc_ridge, fc_2srr, lam_ridge, lam2_2srr,
     file="forecasts/coulombe_forecasts.rda")
save(betas_2srr,  file="forecasts/coulombe_betas_2SRR.rda")
save(betas_ridge, file="forecasts/coulombe_betas_ridge.rda")
cat("RDAs salvos em forecasts/\n")

# --- CSVs de forecasts por horizonte ---
oos_idx <- (tau+1L):bigt
for (hi in seq_along(hor)) {
  h    <- hor[hi]
  real <- forecast_vars[oos_idx, hi]
  fr   <- fc_ridge[oos_idx, hi]
  f2   <- fc_2srr[oos_idx, hi]

  # Converte forecasts de nivel acumulado para erro de previsao
  # e.u. = forecast - realized
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
  fname <- sprintf("forecasts/coulombe_fc_h%02d.csv", h)
  write.csv(df_out, file=fname, row.names=FALSE)
  cat(sprintf("Salvo: %s | Ridge=%d validos | 2SRR=%d validos\n",
              fname, sum(!is.na(fr)), sum(!is.na(f2))))
}

# --- CSVs de betas h=1 (para 07_compare.R) ---
hi1 <- which(hor == 1L)

# Betas Ridge h=1
valid_r <- Filter(Negate(is.null), betas_ridge[[hi1]])
if (length(valid_r) > 0L) {
  df_rh1 <- do.call(rbind, lapply(valid_r, function(b) {
    c(date=as.character(b$date), t=b$t, beta0=b$beta0,
      setNames(b$betas, paste0("b", seq_along(b$betas))))
  }))
  write.csv(as.data.frame(df_rh1), "results/betas_ridge_h1.csv", row.names=FALSE)
  cat(sprintf("Betas Ridge h=1: %d janelas -> results/betas_ridge_h1.csv\n", nrow(df_rh1)))
}

# Betas 2SRR h=1 (ultima linha da trajetoria TVP por janela)
valid_b <- Filter(Negate(is.null), betas_2srr[[hi1]])
if (length(valid_b) > 0L) {
  df_bh1 <- do.call(rbind, lapply(valid_b, function(b) {
    bm   <- b$betas
    bvec <- if (is.matrix(bm)) bm[nrow(bm),] else as.numeric(bm)
    c(date=as.character(b$date), t=b$t,
      setNames(bvec, paste0("b", seq_along(bvec))))
  }))
  write.csv(as.data.frame(df_bh1), "results/betas_2srr_h1.csv", row.names=FALSE)
  cat(sprintf("Betas 2SRR h=1: %d janelas -> results/betas_2srr_h1.csv\n", nrow(df_bh1)))
}

cat("\n=== 06_coulombe_2SRR_pipeline.R v6 --- COMPLETO ===\n")
