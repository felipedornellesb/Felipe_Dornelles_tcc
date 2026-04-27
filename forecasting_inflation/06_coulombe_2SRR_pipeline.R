# ============================================================
# 06_coulombe_2SRR_pipeline.R   (v6.1)
#
# Script UNICO e AUTONOMO — nao precisa de 06a nem 06b.
# Incorpora tudo que o 06a fazia, com correcoes.
#
# NOVIDADES v6.1 vs v6:
#   [N1] rm_const(): remove colunas constantes antes do prcomp
#        (prcomp(..., scale.=TRUE) falha com var=0)
#        Incorporado do 06a_coulombe_setup.R
#   [N2] Alinhamento de tau com yout.rda do Medeiros:
#        se yout.rda existir, usa nrow(yout) como n_oos;
#        caso contrario, cai no hardcode nwindows=312
#        Incorporado do 06a_coulombe_setup.R
#   [N3] factor() PCA definida antes dos source()
#        Evita "unused argument (n_fac)" no EM_sw.R
#
# BUGS CORRIGIDOS em v6 (mantidos):
#   [BUG1] EM_sw fora do loop (1x total, nao 1x por iteracao)
#   [BUG2] prcomp sem duplicacao por iteracao
#   [BUG3] last e forecast em nivel acumulado correto
#   [BUG4] Checkpoint a cada 50 iteracoes
#
# FLUXO COMPLETO:
#   1. Pacotes
#   2. factor() PCA (antes dos source)
#   3. Funcoes Coulombe
#   4. Filtro outliers OF()
#   5. Base Medeiros
#   6. Alinhamento tau com yout.rda
#   7. Imputacao EM (1x, fora do loop)
#   8. Variavel Y acumulada
#   9. Parametros POOS
#  10. Loop POOS (Ridge + 2SRR)
#  11. Salva resultados
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
# 1. factor() PCA ANTES dos source()  [N3 do 06a]
#
# EM_sw.R chama internamente: factor(X, n_fac=n)
# A base::factor nao aceita n_fac e lanca erro.
# Esta versao PCA precisa estar no global ANTES do source(EM_sw.R).
# ============================================================
factor <- function(X, n_fac) {
  X    <- as.matrix(X)
  Tobs <- nrow(X)
  S    <- (1/Tobs) * t(X) %*% X
  eig  <- eigen(S, symmetric=TRUE)
  lam  <- eig$vectors[, 1:n_fac, drop=FALSE]
  fac  <- X %*% lam
  fit  <- fac %*% t(lam)
  mse  <- mean((X - fit)^2, na.rm=TRUE)
  list(factors=fac, lambda=lam, mse=mse)
}
cat("[OK] factor() PCA definida (compativel com EM_sw)\n")

# ============================================================
# 2. CARREGA FUNCOES DO COULOMBE
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
  source(p, local=FALSE)  # local=FALSE: deposita no global (necessario para EM_sw achar factor)
  cat(sprintf("  [OK] %s\n", f))
}
cat("=== Carregando funcoes Coulombe ===\n")
cs("EM_sw.R"); cs("ICp2.R"); cs("Xgenerators_v190127.R")
cs("dualGRRmdA_v190215.R"); cs("CVGSBHK_v181127.R"); cs("zfun_v190304.R")
cs("factor.R"); cs("TVPRRcosso_v181120.R"); cs("TVPRR_v181111.R")
cs("fastZrot_v181125.R"); cs("CVKFMV_v190214.R")

# Verifica funcoes criticas
for (fn in c("make_reg_matrix", "TVPRR_cosso", "EM_sw")) {
  if (!exists(fn)) stop(sprintf("[CRITICO] %s nao encontrada apos source()", fn))
}
cat("Todas as funcoes carregadas.\n\n")

# ============================================================
# 3. FILTRO DE OUTLIERS (Hugo, Empirical_v2.R)
# ============================================================
OF <- function(pred, y, tol=2, go.to.pred) {
  newx <- pred
  cm   <- (newx - mean(y)) > tol*(max(y) - mean(y))
  cmi  <- (newx - mean(y)) < tol*(min(y) - mean(y))
  newx[cm]  <- go.to.pred[cm]
  newx[cmi] <- go.to.pred[cmi]
  newx
}

# ============================================================
# 4. HELPER: remove colunas constantes  [N1 do 06a]
#
# prcomp(..., scale.=TRUE) falha se qualquer coluna tiver
# variancia zero. Ocorre especialmente em janelas pequenas
# com muitos preditores (116 variaveis FRED-MD).
# ============================================================
rm_const <- function(X) {
  keep <- apply(X, 2, function(col) {
    v <- var(col, na.rm=TRUE)
    !is.na(v) && v > .Machine$double.eps
  })
  if (sum(!keep) > 0)
    cat(sprintf("    [rm_const] %d coluna(s) constante(s) removida(s)\n", sum(!keep)))
  X[, keep, drop=FALSE]
}

# ============================================================
# 5. CARREGA BASE DO MEDEIROS
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
# 6. ALINHAMENTO tau COM yout.rda DO MEDEIROS  [N2 do 06a]
#
# Prioridade 1: usa nrow(yout) de forecasts/yout.rda
#               (garante janela OOS identica ao Medeiros)
# Prioridade 2: fallback para nwindows=312 se nao existir
# ============================================================
cat("\nAlinhamento OOS com o Medeiros:\n")
if (file.exists("forecasts/yout.rda")) {
  load("forecasts/yout.rda")
  n_oos <- nrow(yout)
  tau   <- bigt - n_oos
  cat(sprintf("  yout.rda encontrado -> n_oos=%d (alinhado com Medeiros)\n", n_oos))
} else {
  nwindows <- 312
  n_oos    <- nwindows
  tau      <- bigt - n_oos
  cat(sprintf("  yout.rda nao encontrado -> usando nwindows=%d (hardcode)\n", nwindows))
}
if (tau < 50) stop(sprintf("tau=%d muito pequeno. Verifique os dados.", tau))
cat(sprintf("  tau=%d | OOS: %s -> %s | n_oos=%d\n\n",
            tau, as.character(dates[tau+1]), as.character(dates[bigt]), n_oos))

# ============================================================
# 7. IMPUTACAO EM — UMA VEZ, FORA DO LOOP  [BUG1 v6]
#
# Hugo roda EM_sw uma vez antes do POOS.
# factor() PCA foi definida no passo 1, compativel com EM_sw.
# ============================================================
cat("Imputacao EM Stock & Watson (n=8, it_max=1000)...\n")
X_imp <- tryCatch({
  em_out <- EM_sw(data=as.data.frame(X_raw), n=8, it_max=1000)
  cat(sprintf("  EM_sw OK | NAs restantes: %d\n", sum(is.na(em_out$data))))
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
# 8. VARIAVEL Y ACUMULADA h-PASSOS
#    y_cum[t] = y[t+1] + ... + y[t+h]  (direct multi-step)
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
# 9. PARAMETROS POOS
# ============================================================
nf        <- 8    # fatores PCA no regressor
ly        <- 2    # lags de y
lf        <- 2    # lags dos fatores
lambdavec <- exp(pracma::linspace(-2, 12, n=15))  # grid lambda (Hugo)
silent    <- 1

cat(sprintf("POOS: bigt=%d | tau=%d | n_oos=%d | nf=%d | ly=%d | lf=%d\n",
            bigt, tau, n_oos, nf, ly, lf))
cat(sprintf("In-sample : %s a %s\n", as.character(dates[1]), as.character(dates[tau])))
cat(sprintf("OOS       : %s a %s\n\n", as.character(dates[tau+1]), as.character(dates[bigt])))

dir.create("forecasts",   showWarnings=FALSE)
dir.create("results",     showWarnings=FALSE)
dir.create("checkpoints", showWarnings=FALSE)

# ============================================================
# 10. ARRAYS DE RESULTADO
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
# 11. LOOP POOS PRINCIPAL
# ============================================================
closeAllConnections()  # fix Windows: limite de conexoes

cat("=== INICIANDO LOOP POOS ===\n")
t0_total <- proc.time()

for (hi in seq_along(hor)) {
  h <- hor[hi]
  cat(sprintf("\n========================================\n"))
  cat(sprintf("  Horizonte h = %d meses\n", h))
  cat(sprintf("========================================\n"))
  y_target <- forecast_vars[, hi]

  for (t in tau:(bigt-1)) {
    idx   <- t - tau + 1L
    T_end <- t - h          # information set ate t-h (sem look-ahead)
    if (T_end < (ly + lf + 20L)) next

    # PCA incremental sobre X_imp[1:T_end]  [BUG2 v6]
    # rm_const evita erro de prcomp com colunas constantes  [N1]
    X_is <- rm_const(X_imp[1:T_end, , drop=FALSE])
    pc   <- tryCatch(prcomp(X_is, center=TRUE, scale.=TRUE), error=function(e) NULL)
    if (is.null(pc)) next
    fac  <- pc$x[, 1:min(nf, ncol(pc$x)), drop=FALSE]

    # Variavel dependente e nivel
    y_h <- y_target[1:T_end]
    Y_h <- y_raw[1:T_end]

    first_valid <- which(!is.na(y_h))[1]
    if (is.na(first_valid)) next
    si <- first_valid
    if ((T_end - si + 1L) < (ly + lf + 10L)) next

    y_is <- as.matrix(y_h[si:T_end])
    Y_is <- as.matrix(Y_h[si:T_end])
    f_is <- fac[si:T_end, , drop=FALSE]

    # Matriz de regressao
    # Colunas: [y_cum | Y_{t-1} | Y_{t-2} | f1_{t-1} | f1_{t-2} | ...]
    reg <- tryCatch(
      make_reg_matrix(y=y_is, Y=Y_is, factors=f_is, h=h, ly=ly, lf=lf),
      error=function(e) NULL
    )
    if (is.null(reg) || nrow(reg) < (ly + lf + h + 5L)) next

    # last = ultima linha (obs t): [y_cum_t | regressores_t]
    # last[-1] = vetor de regressores para prever em t+h  [BUG3 v6]
    last <- as.numeric(reg[nrow(reg), ])
    reg  <- reg[1:(nrow(reg)-1L), , drop=FALSE]
    ml   <- max(ly, lf)
    if (nrow(reg) <= ml) next
    reg  <- reg[(ml+1L):nrow(reg), , drop=FALSE]
    reg  <- reg[complete.cases(reg), , drop=FALSE]
    if (nrow(reg) < 15L || ncol(reg) < 2L) next

    xnew <- matrix(last[-1], nrow=1)  # regressores para previsao

    # ----------------------------------------------------------
    # m=1: RIDGE PLANO
    # ----------------------------------------------------------
    CV <- tryCatch(
      cv.glmnet(x=reg[,-1,drop=FALSE], y=reg[,1],
                family="gaussian", alpha=0),
      error=function(e) NULL
    )
    if (is.null(CV)) next

    mdl_r    <- glmnet(x=reg[,-1,drop=FALSE], y=reg[,1],
                       family="gaussian", alpha=0, lambda=CV$lambda.min)
    pred_lin <- as.numeric(predict(mdl_r, newx=xnew))

    fc_ridge[t, hi]  <- pred_lin   # nivel acumulado
    lam_ridge[t, hi] <- CV$lambda.min

    rc <- as.numeric(coef(mdl_r))  # [intercepto, b1, ..., bK]
    betas_ridge[[hi]][[idx]] <- list(
      t=t, date=dates[t], beta0=rc[1], betas=rc[-1]
    )

    # ----------------------------------------------------------
    # m=2: 2SRR (Coulombe — TVPRR_cosso)
    # ----------------------------------------------------------
    aa <- tryCatch(
      TVPRR_cosso(
        y=reg[,1], X=reg[,-1,drop=FALSE],
        lambdavec=lambdavec, sweigths=1, type=2,
        alpha=0.01, silent=silent, kfold=5,
        lambda2=CV$lambda.min, tol=1e-6, maxit=10,
        oosX=as.numeric(xnew)
      ),
      error=function(e) NULL
    )

    if (!is.null(aa)) {
      if (!is.null(aa$fcast)) {
        p_raw <- as.numeric(aa$fcast)
      } else {
        bm  <- aa$grrats$betas_grr
        bl  <- if (is.matrix(bm)) bm[nrow(bm),] else as.numeric(bm)
        p_raw <- as.numeric(bl[1] + sum(bl[-1] * as.numeric(xnew)))
      }

      p_filt <- OF(pred=p_raw, y=reg[,1], go.to.pred=pred_lin)

      fc_2srr[t, hi]   <- p_filt
      lam2_2srr[t, hi] <- if (!is.null(aa$grrats$lambdas)) aa$grrats$lambdas[1] else NA_real_
      betas_2srr[[hi]][[idx]] <- list(
        t=t, date=dates[t],
        betas=aa$grrats$betas_grr  # matriz T x (K+1): trajetoria TVP
      )
    }

    # Progresso a cada 12 obs (~1 ano)
    if (idx %% 12L == 0L) {
      el  <- (proc.time() - t0_total)["elapsed"]
      rem <- el / idx * (n_oos - idx)
      cat(sprintf("  h=%2d | %s | %3d/%d (%3.0f%%) | %.1f min | ~%.0f min restam\n",
                  h, as.character(dates[t]), idx, n_oos,
                  100*idx/n_oos, el/60, rem/60))
    }

    # CHECKPOINT a cada 50 obs  [BUG4 v6]
    if (idx %% 50L == 0L) {
      cp <- list(hi=hi, t=t, idx=idx,
                 fc_ridge=fc_ridge, fc_2srr=fc_2srr,
                 lam_ridge=lam_ridge, lam2_2srr=lam2_2srr,
                 betas_2srr=betas_2srr, betas_ridge=betas_ridge)
      save(cp, file=sprintf("checkpoints/cp_h%d_t%d.rda", h, t))
    }

  }  # fim loop t
  cat(sprintf("  h=%d concluido. Acumulado: %.1f min\n",
              h, (proc.time()-t0_total)["elapsed"]/60))
}  # fim loop hi

cat(sprintf("\nPOOS completo: %.1f min\n\n",
            (proc.time()-t0_total)["elapsed"]/60))

# ============================================================
# 12. SALVA RESULTADOS FINAIS
# ============================================================
save(fc_ridge, fc_2srr, lam_ridge, lam2_2srr,
     file="forecasts/coulombe_forecasts.rda")
save(betas_2srr,  file="forecasts/coulombe_betas_2SRR.rda")
save(betas_ridge, file="forecasts/coulombe_betas_ridge.rda")
cat("RDAs salvos em forecasts/\n")

# CSVs por horizonte
oos_idx <- (tau+1L):bigt
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
  fname <- sprintf("forecasts/coulombe_fc_h%02d.csv", h)
  write.csv(df_out, file=fname, row.names=FALSE)
  cat(sprintf("Salvo: %s | Ridge=%d validos | 2SRR=%d validos\n",
              fname, sum(!is.na(fr)), sum(!is.na(f2))))
}

# CSVs de betas h=1 (para analise comparativa)
hi1 <- which(hor == 1L)

valid_r <- Filter(Negate(is.null), betas_ridge[[hi1]])
if (length(valid_r) > 0L) {
  df_rh1 <- do.call(rbind, lapply(valid_r, function(b)
    c(date=as.character(b$date), t=b$t, beta0=b$beta0,
      setNames(b$betas, paste0("b", seq_along(b$betas))))))
  write.csv(as.data.frame(df_rh1), "results/betas_ridge_h1.csv", row.names=FALSE)
  cat(sprintf("Betas Ridge h=1: %d janelas -> results/betas_ridge_h1.csv\n", nrow(df_rh1)))
}

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

cat("\n=== 06_coulombe_2SRR_pipeline.R v6.1 --- COMPLETO ===\n")
