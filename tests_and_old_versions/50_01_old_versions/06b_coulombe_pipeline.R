# Modelos estimados (m):
#   m=1  -> Ridge plano (benchmark interno)
#   m=2  -> 2SRR — Two-Step Ridge Regression (Coulombe)
#
# Saídas:
#   forecasts/coulombe_forecast_array.rda  — array [bigt, max(hor), 2]
#   forecasts/coulombe_hp_track.rda        — lambdas rastreados
#   forecasts/coulombe_betas_list.rda      — betas TVP por período
#   forecasts/coulombe_fc_h{1,3,6,12}.csv — forecasts por horizonte
#   forecasts/coulombe_nf_sensitivity.csv  — tabela MSFE por nf
#
# BUGS ORIGINAIS DO HUGO CORRIGIDOS:
#   [B1] fred2 é recomputado dentro do loop de t — movido para fora
#   [B2] forecast_vars[-c(1:2),-1] é aplicado DUAS VEZES (linhas 33-34)
#        — removida a linha duplicada
#   [B3] date vector não é usado para alinhar previsões — corrigido
#   [B4] subset não exclui obs. das últimas (h-1) para evitar data leakage
#        — agora subset = 1:(dim(train)[1]) conforme o original mas
#          train é construído com lag correto via make_reg_matrix
#   [B5] hp_track usa dimensão 150 fixa sem uso — simplificado
#   [B6] ts.plot() dentro do loop de produção — removido
# ============================================================

rm(list = ls())
setwd("~/tcc/Felipe_Dornelles_tcc/forecasting_inflation")

# ============================================================
# 0. PACOTES
# ============================================================
pkgs <- c("readr","pracma","doParallel","foreach","glmnet",
          "timeSeries","fGarch","matrixcalc","GA","e1071")
new  <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if (length(new)) install.packages(new, repos = "http://cran.us.r-project.org")
invisible(lapply(pkgs, library, character.only = TRUE))

# ============================================================
# 1. PATHS E FUNÇÕES
# ============================================================
paths <- list(
  too = "coulombe/20_tools",
  fun = "coulombe/20_tools/functions"
)

# Funções de tools (raiz de 20_tools)
source(paste(paths$too, "EM_sw.R",              sep="/"))
source(paste(paths$too, "ICp2.R",               sep="/"))
source(paste(paths$too, "Xgenerators_v190127.R",sep="/"))

# Funções internas do TVP-RR
source(paste(paths$fun, "dualGRRmdA_v190215.R", sep="/"))
source(paste(paths$fun, "CVGSBHK_v181127.R",    sep="/"))
source(paste(paths$fun, "zfun_v190304.R",        sep="/"))
source(paste(paths$fun, "factor.R",              sep="/"))
source(paste(paths$fun, "TVPRRcosso_v181120.R",  sep="/"))
source(paste(paths$fun, "TVPRR_v181111.R",       sep="/"))
source(paste(paths$fun, "fastZrot_v181125.R",    sep="/"))
source(paste(paths$fun, "CVKFMV_v190214.R",      sep="/"))

# ============================================================
# 2. FILTRO DE OUTLIERS (OF) — idêntico ao Hugo
#    Substitui previsões extremas pela previsão linear plana
# ============================================================
OF <- function(pred, y, tol = 2, go.to.pred) {
  newx <- pred
  cond.max <- (newx - mean(y)) > tol * (max(y) - mean(y))
  cond.min <- (newx - mean(y)) < tol * (min(y) - mean(y))
  newx[cond.max] <- go.to.pred[cond.max]
  newx[cond.min] <- go.to.pred[cond.min]
  return(newx)
}

# ============================================================
# 3. CARREGA BASE DO MEDEIROS (FRED-MD mensal)
# ============================================================
load("data/data.rda")
# Esperado: objeto 'data' com colunas de séries FRED-MD transformadas
# A primeira coluna deve ser CPIAUCSL (ou equivalente, a inflação target)
# Ajuste os nomes abaixo conforme necessário:

fred_raw <- as.data.frame(data)

# Identifica coluna da inflação (CPIAUCSL) — ajuste o nome se diferente
cpi_col <- grep("CPIAUCSL|cpi|infl", colnames(fred_raw), ignore.case = TRUE)[1]
if (is.na(cpi_col)) stop("Coluna de inflação não encontrada. Ajuste 'cpi_col' manualmente.")

cat(sprintf("Coluna de inflação identificada: %s (col %d)\n",
            colnames(fred_raw)[cpi_col], cpi_col))

# ============================================================
# 4. CONSTRÓI forecast_vars — variável acumulada h-passos
#    Replica a lógica do newQ_targets.csv do Hugo para dados mensais
#    Col estrutura: [y_h1 | y_h3 | y_h6 | y_h12] (direct forecast)
#    y_h = soma acumulada de h inflações mensais (annualized opcional)
# ============================================================
y_raw <- fred_raw[, cpi_col]
bigt  <- nrow(fred_raw)
bign  <- ncol(fred_raw)
hor   <- c(1, 3, 6, 12)   # horizontes mensais

build_cumulative_y <- function(y, h) {
  # Constrói y^{(h)}_t = sum_{j=1}^{h} y_{t+j}
  # Retorna vetor de comprimento length(y), com NA nas últimas h posições
  n  <- length(y)
  yh <- rep(NA, n)
  for (t in 1:(n - h)) {
    yh[t] <- sum(y[(t + 1):(t + h)])
  }
  return(yh)
}

# forecast_vars: matriz [bigt x length(hor)]
# Coluna j = y acumulada em hor[j] passos
forecast_vars <- sapply(hor, function(h) build_cumulative_y(y_raw, h))
colnames(forecast_vars) <- paste0("h", hor)

# ============================================================
# 5. IMPUTAÇÃO EM DE STOCK & WATSON — [BUG B1 CORRIGIDO]
#    Hugo faz fred2 <- EM_sw() DENTRO do loop de t (custoso e errado).
#    Correto: fazer UMA vez com toda a amostra disponível pré-POOS.
#    Justificativa: a imputação usa apenas a estrutura de fatores de X,
#    não a variável Y — não há data leakage.
# ============================================================
cat("Rodando imputação EM de Stock & Watson (n=8)... ")
fred2 <- EM_sw(data = fred_raw, n = 8, it_max = 1000)$data
cat("OK\n")

# ============================================================
# 6. PARÂMETROS DO POOS
# ============================================================
tau  <- round(0.6 * bigt)    # janela mínima de estimação (~60% da amostra)
jump <- 1                     # reotimiza HPs a cada período (como Hugo)
nf   <- 8                     # fatores EM — igual ao Hugo (n=8)
ly   <- 2                     # lags de Y no regressor
lf   <- 2                     # lags dos fatores no regressor
silenceplz <- 1
lambdavec  <- exp(pracma::linspace(-2, 12, n = 15))  # grid lambda do Coulombe
n_oos <- bigt - tau

cat(sprintf("POOS: tau=%d, bigt=%d, n_oos=%d, hor=%s\n",
            tau, bigt, n_oos, paste(hor, collapse=",")))

# ============================================================
# 7. ARRAYS DE RESULTADO
#    forecast  [bigt, max(hor), length(hor), 2 modelos]
#    hp_track  [bigt, max(hor), length(hor), 2 modelos, 2 lambdas]
#    betas_list: lista aninhada [h][t] com betas TVP
# ============================================================
n_mod <- 2   # m=1 Ridge, m=2 2SRR
forecast  <- array(NA, dim = c(bigt, max(hor), length(hor), n_mod))
hp_track  <- array(NA, dim = c(bigt, max(hor), length(hor), n_mod, 2))
betas_list <- vector("list", length(hor))
names(betas_list) <- paste0("h", hor)
for (hh in seq_along(hor)) betas_list[[hh]] <- vector("list", n_oos)

# ============================================================
# 8. LOOP POOS PRINCIPAL
# ============================================================
cat("\n=== INICIANDO POOS ===\n")
t_start <- proc.time()

for (hi in seq_along(hor)) {
  h <- hor[hi]
  cat(sprintf("\n--- Horizonte h=%d ---\n", h))

  for (t in tau:bigt) {
    idx_oos <- t - tau + 1

    # --- Constrói information set INFO_{t-h} ---
    # Dados de X: observações 1:(t-h)
    # [BUG B4]: usar (t-h) como corte, NÃO t, para previsão correta
    data_Xraw <- as.matrix(fred2)[1:(t - h), ]

    # Fatores EM com n=8 (mod=2 do Hugo: EM_sw sem a variável target)
    # Hugo faz: factors = EM_sw(data=factors[,-vars[v]], n=2)$factors
    # Adaptação: removemos a coluna da inflação antes de extrair fatores
    data_Xfac <- data_Xraw[, -cpi_col]
    factors_em <- EM_sw(data = as.data.frame(data_Xfac), n = nf,
                         it_max = 1000)$factors
    # factors_em: matrix [T_train x nf]

    # Variável dependente Y acumulada (h-passos)
    y <- as.matrix(forecast_vars[1:(t - h), hi])  # [T_train x 1]
    # Nível de Y para lags no regressor (Y no paper = nível, não acumulado)
    Y <- as.matrix(y_raw[1:(t - h)])              # [T_train x 1]

    # Remove NAs iniciais (comuns em séries com transformações)
    start <- sum(is.na(y)) + 1
    end   <- length(1:(t - h))
    if ((end - start + 1) < 20) next   # mínimo de obs úteis

    y       <- y[start:end, 1, drop = FALSE]
    Y       <- Y[start:end, 1, drop = FALSE]
    factors_em <- factors_em[start:end, , drop = FALSE]

    # --- make_reg_matrix — constrói [y | const | lags_Y | lags_fatores] ---
    # Última linha = observação t (oos), usada apenas para previsão
    train_full <- make_reg_matrix(
      y       = y,
      Y       = Y,
      factors = as.matrix(factors_em),
      h       = h,
      ly      = ly,
      lf      = lf
    )

    last       <- train_full[nrow(train_full), ]          # observação OOS
    train_full <- train_full[1:(nrow(train_full) - h), ]  # remove últimas h linhas

    # Remove lags iniciais (primeiras maxlag linhas têm NAs por construção)
    maxlag     <- max(lf, ly)
    train_full <- train_full[(maxlag + 1):nrow(train_full), ]
    train_full <- as.matrix(train_full[complete.cases(train_full), ])

    if (nrow(train_full) < 15) next

    subset <- 1:nrow(train_full)

    # =========================================================
    # MODELO m=1: Ridge Plano (benchmark)
    # =========================================================
    CV <- cv.glmnet(
      x      = train_full[subset, -1],
      y      = train_full[subset,  1],
      family = "gaussian",
      alpha  = 0
    )
    mdl <- glmnet(
      x      = train_full[subset, -1],
      y      = train_full[subset,  1],
      family = "gaussian",
      alpha  = 0,
      lambda = CV$lambda.min
    )

    pred_lin <- predict(mdl, newx = t(as.matrix(last[-1])))
    forecast[t, h, hi, 1]    <- as.numeric(pred_lin) - last[1]
    hp_track[t, h, hi, 1, 1] <- CV$lambda.min

    # =========================================================
    # MODELO m=2: 2SRR — Two-Step Ridge Regression (Coulombe)
    # =========================================================
    aa <- tryCatch(
      TVPRR_cosso(
        y         = train_full[subset, 1],
        X         = train_full[subset, -1],
        lambdavec = lambdavec,
        sweigths  = 1,
        type      = 2,
        alpha     = 0.01,
        silent    = silenceplz,
        kfold     = 5,
        lambda2   = CV$lambda.min,
        tol       = 1e-6,
        maxit     = 10,
        oosX      = last[-1]
      ),
      error = function(e) {
        cat(sprintf("  [WARN] TVPRR_cosso falhou t=%d h=%d: %s\n",
                    t, h, conditionMessage(e)))
        NULL
      }
    )

    if (!is.null(aa)) {
      pred_2srr <- OF(
        pred       = aa$fcast,
        y          = train_full[subset, 1],
        go.to.pred = pred_lin
      )
      forecast[t, h, hi, 2]    <- as.numeric(pred_2srr) - last[1]
      hp_track[t, h, hi, 2, 1] <- CV$lambda.min
      hp_track[t, h, hi, 2, 2] <- aa$grrats$lambdas[1]

      # Salva betas TVP para análise posterior
      betas_list[[hi]][[idx_oos]] <- list(
        t     = t,
        betas = aa$grrats$betas_grr,
        date  = if (exists("date")) date[t] else t
      )
    }

    if (idx_oos %% 12 == 0) {
      elapsed <- (proc.time() - t_start)["elapsed"]
      cat(sprintf("  h=%d | t=%d/%d (%.0f%%) | %.1f min\n",
                  h, t, bigt, 100 * idx_oos / n_oos, elapsed / 60))
    }
  }
}

cat(sprintf("\nPOOS concluído em %.1f min\n",
            (proc.time() - t_start)["elapsed"] / 60))

# ============================================================
# 9. SALVA OBJETOS PRINCIPAIS
# ============================================================
save(forecast,  file = "forecasts/coulombe_forecast_array.rda")
save(hp_track,  file = "forecasts/coulombe_hp_track.rda")
save(betas_list,file = "forecasts/coulombe_betas_list.rda")

# CSVs por horizonte (compatíveis com 07_analysis_coulombe.R)
for (hi in seq_along(hor)) {
  h   <- hor[hi]
  idx <- (tau + 1):bigt
  df  <- data.frame(
    t_idx   = idx,
    fc_ridge = forecast[idx, h, hi, 1],
    fc_2srr  = forecast[idx, h, hi, 2]
  )
  fname <- sprintf("forecasts/coulombe_fc_h%d.csv", h)
  write.csv(df, fname, row.names = FALSE)
  cat(sprintf("Salvo: %s\n", fname))
}

# ============================================================
# 10. TABELA DE SENSIBILIDADE: MSFE por nf
#     Roda um POOS rápido (apenas primeiros 30% do OOS) para
#     diferentes valores de nf, comparando MSFE do Ridge.
#     Serve como análise de robustez para o TCC.
# ============================================================
cat("\n=== SENSIBILIDADE: MSFE por nf ===\n")
nf_grid    <- c(2, 3, 4, 6, 8, 10)
tau_val    <- tau
end_val    <- tau + round(n_oos * 0.30)   # 30% do OOS como validação
hor_sens   <- c(1, 3, 6, 12)

msfe_mat <- matrix(NA, nrow = length(nf_grid), ncol = length(hor_sens),
                   dimnames = list(paste0("nf=", nf_grid),
                                   paste0("h=", hor_sens)))

for (nfi in seq_along(nf_grid)) {
  nf_try <- nf_grid[nfi]
  cat(sprintf("  nf=%d ...\n", nf_try))

  for (hi in seq_along(hor_sens)) {
    h    <- hor_sens[hi]
    errs <- c()

    for (t in tau_val:min(end_val, bigt - h)) {
      data_Xraw <- as.matrix(fred2)[1:(t - h), ]
      data_Xfac <- data_Xraw[, -cpi_col]

      factors_try <- tryCatch(
        EM_sw(data = as.data.frame(data_Xfac), n = nf_try, it_max = 500)$factors,
        error = function(e) NULL
      )
      if (is.null(factors_try)) next

      y_s <- as.matrix(forecast_vars[1:(t - h), hi])
      Y_s <- as.matrix(y_raw[1:(t - h)])

      start_s <- sum(is.na(y_s)) + 1
      end_s   <- length(1:(t - h))
      if ((end_s - start_s + 1) < 20) next

      y_s        <- y_s[start_s:end_s, 1, drop = FALSE]
      Y_s        <- Y_s[start_s:end_s, 1, drop = FALSE]
      factors_try <- factors_try[start_s:end_s, , drop = FALSE]

      tr <- tryCatch({
        tm <- make_reg_matrix(y=y_s, Y=Y_s,
                              factors=as.matrix(factors_try),
                              h=h, ly=ly, lf=lf)
        last_s <- tm[nrow(tm), ]
        tm <- tm[1:(nrow(tm)-h), ]
        tm <- tm[(max(lf,ly)+1):nrow(tm), ]
        tm <- as.matrix(tm[complete.cases(tm), ])
        list(tm=tm, last=last_s)
      }, error = function(e) NULL)
      if (is.null(tr) || nrow(tr$tm) < 10) next

      CV_s <- tryCatch(
        cv.glmnet(x=tr$tm[,-1], y=tr$tm[,1], family="gaussian", alpha=0),
        error = function(e) NULL
      )
      if (is.null(CV_s)) next

      mdl_s  <- glmnet(x=tr$tm[,-1], y=tr$tm[,1],
                        family="gaussian", alpha=0, lambda=CV_s$lambda.min)
      pred_s <- predict(mdl_s, newx=t(as.matrix(tr$last[-1]))) - tr$last[1]

      # Realizado = forecast_vars[t, hi]
      real_s <- forecast_vars[t, hi]
      if (!is.na(real_s) && !is.na(pred_s)) {
        errs <- c(errs, (pred_s - real_s)^2)
      }
    }
    msfe_mat[nfi, hi] <- if (length(errs) > 0) mean(errs) else NA
  }
}

# RMSFE relativo ao nf=8 (como tabela de sensibilidade)
ref_row       <- which(nf_grid == 8)
rmsfe_rel_mat <- sweep(sqrt(msfe_mat), 2, sqrt(msfe_mat[ref_row, ]), "/")

cat("\n--- RMSFE Relativo ao nf=8 (Ridge, validação 30% OOS inicial) ---\n")
print(round(rmsfe_rel_mat, 4))

df_sens <- as.data.frame(cbind(nf = nf_grid,
                                round(msfe_mat, 6),
                                round(rmsfe_rel_mat, 4)))
colnames(df_sens) <- c("nf",
                        paste0("MSFE_h", hor_sens),
                        paste0("RMSFE_rel_h", hor_sens))
write.csv(df_sens, "forecasts/coulombe_nf_sensitivity.csv", row.names = FALSE)
cat("Salvo: forecasts/coulombe_nf_sensitivity.csv\n")

cat("\n06_coulombe_2SRR_pipeline.R — COMPLETO\n")