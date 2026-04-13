# ============================================================
# 02_04_forecast.R
#
# Forecasting com MV-2SRR (Coulombe 2022)
# Adaptado para dados brasileiros 1996-2025
#
# Estrutura baseada em:
#   - Coulombe (Empirical/00_prog/forecasting_table16to17.R)
#   - Nathalia Oreda (Rcode/02_call_model.R)
#
# Output por combinação (M, V, H):
#   forecast[t, h, v, m]     — previsão POOS
#   betas_track[t, h, v, m,] — coeficientes TVP no tempo t
#   hp_track[t, h, v, m,]    — hiperparâmetros otimizados
# ============================================================

rm(list = ls())

# ============================================================
# PACKAGES
# ============================================================

myPKGs <- c("dplyr", "pracma", "doParallel", "foreach",
            "glmnet", "timeSeries", "matrixcalc", "zoo")

InstalledPKGs    <- names(installed.packages()[, "Package"])
InstallThesePKGs <- myPKGs[!myPKGs %in% InstalledPKGs]
if (length(InstallThesePKGs) > 0)
  install.packages(InstallThesePKGs, repos = "http://cran.us.r-project.org")

invisible(lapply(myPKGs, library, character.only = TRUE))

# ============================================================
# SETUP
# ============================================================

wd <- "C:/Users/felip/OneDrive/Documentos/tcc/Felipe_Dornelles_tcc/"
setwd(wd)

paths <- list(
  program   = "00_program",
  data      = "10_data",
  tools     = "20_tools",
  functions = "20_tools/functions",
  output    = "30_output",
  results   = "40_results"
)

# Encontra a pasta de dados mais recente
data_dirs  <- list.dirs(paths$data, full.names = TRUE, recursive = FALSE)
run_folder <- tail(sort(data_dirs), 1)
cat(sprintf("Loading data from: %s\n", run_folder))

# ============================================================
# 1. CARREGA DADOS
# ============================================================

load(file.path(run_folder, "df_model.rda"))
load(file.path(run_folder, "df_targets.rda"))
load(file.path(run_folder, "df_panel_pca.rda"))
load(file.path(run_folder, "targets_br.rda"))
load(file.path(run_folder, "all_options.rda"))

# Converte para matrizes numéricas (sem coluna date)
target_names <- unname(unlist(targets_br))   # c("PIB","IPCA","SELIC","CAMBIO","DESEMPREGO")
szv          <- length(target_names)          # 5

mat_targets <- df_targets |>
  select(all_of(target_names)) |>
  as.matrix()

mat_panel <- df_panel_pca |>
  select(-date) |>
  as.matrix()

dates_vec <- df_model$date

bigt <- nrow(mat_targets)
cat(sprintf("T = %d obs  (%s to %s)\n",
            bigt,
            format(min(dates_vec), "%b/%Y"),
            format(max(dates_vec), "%b/%Y")))
cat(sprintf("N panel variables: %d\n", ncol(mat_panel)))

# ============================================================
# 2. FONTES DO MV2SRR
# ============================================================

# Source das ferramentas do Coulombe (não modificar esses arquivos)
source(file.path(paths$tools,     "EM_sw.R"))
source(file.path(paths$tools,     "ICp2.R"))
source(file.path(paths$tools,     "Xgenerators_v190127.R"))
source(file.path(paths$functions, "dualGRRmdA_v190215.R"))
source(file.path(paths$functions, "CVGSBHK_v181127.R"))
source(file.path(paths$functions, "zfun_v190304.R"))
source(file.path(paths$functions, "factor.R"))
source(file.path(paths$functions, "TVPRRcosso_v181120.R"))
source(file.path(paths$functions, "TVPRRcossoF_v190125.R"))
source(file.path(paths$functions, "TVPRR_v181111.R"))
source(file.path(paths$functions, "fastZrot_v181125.R"))
source(file.path(paths$functions, "CVKFMV_v190214.R"))
source(file.path(paths$functions, "TVPRR_VARF_v190304.R"))
source(file.path(paths$functions, "MV2SRR_v221103.R"))  

# ============================================================
# 3. PARÂMETROS POOS
# ============================================================

# Ponto de início do POOS — 40% dos dados como burn-in mínimo
tau  <- round(bigt * 0.40)
cat(sprintf("POOS start: t = %d  (%s)\n", tau, format(dates_vec[tau], "%b/%Y")))

# Lags máximos
ly <- 2    # lags da variável dependente
lf <- 2    # lags dos fatores

# Fatores máximos para IM_sw
nf <- 3

# Re-otimizar hiperparâmetros a cada `jump` períodos (1 = toda iteração)
jump <- 12   # anual — igual ao Coulombe

# Grid de lambdas (igual ao Coulombe)
lambdavec <- exp(pracma::linspace(-2, 12, n = 15))

# Parâmetros do modelo 4 (Dense TVPs)
scree_ts  <- 0.05
maxf      <- 3
alpha_m4  <- 0.15
sv.param  <- 0

# Silencia outputs intermediários
silenceplz <- 1

# Sequências de otimização e "heranças"
optim_seq <- seq(tau, bigt, by = jump)
compl_seq <- setdiff(tau:bigt, optim_seq)

# ============================================================
# 4. FUNÇÃO DE OUTLIER (idêntica ao Coulombe)
# ============================================================

OF <- function(pred, y, tol = 2, go.to.pred) {
  newx      <- pred
  cond.max  <- (newx - mean(y)) >  tol * (max(y) - mean(y))
  cond.min  <- (newx - mean(y)) <  tol * (min(y) - mean(y))
  newx[cond.max] <- go.to.pred[cond.max]
  newx[cond.min] <- go.to.pred[cond.min]
  newx
}

# ============================================================
# 5. FUNÇÃO POOS — UMA COMBINAÇÃO (it_pos)
# ============================================================

POOS <- function(it_pos) {

  # Recarrega pacotes (necessário no worker paralelo)
  invisible(lapply(c("dplyr","pracma","glmnet","timeSeries",
                     "matrixcalc","zoo"), library,
                   character.only = TRUE))

  # Recarrega sources no worker
  source(file.path(paths$tools,     "EM_sw.R"))
  source(file.path(paths$tools,     "ICp2.R"))
  source(file.path(paths$tools,     "Xgenerators_v190127.R"))
  source(file.path(paths$functions, "dualGRRmdA_v190215.R"))
  source(file.path(paths$functions, "CVGSBHK_v181127.R"))
  source(file.path(paths$functions, "zfun_v190304.R"))
  source(file.path(paths$functions, "factor.R"))
  source(file.path(paths$functions, "TVPRRcosso_v181120.R"))
  source(file.path(paths$functions, "TVPRRcossoF_v190125.R"))
  source(file.path(paths$functions, "TVPRR_v181111.R"))
  source(file.path(paths$functions, "fastZrot_v181125.R"))
  source(file.path(paths$functions, "CVKFMV_v190214.R"))
  source(file.path(paths$functions, "TVPRR_VARF_v190304.R"))
  source(file.path(paths$functions, "MV2SRR_v221103.R"))

  h   <- all_options$H[it_pos]
  mod <- all_options$M[it_pos]
  v   <- all_options$V[it_pos]   # índice numérico 1..5

  cat(sprintf("[it=%d] M=%d  V=%d (%s)  H=%d\n",
              it_pos, mod, v, target_names[v], h))

  # Caminho de saída
  sname <- file.path(
    paths$output,
    sprintf("TVPfcst_V%d_H%d_M%d.RData", v, h, mod)
  )

  # Inicializa arrays de resultados
  # Dimensões: tempo × horizonte_max × vars × modelos × HPs/betas
  n_hp    <- 150
  n_betas <- ncol(mat_panel) * 2 + 10   # margem generosa

  forecast    <- array(NA_real_, dim = c(bigt, max(all_options$H), szv, 4))
  hp_track    <- array(NA_real_, dim = c(bigt, max(all_options$H), szv, 4, n_hp))
  betas_track <- array(NA_real_, dim = c(bigt, max(all_options$H), szv, 4, n_betas))

  # Imputa NaN do painel via EM antes do loop
  panel_full <- tryCatch(
    EM_sw(data = mat_panel, n = 8, it_max = 1000)$data,
    error = function(e) {
      warning(sprintf("EM_sw falhou: %s — usando mat_panel sem imputação.", e$message))
      mat_panel
    }
  )

  # hor_pos: mapeamento de H para posição no vetor de alvos do Coulombe
  hor_pos <- c(1, 2, 3, 3)   # H=1,2,4 → posição 1,2,3,3

  # -------------------------------------------------------
  # LOOP POOS
  # -------------------------------------------------------
  for (t in tau:bigt) {

    # Janela de informação disponível em t-h
    idx_end <- t - h
    if (idx_end < (ly + lf + 5)) next  # burn-in mínimo

    y_full <- mat_targets[, v]          # série completa do target v
    Y_full <- mat_targets[, v]          # regressor AR (mesma variável)
    factors_full <- panel_full

    # Corta até t-h
    y       <- y_full[1:idx_end]
    Y       <- Y_full[1:idx_end]
    factors <- factors_full[1:idx_end, , drop = FALSE]

    # Remove NAs no início da série dependente
    start_ok <- sum(is.na(y)) + 1
    end_ok   <- idx_end

    if (end_ok - start_ok < (ly + lf + 10)) next  # amostra muito curta

    y       <- y[start_ok:end_ok]
    Y       <- Y[start_ok:end_ok]
    factors <- factors[start_ok:end_ok, , drop = FALSE]

    # ---- Seleciona regressores por modelo (seguindo Coulombe) ----

    if (mod == 1) {
      # AR puro: só lags de y
      factors_use <- NULL
      train <- tryCatch(
        make_reg_matrix(y = y, Y = Y,
                        factors = matrix(0, length(y), 1),
                        h = h, ly = ly, lf = 0)[, 1:3],
        error = function(e) NULL
      )
    }
    if (mod == 2) {
      # 2 fatores PCA (excluindo a própria variável)
      fac2 <- tryCatch(
        EM_sw(data = factors[, -v, drop = FALSE], n = 2, it_max = 1000)$factors,
        error = function(e) factors[, 1:min(2, ncol(factors)), drop = FALSE]
      )
      train <- tryCatch(
        make_reg_matrix(y = y, Y = Y, factors = fac2, h = h, ly = ly, lf = lf),
        error = function(e) NULL
      )
    }
    if (mod == 3) {
      # Outros 4 targets como fatores
      other_targets <- mat_targets[1:idx_end, -v, drop = FALSE]
      other_clean   <- other_targets[start_ok:end_ok, , drop = FALSE]
      train <- tryCatch(
        make_reg_matrix(y = y, Y = Y, factors = other_clean,
                        h = h, ly = ly, lf = lf),
        error = function(e) NULL
      )
    }
    if (mod == 4) {
      # Painel completo (modelo "medium" do Coulombe)
      train <- tryCatch(
        make_reg_matrix(y = y, Y = Y, factors = factors,
                        h = h, ly = ly, lf = lf),
        error = function(e) NULL
      )
    }

    if (is.null(train)) next

    # Última linha = vetor de previsão (last)
    last  <- train[nrow(train), , drop = FALSE]
    train <- train[1:(nrow(train) - h), , drop = FALSE]

    # Remove lags iniciais e linhas incompletas
    maxlag <- max(lf, ly)
    if (nrow(train) <= maxlag + h + 2) next
    train  <- train[(maxlag + 1):nrow(train), , drop = FALSE]
    train  <- as.matrix(train[complete.cases(train), , drop = FALSE])

    if (nrow(train) < 20) next   # mínimo de observações

    subset_idx <- 1:nrow(train)

    # -------------------------------------------------------
    # MODELO 1 (M=1) — Ridge plain (cv.glmnet, alpha=0)
    # -------------------------------------------------------
    tryCatch({
      CV  <- cv.glmnet(x = train[subset_idx, -1, drop = FALSE],
                       y = train[subset_idx,  1],
                       family = "gaussian", alpha = 0)
      mdl <- glmnet(x = train[subset_idx, -1, drop = FALSE],
                    y = train[subset_idx,  1],
                    family = "gaussian", alpha = 0,
                    lambda = CV$lambda.min)

      pred_ridge <- as.numeric(
        predict(mdl, newx = as.matrix(last[, -1, drop = FALSE]))
      )
      forecast[t, h, v, 1]    <- pred_ridge - last[1, 1]
      hp_track[t, h, v, 1, 1] <- CV$lambda.min

      # Betas do ridge (coeficientes estáticos)
      beta_ridge <- as.numeric(coef(mdl))
      n_b <- min(length(beta_ridge), n_betas)
      betas_track[t, h, v, 1, 1:n_b] <- beta_ridge[1:n_b]

      pred_lin <- pred_ridge   # referência para OF
    }, error = function(e)
      cat(sprintf("  M1 t=%d v=%d h=%d: %s\n", t, v, h, e$message))
    )

    # -------------------------------------------------------
    # MODELO 2 (M=2) — 2SRR (TVP esparso)
    # -------------------------------------------------------
    lv2 <- if (exists("CV")) CV$lambda.min else 1.0

    tryCatch({
      lv2_use <- if (exists("CV") && !is.null(CV$lambda.min))
                   CV$lambda.min else 1.0

      aa <- TVPRR_cosso(
        y        = train[subset_idx, 1],
        X        = train[subset_idx, -1, drop = FALSE],
        lambdavec = lambdavec,
        sweigths  = 1,
        type      = 2,
        alpha     = 0.01,
        silent    = silenceplz,
        kfold     = 5,
        lambda2   = lv2_use,
        tol       = 1e-6,
        maxit     = 10,
        oosX      = as.numeric(last[, -1])
      )

      pred_2srr <- OF(
        pred       = aa$fcast,
        y          = train[subset_idx, 1],
        go.to.pred = if (exists("pred_lin")) pred_lin else aa$fcast
      )

      forecast[t, h, v, 2]    <- pred_2srr - last[1, 1]
      hp_track[t, h, v, 2, 1] <- lv2_use
      hp_track[t, h, v, 2, 2] <- aa$grrats$lambdas[1]

      # --- EXTRAÇÃO DOS BETAS TVP (Coulombe: aa$grrats$betas_grr) ---
      # betas_grr: matriz (T_train x K) — coeficientes que variam no tempo
      # Pegamos a última linha = coeficiente no instante t
      if (!is.null(aa$grrats) && !is.null(aa$grrats$betas_grr)) {
        betas_t <- as.numeric(aa$grrats$betas_grr[nrow(aa$grrats$betas_grr), ])
        n_b <- min(length(betas_t), n_betas)
        betas_track[t, h, v, 2, 1:n_b] <- betas_t[1:n_b]
      }
      # Betas do passo 2 (menos suavizados)
      # aa$grr$betas_grr: também disponível se quiser comparar
    }, error = function(e)
      cat(sprintf("  M2 t=%d v=%d h=%d: %s\n", t, v, h, e$message))
    )

    # -------------------------------------------------------
    # MODELO 3 (M=3) — MSRRs (TVP esparso multi-step)
    # -------------------------------------------------------
    tryCatch({
      lv3_use <- if (exists("CV") && !is.null(CV$lambda.min))
                   CV$lambda.min else 1.0
      lambdavec_m3 <- if (mod == 4) lambdavec[4:15] else lambdavec

      aa3 <- TVPRR(
        y        = train[subset_idx, 1],
        X        = train[subset_idx, -1, drop = FALSE],
        lambdavec = lambdavec_m3,
        sweigths  = 1,
        type      = 3,
        alpha     = 0.001,
        silent    = silenceplz,
        kfold     = 5,
        lambda2   = lv3_use,
        tol       = 1e-5,
        maxit     = 15,
        oosX      = as.numeric(last[, -1])
      )

      pred_m3 <- OF(
        pred       = aa3$fcast,
        y          = train[subset_idx, 1],
        go.to.pred = if (exists("pred_lin")) pred_lin else aa3$fcast
      )

      forecast[t, h, v, 3]    <- pred_m3 - last[1, 1]
      hp_track[t, h, v, 3, 1] <- lv3_use
      hp_track[t, h, v, 3, 2] <- aa3$lambda1

      if (!is.null(aa3$grra) && !is.null(aa3$grra$sigmasq)) {
        n_s <- min(length(aa3$grra$sigmasq), n_hp - 2)
        hp_track[t, h, v, 3, 3:(2 + n_s)] <- aa3$grra$sigmasq[1:n_s]
      }

      # Betas TVP do MSRRs
      if (!is.null(aa3$grra) && !is.null(aa3$grra$betas_grr)) {
        betas_t3 <- as.numeric(aa3$grra$betas_grr[nrow(aa3$grra$betas_grr), ])
        n_b <- min(length(betas_t3), n_betas)
        betas_track[t, h, v, 3, 1:n_b] <- betas_t3[1:n_b]
      }
    }, error = function(e)
      cat(sprintf("  M3 t=%d v=%d h=%d: %s\n", t, v, h, e$message))
    )

    # -------------------------------------------------------
    # MODELO 4 (M=4) — MSRRd (Dense TVPs via VAR-F)
    # -------------------------------------------------------
    tryCatch({
      lv4_use <- if (exists("CV") && !is.null(CV$lambda.min))
                   CV$lambda.min else 1.0
      lambdavec_m4 <- lambdavec[4:15]

      aa4 <- TVPRR_VARF(
        Y          = train[subset_idx, 1],
        X          = train[subset_idx, -1, drop = FALSE],
        orthoFac   = TRUE,
        lambdavec  = lambdavec_m4,
        sweigths   = 1,
        type       = 2,
        fp.model   = 1,
        sv.param   = sv.param,
        alpha      = alpha_m4,
        silent     = silenceplz,
        kfold      = 5,
        lambda2    = lv4_use,
        max.step.cv = 8,
        adaptive   = 1,
        aparam     = -0.5,
        tol        = 1e-10,
        maxit      = 20,
        lambdabooster = 1,
        var.share  = scree_ts,
        override   = maxf,
        id         = 1,
        oosX       = as.numeric(last[, -1])
      )

      pred_m4 <- OF(
        pred       = aa4$fcast,
        y          = train[subset_idx, 1],
        go.to.pred = if (exists("pred_lin")) pred_lin else aa4$fcast
      )

      forecast[t, h, v, 4] <- pred_m4 - last[1, 1]

      if (!is.null(aa4$HPs)) {
        n_hp4 <- min(length(aa4$HPs), n_hp)
        hp_track[t, h, v, 4, 1:n_hp4] <- aa4$HPs[1:n_hp4]
      }

      # Betas TVP do Dense — aa4$betas ou aa4$starter_pack$betas
      betas_src <- NULL
      if (!is.null(aa4$betas)) {
        betas_src <- aa4$betas
      } else if (!is.null(aa4$starter_pack) && !is.null(aa4$starter_pack$betas)) {
        betas_src <- aa4$starter_pack$betas
      }
      if (!is.null(betas_src)) {
        betas_t4 <- if (is.matrix(betas_src)) {
          as.numeric(betas_src[nrow(betas_src), ])
        } else {
          as.numeric(betas_src)
        }
        n_b <- min(length(betas_t4), n_betas)
        betas_track[t, h, v, 4, 1:n_b] <- betas_t4[1:n_b]
      }
    }, error = function(e)
      cat(sprintf("  M4 t=%d v=%d h=%d: %s\n", t, v, h, e$message))
    )

    if (t %% 12 == 0)
      cat(sprintf("  ... t=%d (%s) OK\n", t, format(dates_vec[t], "%b/%Y")))

  } # fim loop t

  # -------------------------------------------------------
  # SALVA RESULTADOS
  # -------------------------------------------------------
  save(forecast, hp_track, betas_track,
       file = sname)
  cat(sprintf("[it=%d] Salvo em: %s\n", it_pos, sname))

  return(invisible(NULL))
}

# ============================================================
# 6. EXECUTA EM PARALELO
# ============================================================

if (!dir.exists(paths$output)) dir.create(paths$output, recursive = TRUE)

ncores <- max(1L, parallel::detectCores() - 1L)
cat(sprintf("Usando %d cores\n", ncores))

if (ncores > 1) {
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  foreach(it_pos = 1:nrow(all_options),
          .packages = c("dplyr","pracma","glmnet","timeSeries","matrixcalc","zoo"),
          .export   = c("paths","mat_targets","mat_panel","panel_full",
                        "target_names","szv","dates_vec","bigt",
                        "all_options","tau","ly","lf","nf","jump",
                        "lambdavec","scree_ts","maxf","alpha_m4",
                        "sv.param","silenceplz","n_hp","n_betas",
                        "OF","POOS")) %dopar% POOS(it_pos)
  stopCluster(cl)
} else {
  for (it_pos in 1:nrow(all_options)) POOS(it_pos)
}

cat("\n=== Forecast concluido ===\n")
cat(sprintf("Resultados em: %s/\n", paths$output))

# ============================================================
# 7. CONSOLIDA — JUNTA TODOS OS .RData EM UM ÚNICO OBJETO
# ============================================================

cat("\n=== Consolidando resultados ===\n")

forecast_all    <- array(NA_real_, dim = c(bigt, max(all_options$H), szv, 4))
betas_all       <- array(NA_real_, dim = c(bigt, max(all_options$H), szv, 4, 200))

for (it_pos in 1:nrow(all_options)) {
  h   <- all_options$H[it_pos]
  mod <- all_options$M[it_pos]
  v   <- all_options$V[it_pos]

  sname <- file.path(
    paths$output,
    sprintf("TVPfcst_V%d_H%d_M%d.RData", v, h, mod)
  )
  if (!file.exists(sname)) next

  local({
    e <- new.env()
    load(sname, envir = e)
    forecast_all[, h, v, mod]       <<- e$forecast[, h, v, mod]
    n_b <- min(dim(e$betas_track)[5], 200)
    betas_all[, h, v, mod, 1:n_b]  <<- e$betas_track[, h, v, mod, 1:n_b]
  })
}

# ============================================================
# 8. MONTA DATA.FRAME DE BETAS TVP (para gráficos do TCC)
# ============================================================
# Estrutura: uma linha por (t, h, v, m, k)
# t = data, h = horizonte, v = variável, m = modelo, k = coeficiente

cat("  Construindo df_betas_long...\n")

beta_rows <- list()
k_max     <- 10   # captura os 10 primeiros coeficientes (intercepto + 9 regressores)

for (it_pos in 1:nrow(all_options)) {
  h   <- all_options$H[it_pos]
  mod <- all_options$M[it_pos]
  v   <- all_options$V[it_pos]

  for (t in tau:bigt) {
    betas_t <- betas_all[t, h, v, mod, 1:k_max]
    if (all(is.na(betas_t))) next

    beta_rows[[length(beta_rows) + 1]] <- data.frame(
      date    = dates_vec[t],
      h       = h,
      var     = target_names[v],
      model   = mod,
      coef_id = 0:(k_max - 1),
      beta    = betas_t
    )
  }
}

df_betas_long <- do.call(rbind, beta_rows)
rownames(df_betas_long) <- NULL

cat(sprintf("  df_betas_long: %d linhas\n", nrow(df_betas_long)))

# ============================================================
# 9. MONTA DATA.FRAME DE FORECASTS
# ============================================================

cat("  Construindo df_forecasts...\n")

fc_rows <- list()

for (it_pos in 1:nrow(all_options)) {
  h   <- all_options$H[it_pos]
  mod <- all_options$M[it_pos]
  v   <- all_options$V[it_pos]

  for (t in tau:bigt) {
    fc_val <- forecast_all[t, h, v, mod]
    if (is.na(fc_val)) next

    # Valor realizado (t+h, se disponível)
    t_real <- t + h
    actual_val <- if (t_real <= bigt) mat_targets[t_real, v] else NA_real_

    fc_rows[[length(fc_rows) + 1]] <- data.frame(
      date_forecast = dates_vec[t],
      date_target   = if (t_real <= bigt) dates_vec[t_real] else NA,
      h             = h,
      var           = target_names[v],
      model         = mod,
      forecast      = fc_val,
      actual        = actual_val,
      error         = actual_val - fc_val
    )
  }
}

df_forecasts <- do.call(rbind, fc_rows)
rownames(df_forecasts) <- NULL

cat(sprintf("  df_forecasts: %d linhas\n", nrow(df_forecasts)))

# ============================================================
# 10. SALVA TUDO
# ============================================================

results_folder <- paths$results
if (!dir.exists(results_folder)) dir.create(results_folder, recursive = TRUE)

save(forecast_all,    file = file.path(results_folder, "forecast_all.rda"))
save(betas_all,       file = file.path(results_folder, "betas_all.rda"))
save(df_forecasts,    file = file.path(results_folder, "df_forecasts.rda"))
save(df_betas_long,   file = file.path(results_folder, "df_betas_long.rda"))

write.csv(df_forecasts,  file.path(results_folder, "df_forecasts.csv"),  row.names = FALSE)
write.csv(df_betas_long, file.path(results_folder, "df_betas_long.csv"), row.names = FALSE)

cat("\n=== Objetos salvos em 40_results/ ===\n")
cat("  forecast_all.rda   — array [T x H x V x M] de previsoes\n")
cat("  betas_all.rda      — array [T x H x V x M x K] de betas TVP\n")
cat("  df_forecasts.csv   — long format: data, h, var, model, forecast, actual, error\n")
cat("  df_betas_long.csv  — long format: data, h, var, model, coef_id, beta\n")