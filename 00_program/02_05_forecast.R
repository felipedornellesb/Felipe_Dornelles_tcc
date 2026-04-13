# ============================================================
# 02_05_forecast.R
#
# Script de previsão POOS com modelos TVP (Coulombe 2SRR)
# aplicado a séries macroeconômicas brasileiras
#
# Correções aplicadas:
#   1) panel_full, n_hp e n_betas criados no global ANTES do foreach
#   2) bloco EM_sw removido de dentro do POOS()
#   3) make_reg_matrix substituído pela lógica correta de
#      montagem de regressores (lag manual)
#   4) todos os objetos em .export existem no global
#   5) clusterExport explícito como segurança extra no Windows
# ============================================================

rm(list = ls())

# ============================================================
# 1. PACOTES
# ============================================================

myPKGs <- c("dplyr", "pracma", "doParallel", "foreach",
            "glmnet", "timeSeries", "matrixcalc", "zoo")

InstalledPKGs    <- names(installed.packages()[, "Package"])
InstallThesePKGs <- myPKGs[!myPKGs %in% InstalledPKGs]
if (length(InstallThesePKGs) > 0)
  install.packages(InstallThesePKGs, repos = "http://cran.us.r-project.org")

invisible(lapply(myPKGs, library, character.only = TRUE))

# ============================================================
# 2. CAMINHOS
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

# Pasta de dados mais recente
data_dirs  <- list.dirs(paths$data, full.names = TRUE, recursive = FALSE)
run_folder <- tail(sort(data_dirs), 1)
cat(sprintf("Loading data from: %s\n", run_folder))

# ============================================================
# 3. CARREGA DADOS
# ============================================================

load(file.path(run_folder, "df_model.rda"))
load(file.path(run_folder, "df_targets.rda"))
load(file.path(run_folder, "df_panel_pca.rda"))
load(file.path(run_folder, "targets_br.rda"))
load(file.path(run_folder, "all_options.rda"))

target_names <- unname(unlist(targets_br))   # c("PIB","IPCA","SELIC","CAMBIO","DESEMPREGO")
szv          <- length(target_names)          # 5

mat_targets <- df_targets |>
  select(all_of(target_names)) |>
  as.matrix()

mat_panel <- df_panel_pca |>
  select(-date) |>
  as.matrix()

dates_vec <- df_model$date
bigt      <- nrow(mat_targets)

cat(sprintf("T = %d obs  (%s to %s)\n",
            bigt,
            format(min(dates_vec), "%b/%Y"),
            format(max(dates_vec), "%b/%Y")))
cat(sprintf("N panel variables: %d\n", ncol(mat_panel)))

# ============================================================
# 4. SOURCE DAS FERRAMENTAS DO COULOMBE
# ============================================================

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
# 5. PARÂMETROS GLOBAIS
# ============================================================

tau        <- round(bigt * 0.40)
cat(sprintf("POOS start: t = %d  (%s)\n", tau, format(dates_vec[tau], "%b/%Y")))

ly         <- 2      # lags da variável dependente
lf         <- 2      # lags dos fatores
nf         <- 3      # fatores máximos para EM_sw
jump       <- 12     # re-otimiza HPs a cada 12 meses

lambdavec  <- exp(pracma::linspace(-2, 12, n = 15))

scree_ts   <- 0.05
maxf       <- 3
alpha_m4   <- 0.15
sv.param   <- 0
silenceplz <- 1

# ---- n_hp e n_betas criados AQUI, no global ---------------
n_hp    <- 150
n_betas <- ncol(mat_panel) * 2 + 10

# ============================================================
# 6. PRÉ-COMPUTA EM_sw NO GLOBAL (antes do paralelo)
# ============================================================
# CORREÇÃO PRINCIPAL: panel_full deve existir no ambiente global
# antes de qualquer chamada ao foreach/.export

cat("Imputando painel via EM_sw (pré-cômputo global)...\n")

panel_full <- tryCatch(
  EM_sw(data = mat_panel, n = nf, it_max = 1000)$data,
  error = function(e) {
    warning(sprintf("EM_sw global falhou: %s — usando mat_panel bruto.", e$message))
    mat_panel
  }
)

cat(sprintf("panel_full: %d x %d\n", nrow(panel_full), ncol(panel_full)))

# ============================================================
# 7. FUNÇÕES AUXILIARES
# ============================================================

# Monta matriz de regressores [y_h | lags_y | lags_factors]
# y     : vetor dependente (comprimento T)
# Y     : vetor usado para lags AR (geralmente igual a y)
# F     : matriz de fatores (T x K); pode ser NULL
# h     : horizonte
# ly    : nº de lags de y
# lf    : nº de lags dos fatores
# Retorna matriz (T x (1 + ly + K*lf)) — última linha é o vetor de previsão
build_reg <- function(y, Y, F = NULL, h, ly, lf) {

  TT <- length(y)

  # Target h-step-ahead: y deslocado h períodos à frente
  y_h <- c(y[(h + 1):TT], rep(NA_real_, h))

  # Lags de Y
  lag_mat_y <- matrix(NA_real_, TT, ly)
  for (j in 1:ly) {
    lag_mat_y[(j + 1):TT, j] <- Y[1:(TT - j)]
  }

  if (!is.null(F) && ncol(F) > 0 && lf > 0) {
    # Lags dos fatores
    lag_mat_f <- matrix(NA_real_, TT, ncol(F) * lf)
    col_idx <- 1
    for (k in 1:ncol(F)) {
      for (j in 1:lf) {
        lag_mat_f[(j + 1):TT, col_idx] <- F[1:(TT - j), k]
        col_idx <- col_idx + 1
      }
    }
    reg <- cbind(y_h, lag_mat_y, lag_mat_f)
  } else {
    reg <- cbind(y_h, lag_mat_y)
  }

  colnames(reg)[1] <- "y_h"
  reg
}

# Função de outlier-filter: substitui previsões absurdas pela referência linear
OF <- function(pred, y, tol = 2, go.to.pred) {
  newx     <- pred
  cond.max <- (newx - mean(y)) >  tol * (max(y) - mean(y))
  cond.min <- (newx - mean(y)) <  tol * (min(y) - mean(y))
  newx[cond.max] <- go.to.pred[cond.max]
  newx[cond.min] <- go.to.pred[cond.min]
  newx
}

# ============================================================
# 8. FUNÇÃO POOS
#    — panel_full, n_hp, n_betas chegam via .export do global
#    — NÃO recria EM_sw aqui dentro
# ============================================================

POOS <- function(it_pos) {

  # Recarrega pacotes no worker
  invisible(lapply(c("dplyr", "pracma", "glmnet", "timeSeries",
                     "matrixcalc", "zoo"),
                   library, character.only = TRUE))

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
  v   <- all_options$V[it_pos]

  cat(sprintf("[it=%d] M=%d  V=%d (%s)  H=%d\n",
              it_pos, mod, v, target_names[v], h))

  # Pasta de output por data
  out_dir <- file.path(
    paths$output,
    paste0("outputs_", format(Sys.Date(), "%m_%d_%Y"))
  )
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  sname <- file.path(
    out_dir,
    sprintf("TVPfcst_V%d_H%d_M%d.RData", v, h, mod)
  )

  # Arrays de resultados
  forecast    <- array(NA_real_, dim = c(bigt, max(all_options$H), szv, 4))
  hp_track    <- array(NA_real_, dim = c(bigt, max(all_options$H), szv, 4, n_hp))
  betas_track <- array(NA_real_, dim = c(bigt, max(all_options$H), szv, 4, n_betas))

  # panel_full já vem do .export — não recria aqui
  factors_full <- panel_full   # (bigt x K), imputado pelo EM_sw global

  # -------------------------------------------------------
  # LOOP POOS
  # -------------------------------------------------------
  for (t in tau:bigt) {

    idx_end <- t - h
    if (idx_end < (max(ly, lf) + h + 5)) next

    y_full <- mat_targets[, v]

    # Corta até idx_end
    y_win       <- y_full[1:idx_end]
    factors_win <- factors_full[1:idx_end, , drop = FALSE]

    # Remove NAs iniciais
    start_ok <- which(!is.na(y_win))[1]
    if (is.na(start_ok)) next
    end_ok <- idx_end

    if ((end_ok - start_ok + 1) < (max(ly, lf) + h + 10)) next

    y_win       <- y_win[start_ok:end_ok]
    factors_win <- factors_win[start_ok:end_ok, , drop = FALSE]

    # ---- Monta regressores por modelo ----

    reg <- NULL

    if (mod == 1) {
      # AR puro — sem fatores externos
      reg <- tryCatch(
        build_reg(y = y_win, Y = y_win,
                  F = NULL, h = h, ly = ly, lf = 0),
        error = function(e) NULL
      )
    }

    if (mod == 2) {
      # 2 fatores PCA do painel (excluindo a própria variável se possível)
      col_idx_panel <- seq_len(ncol(factors_win))
      fac2 <- tryCatch({
        EM_sw(data = factors_win[, col_idx_panel, drop = FALSE],
              n = 2, it_max = 1000)$factors
      }, error = function(e) {
        factors_win[, 1:min(2, ncol(factors_win)), drop = FALSE]
      })
      reg <- tryCatch(
        build_reg(y = y_win, Y = y_win,
                  F = fac2, h = h, ly = ly, lf = lf),
        error = function(e) NULL
      )
    }

    if (mod == 3) {
      # Outros targets como fatores
      other_idx    <- setdiff(seq_len(szv), v)
      other_full   <- mat_targets[1:idx_end, other_idx, drop = FALSE]
      other_win    <- other_full[start_ok:end_ok, , drop = FALSE]
      reg <- tryCatch(
        build_reg(y = y_win, Y = y_win,
                  F = other_win, h = h, ly = ly, lf = lf),
        error = function(e) NULL
      )
    }

    if (mod == 4) {
      # Painel completo
      reg <- tryCatch(
        build_reg(y = y_win, Y = y_win,
                  F = factors_win, h = h, ly = ly, lf = lf),
        error = function(e) NULL
      )
    }

    if (is.null(reg)) next

    # Última linha = vetor de previsão (sem y_h disponível)
    last  <- reg[nrow(reg), , drop = FALSE]
    train <- reg[1:(nrow(reg) - h), , drop = FALSE]

    # Remove linhas com NA
    train <- as.matrix(train[complete.cases(train), , drop = FALSE])
    if (nrow(train) < 20) next

    subset_idx <- seq_len(nrow(train))

    X_train <- train[subset_idx, -1, drop = FALSE]
    y_train <- train[subset_idx,  1]
    x_last  <- as.numeric(last[, -1])

    # -------------------------------------------------------
    # M1 — Ridge (cv.glmnet, alpha = 0)
    # -------------------------------------------------------
    pred_lin <- NA_real_

    tryCatch({
      CV  <- cv.glmnet(x = X_train, y = y_train,
                       family = "gaussian", alpha = 0)
      mdl <- glmnet(x = X_train, y = y_train,
                    family = "gaussian", alpha = 0,
                    lambda = CV$lambda.min)

      pred_ridge <- as.numeric(
        predict(mdl, newx = matrix(x_last, nrow = 1))
      )
      pred_lin <- pred_ridge   # referência para OF dos modelos TVP

      forecast[t, h, v, 1]    <- pred_ridge - last[1, 1]
      hp_track[t, h, v, 1, 1] <- CV$lambda.min

      beta_ridge <- as.numeric(coef(mdl))
      n_b <- min(length(beta_ridge), n_betas)
      betas_track[t, h, v, 1, 1:n_b] <- beta_ridge[1:n_b]

    }, error = function(e)
      cat(sprintf("  M1 t=%d v=%d h=%d: %s\n", t, v, h, e$message))
    )

    lv_ref <- if (!is.na(hp_track[t, h, v, 1, 1])) hp_track[t, h, v, 1, 1] else 1.0

    # -------------------------------------------------------
    # M2 — 2SRR (TVP esparso, TVPRR_cosso)
    # -------------------------------------------------------
    tryCatch({
      aa2 <- TVPRR_cosso(
        y         = y_train,
        X         = X_train,
        lambdavec = lambdavec,
        sweigths  = 1,
        type      = 2,
        alpha     = 0.01,
        silent    = silenceplz,
        kfold     = 5,
        lambda2   = lv_ref,
        tol       = 1e-6,
        maxit     = 10,
        oosX      = x_last
      )

      pred_2srr <- OF(
        pred       = aa2$fcast,
        y          = y_train,
        go.to.pred = if (!is.na(pred_lin)) pred_lin else aa2$fcast
      )

      forecast[t, h, v, 2]    <- pred_2srr - last[1, 1]
      hp_track[t, h, v, 2, 1] <- lv_ref
      if (!is.null(aa2$grrats$lambdas))
        hp_track[t, h, v, 2, 2] <- aa2$grrats$lambdas[1]

      if (!is.null(aa2$grrats$betas_grr)) {
        betas_t2 <- as.numeric(
          aa2$grrats$betas_grr[nrow(aa2$grrats$betas_grr), ]
        )
        n_b <- min(length(betas_t2), n_betas)
        betas_track[t, h, v, 2, 1:n_b] <- betas_t2[1:n_b]
      }

    }, error = function(e)
      cat(sprintf("  M2 t=%d v=%d h=%d: %s\n", t, v, h, e$message))
    )

    # -------------------------------------------------------
    # M3 — MSRRs (TVP esparso multi-step, TVPRR)
    # -------------------------------------------------------
    tryCatch({
      lambdavec_m3 <- if (mod == 4) lambdavec[4:15] else lambdavec

      aa3 <- TVPRR(
        y         = y_train,
        X         = X_train,
        lambdavec = lambdavec_m3,
        sweigths  = 1,
        type      = 3,
        alpha     = 0.001,
        silent    = silenceplz,
        kfold     = 5,
        lambda2   = lv_ref,
        tol       = 1e-5,
        maxit     = 15,
        oosX      = x_last
      )

      pred_m3 <- OF(
        pred       = aa3$fcast,
        y          = y_train,
        go.to.pred = if (!is.na(pred_lin)) pred_lin else aa3$fcast
      )

      forecast[t, h, v, 3]    <- pred_m3 - last[1, 1]
      hp_track[t, h, v, 3, 1] <- lv_ref
      if (!is.null(aa3$lambda1))
        hp_track[t, h, v, 3, 2] <- aa3$lambda1

      if (!is.null(aa3$grra$sigmasq)) {
        n_s <- min(length(aa3$grra$sigmasq), n_hp - 2)
        hp_track[t, h, v, 3, 3:(2 + n_s)] <- aa3$grra$sigmasq[1:n_s]
      }

      if (!is.null(aa3$grra$betas_grr)) {
        betas_t3 <- as.numeric(
          aa3$grra$betas_grr[nrow(aa3$grra$betas_grr), ]
        )
        n_b <- min(length(betas_t3), n_betas)
        betas_track[t, h, v, 3, 1:n_b] <- betas_t3[1:n_b]
      }

    }, error = function(e)
      cat(sprintf("  M3 t=%d v=%d h=%d: %s\n", t, v, h, e$message))
    )

    # -------------------------------------------------------
    # M4 — MSRRd (Dense TVPs via VAR-F, TVPRR_VARF)
    # -------------------------------------------------------
    tryCatch({
      aa4 <- TVPRR_VARF(
        Y             = y_train,
        X             = X_train,
        orthoFac      = TRUE,
        lambdavec     = lambdavec[4:15],
        sweigths      = 1,
        type          = 2,
        fp.model      = 1,
        sv.param      = sv.param,
        alpha         = alpha_m4,
        silent        = silenceplz,
        kfold         = 5,
        lambda2       = lv_ref,
        max.step.cv   = 8,
        adaptive      = 1,
        aparam        = -0.5,
        tol           = 1e-10,
        maxit         = 20,
        lambdabooster = 1,
        var.share     = scree_ts,
        override      = maxf,
        id            = 1,
        oosX          = x_last
      )

      pred_m4 <- OF(
        pred       = aa4$fcast,
        y          = y_train,
        go.to.pred = if (!is.na(pred_lin)) pred_lin else aa4$fcast
      )

      forecast[t, h, v, 4] <- pred_m4 - last[1, 1]

      if (!is.null(aa4$HPs)) {
        n_hp4 <- min(length(aa4$HPs), n_hp)
        hp_track[t, h, v, 4, 1:n_hp4] <- aa4$HPs[1:n_hp4]
      }

      betas_src <- NULL
      if (!is.null(aa4$betas)) {
        betas_src <- aa4$betas
      } else if (!is.null(aa4$starter_pack$betas)) {
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

  }  # fim loop t

  # Salva resultados
  save(forecast, hp_track, betas_track, file = sname)
  cat(sprintf("[it=%d] Salvo em: %s\n", it_pos, sname))

  return(invisible(NULL))
}

# ============================================================
# 9. EXECUTA EM PARALELO
# ============================================================

out_dir_main <- file.path(
  paths$output,
  paste0("outputs_", format(Sys.Date(), "%m_%d_%Y"))
)
if (!dir.exists(out_dir_main)) dir.create(out_dir_main, recursive = TRUE)

ncores <- max(1L, parallel::detectCores() - 1L)
cat(sprintf("Usando %d cores\n", ncores))

# Lista de todos os objetos globais a exportar
global_exports <- c(
  # dados
  "paths", "mat_targets", "mat_panel", "panel_full",
  "target_names", "szv", "dates_vec", "bigt",
  "all_options",
  # parâmetros
  "tau", "ly", "lf", "nf", "jump",
  "lambdavec", "scree_ts", "maxf", "alpha_m4",
  "sv.param", "silenceplz",
  "n_hp", "n_betas",
  # funções locais
  "OF", "build_reg", "POOS"
)

if (ncores > 1) {

  cl <- makeCluster(ncores)
  registerDoParallel(cl)

  # clusterExport garante que os objetos cheguem aos workers no Windows
  clusterExport(cl, varlist = global_exports, envir = environment())

  foreach(
    it_pos    = 1:nrow(all_options),
    .packages = c("dplyr", "pracma", "glmnet",
                  "timeSeries", "matrixcalc", "zoo"),
    .export   = global_exports
  ) %dopar% POOS(it_pos)

  stopCluster(cl)

} else {

  for (it_pos in 1:nrow(all_options)) POOS(it_pos)

}

cat("\n=== POOS concluído ===\n")
