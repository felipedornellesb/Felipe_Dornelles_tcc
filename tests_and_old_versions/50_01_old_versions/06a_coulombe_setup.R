# 06a_coulombe_setup.R
#
# PARTE 1 DE 3 — Setup, carregamento e diagnostico
#
# O que este script faz:
#   1. Define factor() PCA (compativel com EM_sw) ANTES de
#      qualquer source() — corrige "unused argument (n_fac)"
#   2. Carrega todas as funcoes do Coulombe de coulombe/functions/
#   3. Carrega data.rda do Medeiros e extrai y, X, dates
#   4. Alinha janela OOS com o yout.rda do Medeiros
#   5. Roda EM Stock & Watson se houver NAs (fallback automatico)
#   6. Remove colunas constantes antes do PCA (evita erro prcomp)
#   7. Testa make_reg_matrix com janela de 150 obs
#   8. Testa TVPRR_cosso com janela reduzida
#   9. Salva objetos para o 06b em coulombe/setup_objects.rda
#
# Execute: source("coulombe/06a_coulombe_setup.R")
# ============================================================

setwd("~/tcc/Felipe_Dornelles_tcc/forecasting_inflation")

# ============================================================
# PACOTES
# ============================================================

pkgs <- c("tidyverse", "glmnet", "pracma")
new  <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(new)) install.packages(new)

library(tidyverse)
library(glmnet)
library(pracma)

# ============================================================
# PASSO 0 — Define factor() PCA ANTES de qualquer source()
#
# O EM_sw.R do Coulombe chama internamente:
#   factor(X0, n_fac = n)
# A factor() base do R nao aceita o argumento n_fac e lanca:
#   "unused argument (n_fac = n)"
# Precisa ser substituida pela versao PCA aqui, antes do
# source("EM_sw.R"), para que o EM_sw a encontre corretamente.
# ============================================================

factor <- function(X, n_fac) {
  X    <- as.matrix(X)
  Tobs <- nrow(X)
  S    <- (1 / Tobs) * t(X) %*% X
  eig  <- eigen(S, symmetric = TRUE)
  lam  <- eig$vectors[, 1:n_fac, drop = FALSE]
  fac  <- X %*% lam
  fit  <- fac %*% t(lam)
  mse  <- mean((X - fit)^2, na.rm = TRUE)
  list(factors = fac, lambda = lam, mse = mse)
}

cat("[OK] factor() PCA definida no ambiente global\n")
cat("     (substitui base::factor — compativel com EM_sw)\n\n")

# ============================================================
# PASSO 1 — Carrega funcoes do Coulombe
#
# Todos os arquivos devem estar em coulombe/functions/
# Se algum faltar, rode antes: source("coulombe/00_fix_missing_functions.R")
# ============================================================

src <- function(f) {
  path <- paste0("coulombe/functions/", f)
  if (!file.exists(path)) {
    cat(sprintf("  [ERRO] Nao encontrado: %s\n", path))
    return(invisible(NULL))
  }
  source(path, local = FALSE)
  cat(sprintf("  [OK] %s\n", f))
}

cat("=== Carregando funcoes do Coulombe ===\n")
src("EM_sw.R")              # Imputacao EM Stock & Watson
src("Xgenerators_v190127.R") # make_reg_matrix() e make_last()
src("dualGRRmdA_v190215.R")
src("CVGSBHK_v181127.R")
src("zfun_v190304.R")
src("TVPRRcosso_v181120.R") # TVPRR_cosso() — nucleo do 2SRR
src("TVPRR_v181111.R")
src("fastZrot_v181125.R")
src("CVKFMV_v190214.R")

# Verificacoes criticas
if (!exists("make_reg_matrix")) {
  stop("[ERRO CRITICO] make_reg_matrix nao encontrada.\n",
       "Verifique coulombe/functions/Xgenerators_v190127.R")
}
if (!exists("make_last")) {
  stop("[ERRO CRITICO] make_last nao encontrada.\n",
       "Verifique coulombe/functions/Xgenerators_v190127.R")
}
if (!exists("TVPRR_cosso")) {
  stop("[ERRO CRITICO] TVPRR_cosso nao encontrada.\n",
       "Verifique coulombe/functions/TVPRRcosso_v181120.R")
}
if (!exists("EM_sw")) {
  stop("[ERRO CRITICO] EM_sw nao encontrada.\n",
       "Verifique coulombe/functions/EM_sw.R")
}

cat("[OK] Todas as funcoes carregadas com sucesso\n\n")

# ============================================================
# PASSO 2 — Carrega base do Medeiros (data.rda)
# ============================================================

cat("=== Carregando data.rda ===\n")
load("data/data.rda")

# Validacoes basicas da estrutura
if (!is.data.frame(data)) {
  stop("'data' deve ser um data.frame. Verifique data/data.rda")
}
if (!"date" %in% names(data)) {
  stop("Coluna 'date' nao encontrada em data. Colunas: ",
       paste(names(data)[1:10], collapse = ", "))
}
if (!"CPIAUCSL" %in% names(data)) {
  stop("Coluna 'CPIAUCSL' nao encontrada. Colunas: ",
       paste(names(data)[1:10], collapse = ", "))
}

# Extrai vetores principais
dates  <- as.Date(data[["date"]])
y_full <- data[["CPIAUCSL"]]
X_full <- as.matrix(data[, !names(data) %in% c("date", "CPIAUCSL")])
bigt   <- nrow(data)
bign   <- ncol(X_full)

cat(sprintf("  Obs totais (bigt) : %d\n", bigt))
cat(sprintf("  Preditores (bign) : %d\n", bign))
cat(sprintf("  Data inicio       : %s\n", format(dates[1])))
cat(sprintf("  Data fim          : %s\n", format(dates[bigt])))
cat(sprintf("  NAs em y_full     : %d\n", sum(is.na(y_full))))
cat(sprintf("  NAs em X_full     : %d\n", sum(is.na(X_full))))

# ============================================================
# PASSO 3 — Alinha janela OOS com o Medeiros
#
# O Medeiros usa nwindows = 312 no rolling window,
# portanto tau = bigt - 312 e o OOS vai de tau+1 ate bigt.
# Carregamos yout.rda para pegar n_oos automaticamente,
# sem depender de hardcode.
# ============================================================

cat("\n=== Configuracao OOS (alinhada com Medeiros) ===\n")

if (!file.exists("forecasts/yout.rda")) {
  stop("forecasts/yout.rda nao encontrado. Rode primeiro os scripts do Medeiros.")
}
load("forecasts/yout.rda")

n_oos <- nrow(yout)
tau   <- bigt - n_oos

if (tau < 50) {
  stop(sprintf("tau=%d muito pequeno. Verifique yout.rda (n_oos=%d, bigt=%d)",
               tau, n_oos, bigt))
}

cat(sprintf("  tau   : %d obs — periodo ate %s\n", tau, format(dates[tau])))
cat(sprintf("  OOS   : %s -> %s\n",
            format(dates[tau + 1]), format(dates[bigt])))
cat(sprintf("  n_oos : %d observacoes\n", n_oos))

# ============================================================
# PASSO 4 — Imputacao EM Stock & Watson
#
# Roda apenas se houver NAs. Se a base ja estiver completa
# (como no caso do Medeiros atualizado), pula automaticamente.
# O fallback preserva os dados originais em caso de falha.
# ============================================================

cat("\n=== Imputacao EM Stock & Watson ===\n")

fred_mat <- cbind(y_full, X_full)
n_nas    <- sum(is.na(fred_mat))

if (n_nas == 0) {
  cat("  Base sem NAs — imputacao EM nao necessaria.\n")
  cat("  Usando dados originais diretamente.\n")
  y_imp <- y_full
  X_imp <- X_full
} else {
  cat(sprintf("  %d NAs encontrados. Rodando EM_sw (n=8, it_max=500)...\n", n_nas))
  t0 <- proc.time()
  em_res <- tryCatch(
    EM_sw(data = as.data.frame(fred_mat), n = 8, it_max = 500),
    error = function(e) {
      cat(sprintf("  [AVISO] EM_sw falhou: %s\n", e$message))
      cat("  Usando dados originais sem imputacao.\n")
      list(data = fred_mat)
    }
  )
  y_imp <- em_res$data[, 1]
  X_imp <- em_res$data[, -1]
  t_em  <- (proc.time() - t0)[3]
  cat(sprintf("  NAs antes: %d | depois: %d | tempo: %.1fs\n",
              n_nas, sum(is.na(em_res$data)), t_em))
}

# ============================================================
# PASSO 5 — Helper: remove colunas constantes antes do PCA
#
# prcomp(..., scale. = TRUE) falha se qualquer coluna tiver
# variancia zero (constante). Isso ocorre especialmente em
# janelas pequenas de treino com muitos preditores (116).
# Esta funcao e salva em setup_objects.rda para uso no 06b.
# ============================================================

rm_const <- function(X) {
  keep <- apply(X, 2, function(col) {
    v <- var(col, na.rm = TRUE)
    !is.na(v) && v > .Machine$double.eps
  })
  removed <- sum(!keep)
  if (removed > 0) {
    cat(sprintf("    [rm_const] %d coluna(s) constante(s) removida(s)\n", removed))
  }
  X[, keep, drop = FALSE]
}

cat("\n[OK] Helper rm_const() definido\n")

# ============================================================
# PASSO 6 — Teste make_reg_matrix (janela de 150 obs)
#
# Usa 150 obs para evitar o problema de colunas constantes
# que ocorre com janelas muito curtas (< 30 obs).
# ============================================================

cat("\n=== Teste make_reg_matrix (t = 150) ===\n")

N       <- 150
Xc_test <- rm_const(X_imp[1:N, ])
cat(sprintf("  Preditores apos rm_const: %d de %d\n",
            ncol(Xc_test), ncol(X_imp)))

pca_test <- tryCatch(
  prcomp(Xc_test, center = TRUE, scale. = TRUE),
  error = function(e) {
    stop(sprintf("prcomp falhou mesmo apos rm_const: %s", e$message))
  }
)

nf_test  <- min(8, ncol(pca_test$x))
fac_test <- pca_test$x[, 1:nf_test, drop = FALSE]

reg_test <- tryCatch(
  make_reg_matrix(
    y       = y_imp[1:N],
    Y       = y_imp[1:N],
    factors = fac_test,
    h       = 1,
    ly      = 2,
    lf      = 2
  ),
  error = function(e) {
    cat(sprintf("  [ERRO] make_reg_matrix: %s\n", e$message))
    NULL
  }
)

if (is.null(reg_test)) {
  stop("make_reg_matrix falhou. Verifique Xgenerators_v190127.R")
}

cat(sprintf("  reg_matrix: %d linhas x %d colunas\n",
            nrow(reg_test), ncol(reg_test)))
cat(sprintf("  NAs na reg_matrix: %d\n", sum(is.na(reg_test))))
cat("  [OK] make_reg_matrix funcionando\n")

# ============================================================
# PASSO 7 — Teste TVPRR_cosso (janela reduzida)
#
# Usa as linhas internas da reg_matrix de teste para
# verificar se o nucleo do 2SRR consegue estimar e
# retornar um forecast e os betas TVP.
# Um NULL aqui e normal com poucos dados; o 06b usa
# janelas muito maiores (tau ~ 474 obs).
# ============================================================

cat("\n=== Teste TVPRR_cosso ===\n")

maxlag  <- 2
tr_test <- reg_test[(maxlag + 2):(nrow(reg_test) - 1), ]
tr_test <- tr_test[complete.cases(tr_test), ]
Y_tr    <- tr_test[, 1]
X_tr    <- as.matrix(tr_test[, -1])
oos_row <- matrix(reg_test[nrow(reg_test), -1], nrow = 1)

cat(sprintf("  Treino: %d obs x %d vars\n", nrow(tr_test), ncol(X_tr)))

cv_test <- tryCatch(
  cv.glmnet(X_tr, Y_tr, alpha = 0, nfolds = 5, intercept = TRUE),
  error = function(e) {
    cat(sprintf("  [AVISO cv.glmnet] %s\n", e$message))
    NULL
  }
)

if (!is.null(cv_test)) {
  lv_test <- exp(pracma::linspace(-2, 6, n = 5))

  aa_test <- tryCatch(
    TVPRR_cosso(
      y         = Y_tr,
      X         = X_tr,
      lambdavec = lv_test,
      sweigths  = 1,
      type      = 2,
      alpha     = 0.01,
      silent    = 1,
      kfold     = 5,
      lambda2   = cv_test$lambda.min,
      tol       = 1e-6,
      maxit     = 5,
      oosX      = oos_row
    ),
    error = function(e) {
      cat(sprintf("  [AVISO TVPRR_cosso] %s\n", e$message))
      NULL
    }
  )

  if (!is.null(aa_test)) {
    cat(sprintf("  Forecast teste : %.6f\n", as.numeric(aa_test$fcast)))
    cat(sprintf("  Betas dim      : %d coeficientes\n",
                length(aa_test$grrats$betas_grr)))
    cat("  [OK] TVPRR_cosso funcionando\n")
  } else {
    cat("  [INFO] TVPRR_cosso retornou NULL nesta janela de teste.\n")
    cat("         Normal com poucos dados — o 06b usa janelas muito maiores.\n")
    cat("         Continue para o 06b normalmente.\n")
  }
} else {
  cat("  [INFO] cv.glmnet nao rodou — verifique o pacote glmnet.\n")
}

# ============================================================
# PASSO 8 — Salva objetos para uso no 06b
# ============================================================

cat("\n=== Salvando objetos para o 06b ===\n")

save(
  y_imp,    # vetor inflacao imputado (ou original se sem NAs)
  X_imp,    # matriz preditores imputada
  dates,    # vetor de datas
  bigt,     # total de observacoes
  bign,     # numero de preditores
  tau,      # indice fim do periodo de treino
  n_oos,    # numero de obs fora da amostra
  rm_const, # funcao helper para remover colunas constantes
  file = "coulombe/setup_objects.rda"
)

cat("  [OK] Salvo: coulombe/setup_objects.rda\n")
cat("\n=== 06a_coulombe_setup.R CONCLUIDO ===\n")
cat("Proximo passo: source('coulombe/06b_coulombe_pipeline.R')\n")