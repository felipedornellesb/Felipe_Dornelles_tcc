# =============================================================================================================
# SCRIPT 03: LOCAL PROJECTIONS — IRF VARIANDO NO TEMPO (TVP via 2SRR)
# =============================================================================================================
# Jordà (2005): y_{t+h} = α_h + β_{h,t} * X_t + ε_{t+h}
# Estimador:    tvp.ridge() de MV2SRR_v221103.R
# Choque:       EMBI+ (JPM366_EMBI366) → resposta em IPCA e EMBI+
# Benchmarks:   OLS estático, Ridge homogêneo (1SRR), TVP-Ridge (2SRR)
# Testes:       Diebold-Mariano, estabilidade de Chow (breakpoints), R² OOS
# =============================================================================================================

rm(list = ls())

wd <- 'C:/Users/felip/OneDrive/Documentos/tcc/Felipe_Dornelles_tcc'
setwd(wd)

paths <- list(
  program   = "00_program",
  data      = "10_data",
  tools     = "20_tools",
  functions = "20_tools/functions",
  output    = "30_output",
  results   = file.path("40_results", "40_analysis",
                        paste0("results_", format(Sys.Date(), "%d_%m_%Y_%H")))
)

# Cria todas as pastas da hierarquia automaticamente
dir.create(paths$results, showWarnings = FALSE, recursive = TRUE)
dir.create(paths$output,  showWarnings = FALSE, recursive = TRUE)

cat(sprintf("Outputs salvos em: %s\n", paths$results))
# =============================================================================================================
# PACOTES
# =============================================================================================================

myPKGs <- c('dplyr', 'glmnet', 'fGarch', 'matrixcalc', 'pracma',
            'ggplot2', 'tidyr', 'reshape2', 'scales', 'patchwork',
            'forecast',   # dm.test
            'strucchange' # breakpoints (Chow)
            )

InstalledPKGs    <- names(installed.packages()[, 'Package'])
InstallThesePKGs <- myPKGs[!myPKGs %in% InstalledPKGs]
if (length(InstallThesePKGs) > 0)
  install.packages(InstallThesePKGs, repos = "http://cran.us.r-project.org")
invisible(lapply(myPKGs, library, character.only = TRUE))

# =============================================================================================================
# FUNÇÕES
# =============================================================================================================

source(paste(paths$tools,     'EM_sw.R',               sep = '/'))
source(paste(paths$tools,     'ICp2.R',                sep = '/'))
source(paste(paths$functions, 'MV2SRR_v221103.R',      sep = '/'))

# =============================================================================================================
# UTILITÁRIOS INTERNOS
# =============================================================================================================

# Salva tabela como CSV e imprime no console com separador visual
save_and_print <- function(df, name, paths, digits = 4) {
  path_csv <- file.path(paths$results, paste0(name, ".csv"))
  write.csv(df, path_csv, row.names = FALSE)
  cat(sprintf("\n╔══ TABELA: %s ══╗\n", name))
  print(df, digits = digits)
  cat(sprintf("╚══ Salvo: %s ══╝\n", path_csv))
  invisible(df)
}

# Tema ggplot padronizado para o TCC
theme_tcc <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title       = element_text(face = "bold", size = base_size + 1),
      plot.subtitle    = element_text(size = base_size - 1, color = "gray40"),
      panel.grid.minor = element_blank(),
      legend.position  = "bottom",
      axis.text.x      = element_text(angle = 30, hjust = 1),
      strip.text       = element_text(face = "bold")
    )
}

# R² fora da amostra (Campbell & Thompson 2008)
oos_r2 <- function(actual, pred, benchmark) {
  msfe_model <- mean((actual - pred)^2,      na.rm = TRUE)
  msfe_bench <- mean((actual - benchmark)^2, na.rm = TRUE)
  1 - msfe_model / msfe_bench
}

# =============================================================================================================
# DADOS
# =============================================================================================================

load(paste(paths$data, "df.rda", sep = '/'))

df_means <- colMeans(df, na.rm = TRUE)
df_sds   <- apply(df, 2, sd, na.rm = TRUE)
df       <- as.matrix(scale(df))

T_total <- nrow(df)
n_vars  <- ncol(df)

# Datas mensais — ajuste se necessário
dates_all <- seq(as.Date("1996-01-01"), by = "month", length.out = T_total)

cat(sprintf("\nDados: %d obs x %d vars | %s a %s\n",
            T_total, n_vars,
            format(min(dates_all), "%b/%Y"),
            format(max(dates_all), "%b/%Y")))

# =============================================================================================================
# CONFIGURAÇÕES
# =============================================================================================================

targets    <- c("PRECOS12_IPCA12", "JPM366_EMBI366")
labels_pt  <- c("IPCA (Inflação 12m)", "EMBI+ (Risco País)")
shock_var  <- "JPM366_EMBI366"

# Filtra variáveis existentes
existe    <- targets %in% colnames(df)
targets   <- targets[existe]
labels_pt <- labels_pt[existe]

shock_pos <- which(colnames(df) == shock_var)
if (length(shock_pos) == 0) stop("Variável de choque não encontrada!")

cat(sprintf("Choque: %s (coluna %d) | Targets: %s\n",
            shock_var, shock_pos, paste(targets, collapse = ", ")))

H           <- 12
lambdavec   <- exp(pracma::linspace(4, 18, n = 15))

# Períodos históricos para cortes
periodos_ref <- list(
  "Crise 2002"   = "2002-10",
  "Crise 2008"   = "2008-09",
  "Recessão 2015" = "2015-01",
  "COVID-19"     = "2020-04"
)

# =============================================================================================================
# LOOP PRINCIPAL POR VARIÁVEL-ALVO
# =============================================================================================================

# Acumula resultados de todos os targets para tabela comparativa final
all_fit_stats   <- list()
all_dm_lp       <- list()

for (v in seq_along(targets)) {

  target_name <- targets[v]
  label       <- labels_pt[v]

  cat(sprintf("\n%s\n[LP-TVP] %s\n%s\n",
              strrep("=", 60), label, strrep("=", 60)))

  # -----------------------------------------------------------------------
  # 1. MATRIZES DE PROJEÇÃO LOCAL
  # -----------------------------------------------------------------------

  T_eff <- T_total - H
  Xmat  <- df[1:T_eff, ]
  Y_vec <- df[, target_name]
  Ymat  <- matrix(NA, nrow = T_eff, ncol = H)
  for (h in 1:H) Ymat[, h] <- Y_vec[(1 + h):(T_eff + h)]

  # -----------------------------------------------------------------------
  # 2. LAMBDA BASE VIA CV (glmnet)
  # -----------------------------------------------------------------------

  cat("  [1/6] Lambda base via CV glmnet...\n")
  cvvec <- numeric(H)
  for (h in 1:H) {
    cv_fit   <- cv.glmnet(Xmat, Ymat[, h], family = "gaussian", alpha = 0)
    cvvec[h] <- cv_fit$lambda.min
  }
  lambda2_base <- mean(cvvec) / 2
  cat(sprintf("        lambda2_base = %.6f\n", lambda2_base))

  # -----------------------------------------------------------------------
  # 3. BENCHMARKS: OLS e RIDGE ESTÁTICO (1SRR)
  # -----------------------------------------------------------------------

  cat("  [2/6] Estimando benchmarks OLS e Ridge estático...\n")

  yhat_ols   <- matrix(NA, T_eff, H)
  yhat_ridge <- matrix(NA, T_eff, H)

  for (h in 1:H) {
    # OLS (sem penalização)
    ols_fit        <- lm(Ymat[, h] ~ Xmat)
    yhat_ols[, h]  <- fitted(ols_fit)

    # Ridge estático (glmnet, alpha=0)
    ridge_fit         <- glmnet(Xmat, Ymat[, h], alpha = 0,
                                lambda = cvvec[h])
    yhat_ridge[, h]   <- as.vector(predict(ridge_fit, newx = Xmat))
  }

  # -----------------------------------------------------------------------
  # 4. TVP-RIDGE 2SRR (modelo principal)
  # -----------------------------------------------------------------------

  cat("  [3/6] Estimando TVP-Ridge (2SRR) para h = 1...", H, "...\n")

  betas_tvp_list <- vector("list", H)
  yhat_tvp       <- matrix(NA, T_eff, H)

  for (h in 1:H) {
    cat(sprintf("        h = %2d\n", h))
    out_h <- tryCatch({
      tvp.ridge(
        X                 = Xmat,
        Y                 = as.matrix(Ymat[, h]),
        lambda.candidates = lambdavec,
        oosX              = c(),
        lambda2           = lambda2_base,
        kfold             = 5,
        CV.plot           = FALSE,
        CV.2SRR           = TRUE,
        block_size        = 8,
        sig.u.param       = 0.75,
        sig.eps.param     = 0.75,
        ols.prior         = 0
      )
    }, error = function(e) {
      cat(sprintf("        [ERRO] h=%d: %s\n", h, e$message))
      NULL
    })
    betas_tvp_list[[h]] <- out_h
    if (!is.null(out_h)) yhat_tvp[, h] <- as.vector(out_h$yhat.2srr)
  }

  # T_interno (pode ser < T_eff por causa do Zfun)
  T_interno <- NULL
  for (h in 1:H) {
    if (!is.null(betas_tvp_list[[h]])) {
      T_interno <- dim(betas_tvp_list[[h]]$betas.2srr)[3]
      break
    }
  }
  if (is.null(T_interno)) {
    cat(sprintf("  [ERRO FATAL] Todos os horizontes falharam para %s.\n", target_name))
    next
  }

  # -----------------------------------------------------------------------
  # 5. TABELA DE FIT IN-SAMPLE: R², MSFE, MAE por horizonte e modelo
  # -----------------------------------------------------------------------

  cat("  [4/6] Calculando estatísticas de ajuste in-sample...\n")

  fit_stats <- do.call(rbind, lapply(1:H, function(h) {
    y_h   <- Ymat[, h]
    bench <- mean(y_h)   # média histórica como benchmark naive

    calc_row <- function(yhat, modelo) {
      resid  <- y_h - yhat
      r2     <- 1 - sum(resid^2, na.rm = TRUE) / sum((y_h - mean(y_h))^2, na.rm = TRUE)
      data.frame(
        variable = target_name,
        modelo   = modelo,
        horizonte = h,
        R2       = round(r2, 4),
        MSFE     = round(mean(resid^2, na.rm = TRUE), 4),
        MAE      = round(mean(abs(resid), na.rm = TRUE), 4),
        OOS_R2   = round(oos_r2(y_h, yhat, rep(bench, length(y_h))), 4)
      )
    }

    rbind(
      calc_row(yhat_ols[, h],   "OLS"),
      calc_row(yhat_ridge[, h], "Ridge-1SRR"),
      calc_row(yhat_tvp[, h],   "TVP-2SRR")
    )
  }))

  save_and_print(fit_stats,
                 sprintf("TAB01_fit_inSample_%s", target_name), paths)
  all_fit_stats[[target_name]] <- fit_stats

  # -----------------------------------------------------------------------
  # 6. MATRIZES DE IRF
  # -----------------------------------------------------------------------

  IRF_raw <- matrix(NA, H, T_interno)
  IRF_cum <- matrix(NA, H, T_interno)

  for (h in 1:H) {
    if (!is.null(betas_tvp_list[[h]])) {
      beta_h        <- betas_tvp_list[[h]]$betas.2srr[1, shock_pos + 1, ]
      len           <- min(length(beta_h), T_interno)
      IRF_raw[h, 1:len] <- beta_h[1:len]
    }
  }
  for (t in 1:T_interno) IRF_cum[, t] <- cumsum(IRF_raw[, t])

  dates_irf <- dates_all[1:T_interno]

  # -----------------------------------------------------------------------
  # 7. TABELA IRF CUMULATIVA MÉDIA POR HORIZONTE
  # -----------------------------------------------------------------------

  irf_summary <- data.frame(
    variable  = target_name,
    horizonte = 1:H,
    media     = round(rowMeans(IRF_cum, na.rm = TRUE), 4),
    mediana   = round(apply(IRF_cum, 1, median, na.rm = TRUE), 4),
    dp        = round(apply(IRF_cum, 1, sd,     na.rm = TRUE), 4),
    p10       = round(apply(IRF_cum, 1, quantile, 0.10, na.rm = TRUE), 4),
    p25       = round(apply(IRF_cum, 1, quantile, 0.25, na.rm = TRUE), 4),
    p75       = round(apply(IRF_cum, 1, quantile, 0.75, na.rm = TRUE), 4),
    p90       = round(apply(IRF_cum, 1, quantile, 0.90, na.rm = TRUE), 4),
    pct_pos   = round(apply(IRF_cum, 1, function(x) mean(x > 0, na.rm = TRUE)), 4)
  )

  save_and_print(irf_summary,
                 sprintf("TAB02_IRF_summary_%s", target_name), paths)

  # -----------------------------------------------------------------------
  # 8. TABELA IRF EM PERÍODOS HISTÓRICOS
  # -----------------------------------------------------------------------

  periodos_ativos <- list()
  for (nome in names(periodos_ref)) {
    idx <- which(format(dates_irf, "%Y-%m") == periodos_ref[[nome]])
    if (length(idx) > 0) periodos_ativos[[nome]] <- idx[1]
  }

  if (length(periodos_ativos) >= 1) {
    irf_periodos <- do.call(rbind, lapply(names(periodos_ativos), function(nome) {
      t_idx <- periodos_ativos[[nome]]
      data.frame(
        variable  = target_name,
        periodo   = nome,
        data      = format(dates_irf[t_idx], "%b/%Y"),
        stringsAsFactors = FALSE,
        as.data.frame(t(round(IRF_cum[, t_idx], 4)))
      )
    }))
    colnames(irf_periodos)[5:ncol(irf_periodos)] <- paste0("h", 1:H)
    save_and_print(irf_periodos,
                   sprintf("TAB03_IRF_periodos_%s", target_name), paths)
  }

  # -----------------------------------------------------------------------
  # 9. TESTE DIEBOLD-MARIANO: TVP-2SRR vs OLS e vs Ridge-1SRR
  # -----------------------------------------------------------------------

  cat("  [5/6] Testes Diebold-Mariano...\n")

  dm_lp <- do.call(rbind, lapply(1:H, function(h) {
    y_h     <- Ymat[, h]
    e_tvp   <- y_h - yhat_tvp[, h]
    e_ols   <- y_h - yhat_ols[, h]
    e_ridge <- y_h - yhat_ridge[, h]

    dm_ols   <- tryCatch(dm.test(e_tvp, e_ols,   alternative = "less", h = h, power = 2),
                         error = function(e) NULL)
    dm_ridge <- tryCatch(dm.test(e_tvp, e_ridge, alternative = "less", h = h, power = 2),
                         error = function(e) NULL)

    data.frame(
      variable       = target_name,
      horizonte      = h,
      # vs OLS
      MSFE_rel_OLS   = round(mean(e_tvp^2, na.rm=TRUE) / mean(e_ols^2,   na.rm=TRUE), 4),
      DM_vs_OLS      = if (!is.null(dm_ols))   round(dm_ols$statistic,   4) else NA,
      p_vs_OLS       = if (!is.null(dm_ols))   round(dm_ols$p.value,     4) else NA,
      sig_vs_OLS     = if (!is.null(dm_ols))   ifelse(dm_ols$p.value < 0.01, "***",
                                               ifelse(dm_ols$p.value < 0.05, "**",
                                               ifelse(dm_ols$p.value < 0.10, "*", ""))) else "",
      # vs Ridge-1SRR
      MSFE_rel_Ridge = round(mean(e_tvp^2, na.rm=TRUE) / mean(e_ridge^2, na.rm=TRUE), 4),
      DM_vs_Ridge    = if (!is.null(dm_ridge)) round(dm_ridge$statistic, 4) else NA,
      p_vs_Ridge     = if (!is.null(dm_ridge)) round(dm_ridge$p.value,   4) else NA,
      sig_vs_Ridge   = if (!is.null(dm_ridge)) ifelse(dm_ridge$p.value < 0.01, "***",
                                               ifelse(dm_ridge$p.value < 0.05, "**",
                                               ifelse(dm_ridge$p.value < 0.10, "*", ""))) else ""
    )
  }))

  save_and_print(dm_lp,
                 sprintf("TAB04_DM_LP_%s", target_name), paths)
  all_dm_lp[[target_name]] <- dm_lp

  # -----------------------------------------------------------------------
  # 10. TESTE DE ESTABILIDADE DO CHOQUE (Breakpoints via strucchange)
  # -----------------------------------------------------------------------

  cat("  [6/6] Teste de estabilidade temporal (breakpoints)...\n")

  bp_results <- do.call(rbind, lapply(c(1, 3, 6, 12), function(h) {
    if (h > H) return(NULL)
    beta_h <- IRF_raw[h, ]
    beta_h <- beta_h[!is.na(beta_h)]
    if (length(beta_h) < 20) return(NULL)

    bp <- tryCatch({
      breakpoints(beta_h ~ 1, h = 0.15)
    }, error = function(e) NULL)

    n_breaks <- if (!is.null(bp) && !all(is.na(bp$breakpoints))) {
      sum(!is.na(bp$breakpoints))
    } else 0

    bp_dates <- if (n_breaks > 0) {
      paste(format(dates_irf[bp$breakpoints[!is.na(bp$breakpoints)]], "%b/%Y"),
            collapse = "; ")
    } else "—"

    data.frame(
      variable    = target_name,
      horizonte   = h,
      n_breaks    = n_breaks,
      break_dates = bp_dates,
      stringsAsFactors = FALSE
    )
  }))

  if (!is.null(bp_results) && nrow(bp_results) > 0)
    save_and_print(bp_results,
                   sprintf("TAB05_breakpoints_%s", target_name), paths)

  # -----------------------------------------------------------------------
  # 11. GRÁFICOS
  # -----------------------------------------------------------------------

  cat("\n  Gerando gráficos...\n")

  horizons_sel <- c(1, 3, 6, 12)

  # --- FIG 1: IRF cumulativa no tempo ---
  plot_df1 <- data.frame(
    date      = rep(dates_irf, length(horizons_sel)),
    irf       = as.vector(t(IRF_cum[horizons_sel, ])),
    horizonte = factor(rep(paste0("h = ", horizons_sel, "m"), each = T_interno),
                       levels = paste0("h = ", sort(horizons_sel), "m"))
  )

  p1 <- ggplot(plot_df1, aes(x = date, y = irf, color = horizonte)) +
    geom_line(linewidth = 0.9) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    scale_x_date(labels = date_format("%Y"), date_breaks = "3 years") +
    scale_color_brewer(palette = "Set1") +
    labs(title    = sprintf("IRF-TVP: %s ← Choque EMBI+", label),
         subtitle = "IRF cumulativa — Projeções Locais com TVP-Ridge (2SRR)",
         x = NULL, y = "Resposta acumulada (dp padronizados)", color = "Horizonte") +
    theme_tcc()

  # --- FIG 2: IRF por período histórico ---
  if (length(periodos_ativos) >= 2) {
    irf_hist <- do.call(rbind, lapply(names(periodos_ativos), function(nome) {
      t_idx <- periodos_ativos[[nome]]
      data.frame(horizonte = 1:H, irf_cum = IRF_cum[, t_idx], periodo = nome)
    }))

    p2 <- ggplot(irf_hist, aes(x = horizonte, y = irf_cum,
                               color = periodo, linetype = periodo)) +
      geom_line(linewidth = 1.1) +
      geom_point(size = 2.5) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      scale_x_continuous(breaks = 1:H) +
      scale_color_brewer(palette = "Dark2") +
      labs(title    = sprintf("IRF em Crises: %s ← Choque EMBI+", label),
           subtitle = "Cada linha = curva IRF cumulativa estimada no período",
           x = "Horizonte h (meses)", y = "Resposta acumulada (dp)",
           color = "Período", linetype = "Período") +
      theme_tcc()
  }

  # --- FIG 3: Heatmap IRF-TVP ---
  irf_melt       <- reshape2::melt(IRF_cum)
  colnames(irf_melt) <- c("horizonte", "t_idx", "irf")
  irf_melt$data  <- dates_irf[irf_melt$t_idx]
  lim_abs        <- ceiling(max(abs(irf_melt$irf), na.rm = TRUE) * 10) / 10

  p3 <- ggplot(irf_melt, aes(x = data, y = factor(horizonte, levels = rev(1:H)), fill = irf)) +
    geom_tile() +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#d6604d",
                         midpoint = 0, limits = c(-lim_abs, lim_abs),
                         name = expression(beta[h*","*t])) +
    scale_x_date(labels = date_format("%Y"), date_breaks = "3 years") +
    labs(title    = sprintf("Heatmap IRF-TVP: %s ← Choque EMBI+", label),
         subtitle = "Vermelho = resposta positiva | Azul = negativa",
         x = NULL, y = "Horizonte h (meses)") +
    theme_tcc() + theme(legend.position = "right", panel.grid = element_blank())

  # --- FIG 4: Betas não-cumulativos ---
  plot_df4 <- data.frame(
    date      = rep(dates_irf, length(horizons_sel)),
    beta      = as.vector(t(IRF_raw[horizons_sel, ])),
    horizonte = factor(rep(paste0("h = ", horizons_sel, "m"), each = T_interno),
                       levels = paste0("h = ", sort(horizons_sel), "m"))
  )

  p4 <- ggplot(plot_df4, aes(x = date, y = beta, color = horizonte)) +
    geom_line(linewidth = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    scale_x_date(labels = date_format("%Y"), date_breaks = "3 years") +
    scale_color_brewer(palette = "Set2") +
    labs(title    = sprintf("Beta Bruto (não-cumulativo): %s ← Choque EMBI+", label),
         subtitle = expression(paste(beta[h*","*t], " = resposta marginal")),
         x = NULL, y = expression(beta[h*","*t]), color = "Horizonte") +
    theme_tcc()

  # --- FIG 5: IRF mediana com dispersão ---
  p5 <- ggplot(irf_summary, aes(x = horizonte)) +
    geom_ribbon(aes(ymin = p10, ymax = p90), fill = "#4575b4", alpha = 0.15) +
    geom_ribbon(aes(ymin = p25, ymax = p75), fill = "#4575b4", alpha = 0.30) +
    geom_line(aes(y = mediana), color = "#2c3e7f", linewidth = 1.3) +
    geom_point(aes(y = mediana), color = "#2c3e7f", size = 2.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    scale_x_continuous(breaks = 1:H) +
    labs(title    = sprintf("IRF Mediana Histórica: %s ← Choque EMBI+", label),
         subtitle = "Linha = mediana | Faixas = IQR (50%) e P10-P90",
         x = "Horizonte h (meses)", y = "Resposta acumulada mediana (dp)") +
    theme_tcc() + theme(legend.position = "none")

  # --- FIG 6: R² in-sample por modelo e horizonte ---
  r2_df <- fit_stats %>%
    select(modelo, horizonte, R2) %>%
    mutate(modelo = factor(modelo, levels = c("OLS", "Ridge-1SRR", "TVP-2SRR")))

  p6 <- ggplot(r2_df, aes(x = horizonte, y = R2, color = modelo, linetype = modelo)) +
    geom_line(linewidth = 1.0) +
    geom_point(size = 2.5) +
    scale_x_continuous(breaks = 1:H) +
    scale_color_manual(values = c("OLS" = "#e41a1c",
                                  "Ridge-1SRR" = "#ff7f00",
                                  "TVP-2SRR"  = "#377eb8")) +
    labs(title    = sprintf("R² In-Sample por Modelo: %s", label),
         subtitle = "OLS vs Ridge-1SRR vs TVP-Ridge-2SRR",
         x = "Horizonte h (meses)", y = expression(R^2),
         color = "Modelo", linetype = "Modelo") +
    theme_tcc()

  # --- FIG 7: MSFE relativo OLS e Ridge (TVP = base 1) ---
  msfe_rel_df <- dm_lp %>%
    select(horizonte, MSFE_rel_OLS, MSFE_rel_Ridge) %>%
    tidyr::pivot_longer(-horizonte, names_to = "vs", values_to = "MSFE_rel") %>%
    mutate(vs = recode(vs,
                       "MSFE_rel_OLS"   = "vs OLS",
                       "MSFE_rel_Ridge" = "vs Ridge-1SRR"))

  p7 <- ggplot(msfe_rel_df, aes(x = horizonte, y = MSFE_rel, color = vs, linetype = vs)) +
    geom_line(linewidth = 1.0) +
    geom_point(size = 2.5) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
    scale_x_continuous(breaks = 1:H) +
    scale_color_manual(values = c("vs OLS" = "#e41a1c", "vs Ridge-1SRR" = "#ff7f00")) +
    annotate("text", x = H * 0.9, y = 1.02,
             label = "MSFE igual ao benchmark", size = 3, color = "gray40") +
    labs(title    = sprintf("MSFE Relativo — TVP-2SRR vs Benchmarks: %s", label),
         subtitle = "< 1 = TVP-2SRR superior | Teste DM na TAB04",
         x = "Horizonte h (meses)", y = "MSFE relativo (TVP÷Benchmark)",
         color = NULL, linetype = NULL) +
    theme_tcc()

  # --- FIG 8: Painel 2x2 resumo ---
  p_painel <- (p5 | p6) / (p7 | p3) +
    plot_annotation(
      title    = sprintf("Painel Resumo LP-TVP: %s ← Choque EMBI+", label),
      subtitle = "IRF mediana | R² por modelo | MSFE relativo | Heatmap",
      theme    = theme(plot.title    = element_text(face = "bold", size = 14),
                       plot.subtitle = element_text(size = 11, color = "gray40"))
    )

  # --- Salvar todos os gráficos ---
  figs <- list(
    list(p = p1,       name = "fig1_IRF_tempo"),
    list(p = p3,       name = "fig3_heatmap"),
    list(p = p4,       name = "fig4_beta_bruto"),
    list(p = p5,       name = "fig5_IRF_mediana"),
    list(p = p6,       name = "fig6_R2_modelos"),
    list(p = p7,       name = "fig7_MSFE_relativo"),
    list(p = p_painel, name = "fig8_painel_resumo")
  )
  if (length(periodos_ativos) >= 2)
    figs <- c(figs, list(list(p = p2, name = "fig2_IRF_periodos")))

  for (fig in figs) {
    for (ext in c("pdf", "png")) {
      fname <- file.path(paths$results,
                         sprintf("%s_%s.%s", fig$name, target_name, ext))
      ggsave(fname, fig$p,
             width  = if (grepl("painel", fig$name)) 14 else 11,
             height = if (grepl("painel", fig$name)) 10 else 5.5,
             dpi    = 300)
    }
    cat(sprintf("    -> %s salvo.\n", fig$name))
  }

  # Salva RDA com tudo
  rda_file <- file.path(paths$output, sprintf("LP_IRF_%s.rda", target_name))
  save(betas_tvp_list, IRF_raw, IRF_cum, irf_summary,
       yhat_ols, yhat_ridge, yhat_tvp,
       fit_stats, dm_lp,
       dates_irf, Ymat, Xmat,
       target_name, label, shock_var, shock_pos,
       lambda2_base, H, cvvec,
       file = rda_file)
  cat(sprintf("  -> RDA salvo: %s\n", rda_file))
}

# =============================================================================================================
# TABELAS COMPARATIVAS FINAIS (todos os targets)
# =============================================================================================================

if (length(all_fit_stats) > 1) {

  # Tabela de R² comparativo
  r2_comp <- do.call(rbind, all_fit_stats) %>%
    filter(modelo == "TVP-2SRR") %>%
    select(variable, horizonte, R2, MSFE, MAE)

  save_and_print(r2_comp, "TAB06_R2_comparativo_targets", paths)

  # Tabela DM comparativo
  dm_comp <- do.call(rbind, all_dm_lp) %>%
    select(variable, horizonte,
           MSFE_rel_OLS,   sig_vs_OLS,
           MSFE_rel_Ridge, sig_vs_Ridge)

  save_and_print(dm_comp, "TAB07_DM_comparativo_targets", paths)
}

# =============================================================================================================
# RESUMO FINAL NO CONSOLE
# =============================================================================================================

cat(sprintf("\n%s\nScript 03 concluído!\n%s\n", strrep("=",60), strrep("=",60)))
cat(sprintf("Resultados em: %s\n\n", paths$results))
cat("TABELAS (CSV):\n")
cat("  TAB01 — R², MSFE, MAE in-sample por modelo e horizonte\n")
cat("  TAB02 — IRF cumulativa: média, mediana, DP, percentis, % positiva\n")
cat("  TAB03 — IRF em períodos históricos (crises)\n")
cat("  TAB04 — DM test: TVP-2SRR vs OLS e vs Ridge-1SRR\n")
cat("  TAB05 — Breakpoints (estabilidade temporal do choque)\n")
cat("  TAB06 — R² comparativo entre targets (se > 1 variável)\n")
cat("  TAB07 — DM comparativo entre targets\n\n")
cat("GRÁFICOS (PDF + PNG):\n")
cat("  fig1 — IRF cumulativa no tempo (h=1,3,6,12)\n")
cat("  fig2 — IRF em períodos de crise histórica\n")
cat("  fig3 — Heatmap IRF-TVP (horizonte × tempo)\n")
cat("  fig4 — Betas não-cumulativos por horizonte\n")
cat("  fig5 — IRF mediana com dispersão temporal\n")
cat("  fig6 — R² in-sample: OLS vs Ridge vs TVP\n")
cat("  fig7 — MSFE relativo TVP vs benchmarks\n")
cat("  fig8 — Painel resumo 2×2\n")
cat(strrep("=", 60), "\n")