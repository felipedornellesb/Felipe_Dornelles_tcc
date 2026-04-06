# =============================================================================================================
# SCRIPT 03: LOCAL PROJECTIONS (FUNÇÕES DE RESPOSTA AO IMPULSO VARIANDO NO TEMPO)
# =============================================================================================================
#
# OBJETIVO DO SCRIPT (Para discussão com o orientador):
# 
# 1. Previsão vs. Análise Estrutural: 
#    O script anterior (02_forecast) tinha um objetivo puramente preditivo: minimizar 
#    o erro fora da amostra (RMSPE). Este script (03) tem um objetivo estrutural: 
#    entender a dinâmica da economia através de Funções de Resposta ao Impulso (IRFs).
#
# 2. Metodologia de Projeções Locais (Jordà, 2005):
#    Em vez de estimar um VAR gigante e iterar os coeficientes (o que acumula erros),
#    as Projeções Locais estimam uma regressão direta para cada horizonte h:
#    $y_{t+h} = \alpha_h + \beta_h X_t + \epsilon_{t+h}$
#    O coeficiente $\beta_h$ é exatamente a resposta ao impulso no horizonte h.
#
# 3. Inovação do Coulombe (TVP-Ridge):
#    Aqui, aplicamos a Ridge Regression com Parâmetros Variando no Tempo (TVP). 
#    Isso significa que o nosso $\beta_h$ ganha um subscrito t ($\beta_{h,t}$). 
#    Isso nos permite ver, por exemplo, se a inflação responde de forma diferente a um 
#    choque no Risco País hoje do que respondia na década de 1990.
#
# 4. Adaptação para o TCC:
#    No artigo original, Coulombe tinha uma série já isolada de "choque de política monetária". 
#    Como não temos essa série estrutural limpa aqui, assumimos o SPREAD (EMBI+) 
#    como a nossa variável de "choque" (inovação exógena), para analisar como o IPCA 
#    e o Desemprego reagem a ele ao longo de 12 meses.
#
# =============================================================================================================

rm(list = ls()) # Limpa o ambiente

# Define o diretório de trabalho
wd = 'C:/Users/felip/OneDrive/Documentos/tcc/Felipe_Dornelles_tcc'
setwd(wd)

paths <- list(program = "00_program",
              data = "10_data",
              tools = "20_tools",
              functions = "20_tools/functions",
              output = "30_output",
              results = "40_results")

# =============================================================================================================
# GESTÃO DE PACOTES
# =============================================================================================================

myPKGs <- c('dplyr', 'glmnet', 'fGarch', 'matrixcalc', 'pracma', 'e1071', 'GA')

InstalledPKGs <- names(installed.packages()[,'Package'])
InstallThesePKGs <- myPKGs[!myPKGs %in% InstalledPKGs]
if (length(InstallThesePKGs) > 0) install.packages(InstallThesePKGs, repos = "http://cran.us.r-project.org")

invisible(lapply(myPKGs, library, character.only = TRUE))

# =============================================================================================================
# CARREGAMENTO DE FUNÇÕES E FERRAMENTAS DO COULOMBE
# =============================================================================================================

source(paste(paths$functions, '00_Nathalia_functions.R', sep='/'))
source(paste(paths$tools, 'EM_sw.R', sep='/'))
source(paste(paths$tools, 'ICp2.R', sep='/'))
source(paste(paths$tools, 'Xgenerators_v190127.R', sep='/'))
source(paste(paths$functions, 'dualGRRmdA_v190215.R', sep='/'))
source(paste(paths$functions, 'CVGSBHK_v181127.R', sep='/'))
source(paste(paths$functions, 'zfun_v190304.R', sep='/'))
source(paste(paths$functions, 'factor.R', sep='/'))
source(paste(paths$functions, 'TVPRRcosso_v181120.R', sep='/'))

# =============================================================================================================
# PREPARAÇÃO DE DADOS PARA AS PROJEÇÕES LOCAIS
# =============================================================================================================

load(paste(paths$data, "df.rda", sep='/'))
df = as.matrix(scale(df)) # Padronização essencial para a Ridge Regression

# Definimos as variáveis que queremos observar a reação (Targets)
targets <- c("PRECOS12_IPCA12", "SEADE12_TDTGSP12", "JPM366_EMBI366")

# Definimos a variável que vai sofrer o "Choque" no tempo t
shock_var <- "JPM366_EMBI366" 
shock_pos <- which(colnames(df) == shock_var) # Pega o índice da coluna do choque na matriz X

H <- 12 # Quantos meses à frente queremos analisar o impacto do choque? (Horizonte da IRF)

# Grid de busca para a penalidade do modelo (o quanto os parâmetros variam no tempo)
lambdavec = exp(pracma::linspace(4, 18, n=15))

# Criamos um "cubo" (array 3D) para guardar a Função de Resposta ao Impulso.
# Dimensões: [Horizonte, Variável Alvo, Tempo (pois a resposta varia no tempo!)]
IRF_TVP <- array(NA, dim=c(H, length(targets), nrow(df) - H))

# =============================================================================================================
# LOOP DE ESTIMAÇÃO (PROJEÇÕES LOCAIS)
# =============================================================================================================

for(v in 1:length(targets)) {
  target_name <- targets[v]
  
  # Proteção: Pula a variável caso ela tenha sido deletada na etapa de limpeza de NAs
  if(!target_name %in% colnames(df)) {
    cat(sprintf("\n[ALERTA] %s ausente no df.rda. Pulando...\n", target_name))
    next
  }
  
  cat(sprintf("\nEstimando Projeções Locais (IRF) para %s...\n", target_name))
  
  # A MECÂNICA DO JORDÀ (2005):
  # A matriz Xmat contém toda a informação disponível hoje (tempo t)
  Xmat <- df[1:(nrow(df) - H), ]
  
  # A matriz Ymat contém o futuro. Cada coluna de Ymat é a variável alvo 'empurrada' h meses para frente ($y_{t+h}$)
  Y_target <- df[, target_name]
  Ymat <- matrix(NA, nrow = nrow(df) - H, ncol = H)
  for(h in 1:H) {
    Ymat[, h] <- Y_target[(1 + h):(nrow(df) - H + h)]
  }
  
  # Heurística do Coulombe: 
  # Para poupar poder computacional, ele estima um Lambda base via Regressão Lasso 
  # simples (sem TVP) para usar como âncora nos modelos dinâmicos.
  cvvec = c()
  for(h in 1:H){
    CV = cv.glmnet(x=Xmat, y=Ymat[,h], family='gaussian', alpha=0)
    cvvec = append(cvvec, CV$lambda.min)
  }
  lambda2_base = mean(cvvec)/2
  
  # Loop para cada horizonte preditivo (de 1 a 12 meses à frente)
  for(h in 1:H) {
    cat(sprintf("  -> Horizonte %d\n", h))
    
    # Tratamento de erro (tryCatch): Se a matriz de covariância ficar singular 
    # devido a forte multicolinearidade nas janelas, o R anota 0 (choque nulo) e não quebra o loop.
    tryCatch({
      
      # Estima a Regressão Ridge com Parâmetros Variando no Tempo (TVP)
      out <- TVPRR_cosso(X=Xmat, y=Ymat[,h], type=2, lambdavec=lambdavec, 
                         lambda2=lambda2_base, kfold=5, silent=1)
      
      # O PULO DO GATO: Extraindo o Beta do Choque
      # A matriz de betas gerada tem dimensão [equação, variável explicativa, tempo].
      # Nós queremos extrair APENAS a linha referente à variável 'shock_var' (ex: SPREAD).
      # Somamos +1 ao 'shock_pos' porque a primeira coluna da matriz de betas é sempre o intercepto ($\alpha$).
      IRF_TVP[h, v, ] <- out$grrats$betas_grr[1, shock_pos + 1, ]
      
    }, error = function(e) {
      message(paste("    Matriz singular. Imputando resposta 0 para o horizonte", h))
      IRF_TVP[h, v, ] <- 0
    })
  }
  
  # Salva o arquivo de IRF para cada variável conforme fica pronto
  filename = paste(paths$output, sprintf("Local_Projections_%s.rda", target_name), sep="/")
  save(Ymat, Xmat, IRF_TVP, file=filename)
}

cat("\nProjeções Locais completas! Arquivos salvos na pasta 30_output.\n")
