# =============================================================================================================
# INITIALIZATION & DIRECTORY SETUP
# =============================================================================================================

rm(list = ls()) # Clear environment

# INSERT YOUR PATH HERE
wd = 'C:/Users/felip/OneDrive/Documentos/tcc/Felipe_Dornelles_tcc'
setwd(wd)

# Define paths relative to the working directory
paths <- list(program = "00_program",
              data = "10_data",
              tools = "20_tools",
              functions = "20_tools/functions",
              output = "30_output",
              results = "40_results")

# =============================================================================================================
# PACKAGE MANAGEMENT
# =============================================================================================================

myPKGs <- c('tidyverse', 'forecast', 'MCS', 'xtable', 'ggplot2', 'ggsci', 'ggpubr')

InstalledPKGs <- names(installed.packages()[,'Package'])
InstallThesePKGs <- myPKGs[!myPKGs %in% InstalledPKGs]
if (length(InstallThesePKGs) > 0) install.packages(InstallThesePKGs, repos = "http://cran.us.r-project.org")

invisible(lapply(myPKGs, library, character.only = TRUE))

# =============================================================================================================
# LOAD DATA & SETUP VARIABLES
# =============================================================================================================

# Load actual values (The "Ground Truth")
actual_values <- readRDS(paste(paths$output, "actual_values.rda", sep = "/"))
test_matrix <- as.matrix(actual_values)

# Variables and Horizons
UNRATE <- "SEADE12_TDTGSP12"
IPCA <- "PRECOS12_IPCA12"
SPREAD <- "JPM366_EMBI366"

variable <-  c(IPCA, UNRATE, SPREAD)
horizon_list <- c(1, 3, 6, 12)

# Create an empty data frame to store RMSPE results
rmspe_data <- data.frame(variable = character(),
                         horizon = numeric(),
                         model = character(),
                         RMSPE = numeric(),
                         stringsAsFactors = FALSE)

# =============================================================================================================
# LOAD FORECASTS
# =============================================================================================================

# Carregue os modelos que você já rodou no 02_forecast.R
coulombe_2srr <- readRDS(paste(paths$output, "coulombe_2srr.rda", sep = "/"))

# Combine all loaded lists into a single master list
model_list <- c(coulombe_2srr) 
models_list_names <- c("coulombe_2srr")

# =============================================================================================================
# CALCULATE ERRORS (RMSPE)
# =============================================================================================================

matrix_list <- list()
squared_error_list <- list()

for (v in variable) {
  for (h in horizon_list) {
    for (model in models_list_names) {
      
      # Filter the list for the specific variable, horizon, and model
      model_list_filtered <- model_list[sapply(model_list, function(x) x$horizon == h & x$variable == v & x$model == model)]
      
      # Check if the model actually ran for this combination
      if (length(model_list_filtered) > 0) {
        
        # Extract predictions
        model_matrix <- sapply(model_list_filtered, function(x) x$pred)
        
        # PEGA APENAS A COLUNA DA VARIÁVEL ALVO NA MATRIZ REAL
        actual_target <- test_matrix[, v]
        
        # Calculate Prediction Error (Actual - Predicted)
        pred_error <- (actual_target - as.vector(model_matrix))
        matrix_list[[model]] <- pred_error[is.finite(pred_error)]^2
        
        # Calculate RMSPE (Root Mean Squared Predictive Error)
        rmspe <- sqrt(mean(pred_error[is.finite(pred_error)]^2, na.rm = TRUE))
        
        # Bind to results table
        rmspe_data <- rbind(rmspe_data, data.frame(variable = v, horizon = h, model = model, RMSPE = rmspe))
      }
    }
    
    # Store the squared error values for potential MCS evaluation later
    if (length(matrix_list) > 0) {
      result_matrix <- do.call(cbind, matrix_list)
      squared_error_list[[paste(v, h, sep="_")]] <- result_matrix
    }
  }
}

# =============================================================================================================
# EXPORT RESULTS (TABLES & LATEX)
# =============================================================================================================

# Format table
rmspe_data$RMSPE <- round(rmspe_data$RMSPE, 4)

# Export LaTeX Tables for each variable
export_latex_table <- function(data, var_name, path) {
  var_data <- data %>%
    filter(variable == var_name) %>%
    select(model, horizon, RMSPE) %>%
    pivot_wider(names_from = horizon, names_prefix = "h=", values_from = RMSPE)
  
  if (nrow(var_data) > 0) {
    latex_code <- print(xtable(var_data, caption = paste("RMSPE for", var_name)), type = 'latex', comment = FALSE)
    writeLines(latex_code, paste(path, paste0("table_RMSPE_", var_name, ".tex"), sep = "/"))
  }
}

export_latex_table(rmspe_data, SPREAD, paths$results)
export_latex_table(rmspe_data, IPCA, paths$results)
export_latex_table(rmspe_data, UNRATE, paths$results)

cat("\nEvaluation complete! LaTeX tables saved to 40_results folder.\n")
print(rmspe_data)

# Criando a imagem do gráfico
plot_rmspe <- ggplot(rmspe_data, aes(x = as.factor(horizon), y = RMSPE, group = variable, color = variable)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 4) +
  scale_color_npg() +
  labs(title = "Erro Quadrático Médio de Previsão (RMSPE)",
       subtitle = "Avaliação do modelo Coulombe 2SRR por horizonte preditivo",
       x = "Horizonte de Previsão (Meses)",
       y = "RMSPE",
       color = "Variável Alvo") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        legend.position = "bottom",
        legend.title = element_text(face = "bold"))

# Visualizar o gráfico na sua tela
print(plot_rmspe)

# Salvar a imagem em alta resolução (PNG) na sua pasta 40_results
ggsave(filename = paste(paths$results, "grafico_RMSPE_Coulombe.png", sep = "/"), 
       plot = plot_rmspe, 
       width = 8, 
       height = 5, 
       dpi = 300, 
       bg = "white")

cat("\nImagem gerada e salva com sucesso na pasta 40_results!\n")

