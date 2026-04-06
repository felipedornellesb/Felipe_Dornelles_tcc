# =============================================================================================================
# INITIALIZATION & DIRECTORY SETUP
# =============================================================================================================

rm(list = ls()) # Clear environment

# INSERT YOUR PATH HERE
wd = 'C:/Users/felip/OneDrive/Documentos/tcc/Felipe_Dornelles_tcc/'
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

myPKGs <- c('dplyr', 'ipeadatar', 'readxl', 'lubridate', 'urca', 'tidyr', 'tseries')

InstalledPKGs <- names(installed.packages()[,'Package'])
InstallThesePKGs <- myPKGs[!myPKGs %in% InstalledPKGs]
if (length(InstallThesePKGs) > 0) install.packages(InstallThesePKGs, repos = "http://cran.us.r-project.org")

# Load libraries
invisible(lapply(myPKGs, library, character.only = TRUE))

# Load data functions
source(paste(paths$functions, '00_Nathalia_functions.R', sep='/'))

# =============================================================================================================
# 1. GET DATA
# =============================================================================================================

# Read dataset definitions
dataset <- read_excel(paste(paths$data, "dataset.xlsx", sep='/'),
                      col_types = c("text", "text", "text", "date", "date", "numeric", "text"))

# Target Variables
UNRATE <- "SEADE12_TDTGSP12"
IPCA <- "PRECOS12_IPCA12"
SPREAD <- "JPM366_EMBI366"

# Acquiring the metadata from each CODE
metadados <- metadata(dataset$codigo)
metadados_spread <- metadata("JPM366_EMBI366")

# Acquiring the values from each CODE
data <- ipeadata(metadados$code)
data_spread <- ipeadata(metadados_spread$code)

# Converting data to wide format
df <- data %>%
  pivot_wider(names_from = "code") %>%
  select(-c(uname, tcode))
df_sorted <- df[order(as.Date(df$date, format = "%Y/%m/%d")), ]

df_spread <- data_spread %>%
  pivot_wider(names_from = "code") %>%
  select(-c(uname, tcode))

# Transform SPREAD from daily data to monthly data using mean
df_spread <- df_spread %>%
  mutate(year_month = format(date, "%Y-%m")) %>%  
  group_by(year_month) %>%                        
  summarise(JPM366_EMBI366 = mean(JPM366_EMBI366)) %>% 
  mutate(date = as.Date(paste(year_month, "-01", sep = ""))) %>%  
  select(date, JPM366_EMBI366) 

# Merge results in the same df
df <- left_join(df_sorted, df_spread, by = "date")

# Filtering by dates
df_filtered <- subset(df, date >= as.Date("1996-01-01") & date < as.Date("2019-06-01"))
df_filtered <- df_filtered  %>% select_if(~ !any(is.na(.)))


# =============================================================================================================
# 2. STATIONARITY TESTS & TRANSFORMATIONS
# =============================================================================================================

non_stationary <-  df

# Identifying the columns types
tipo <- c(4)
for (i in 2:ncol(non_stationary)) {
  if (sum(non_stationary[, i] < 0, na.rm = T) > 0) { 
    if (sum(non_stationary[, i] == 0, na.rm = T) > 0) { 
      tipo <- append(tipo, 3) 
    } else {
      if (sum(non_stationary[, i] < 0, na.rm = T) == dim(non_stationary[, i])[1]) {
        tipo <- append(tipo,5) 
      } else {
        tipo <- append(tipo, 1) }
    }
  } else { 
    if (sum(non_stationary[, i] == 0, na.rm = T) > 0) { 
      tipo <- append(tipo, 2)
    } else {
      tipo <- append(tipo, 0) 
    }
  }
}

# ADF Test
test <- c(4)
for (i in 2:(ncol(non_stationary))) {
  cat(sprintf("Testing variable %d of %d\r", i, ncol(non_stationary)))
  X <- na.exclude(as.matrix(non_stationary[, i]))
  k <- 0
  j <- 0
  status <- "non-stationary"
  
  while (status == "non-stationary") {
    adf_test <- summary(ur.df(X, "none", lags = 12))
    if (k == 2) {
      dale <- colnames(non_stationary)[i]
    }
    
    if (adf_test@teststat[1] <= adf_test@cval[1, 2] && kpss.test(X, null = "T")$p.value >= 0.05) {
      status <- "stationary"
      if (j == 1) {
        k <- 0.5
      }
    } else {
      if (j == 0) {
        X <- preparacao(X, i)
        j <- 1
      } else {
        k <- k + 1
        if (tipo[i] == 0 | tipo[i] == 2 | tipo[i] == 3 | tipo[i]==1) {
          X <- diff(X)
        }
        if (tipo[i] == 5) {
          X <- cresc_discreto(X)
        }
        j <- j + 1
      }
    }
  }
  test <- append(test, k)
}
cat("\nStationarity testing complete.\n")

# Assessing transformation value
transformation <- c(-1)
for (i in 2:ncol(non_stationary)) {
  if (tipo[i] == 0 && test[i] == 0) {
    transformation[i] <- 4
  } else if (tipo[i] == 0 && test[i] == 1) {
    transformation[i] <- 5
  } else if (tipo[i] == 0 && test[i] == 0.5) {
    transformation[i] <- 2
  } else if (tipo[i] == 0 && test[i] == 2) {
    transformation[i] <- 6
  } else if (tipo[i] == 1 && test[i] == 0) {
    transformation[i] <- 1
  } else if (tipo[i] == 1 && test[i] == 1) {
    transformation[i] <- 2
  } else if (tipo[i] == 1 && test[i] == 2) {
    transformation[i] <- 3
  } else if (tipo[i] == 2 && test[i] == 0) {
    transformation[i] <- 1
  } else if (tipo[i] == 3 && test[i] == 0) {
    transformation[i] <- 1
  } else if (tipo[i] == 3 & test[i] == 1) {
    transformation[i] <- 2
  } else if (tipo[i] == 2 & test[i] == 1) {
    transformation[i] <- 2
  } else if (tipo[i] == 5 && test[i] == 1) {
    transformation[i] <- 7}
}

month <- t(as.matrix(month(df_filtered$date)))
year <- t(as.matrix(year(df_filtered$date)))
reference <- available_subjects()
metadados <- metadados %>% inner_join(reference, by = "code")
nome <- names(df_filtered[2:ncol(df_filtered)])
metadados <- metadados[metadados$code %in% nome,]
dataset <- list(
  "data" = as.matrix(df_filtered[2:ncol(df_filtered)]),
  "transformation" = t(as.matrix(transformation[2:ncol(df_filtered)])),
  month = month, year = year, "metadados" = metadados
)
dataset$names <- names(as.data.frame(dataset$data))
names(dataset$data) <- NULL

# =============================================================================================================
# 3. APPLY TRANSFORMATIONS AND EXPORT
# =============================================================================================================

transformed_dataset <- data.frame(matrix(NA, nrow = nrow(dataset$data), ncol = ncol(dataset$data)))
for (i in (1:ncol(dataset$data))){
  transformed_dataset[,i] <- transform_singlestep(dataset$data[,i],dataset$transform[i])
}
names(transformed_dataset) <- dataset$names

# Remove NA values generated by lags/differences
df <- transformed_dataset[-c(1:2),] %>% select_if(~ !any(is.na(.)))

# Save df to data folder
save(df, file = paste(paths$data,"df.rda",sep = "/"))
cat("Data preparation complete. Output saved to 10_data/df.rda\n")

# Load para see the results
load("C:/Users/felip/OneDrive/Documentos/tcc/Felipe_Dornelles_tcc/10_data/df.rda")

View(df)