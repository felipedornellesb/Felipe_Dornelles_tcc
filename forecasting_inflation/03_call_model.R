library(here)
setwd(here("forecasting_inflation"))

install.packages("pak")

### must add package for specific models ###
 library(devtools)
 pak::pak("gabrielrvsc/HDeconometrics")
library(HDeconometrics)
library(glmnet)
library(randomForest)
library(tidyverse)

source("functions/rolling_window.R")
source("functions/functions.R")

#####
## The file with the forecasts will be saved with model_name
model_name <- "T.Factor"
## The function called to run models is model_function, which is a function from functions.R
model_function <- runtfact
#####


load("data/data.rda")
dates <- data$date
data <- data %>%
  select(-date) %>%
  as.matrix()
rownames(data) <- as.character(dates)

####### run rolling window ##########
nwindows <- 312
model_list <- list()
for (i in 1:12) {
  model <- rolling_window(model_function, data, nwindows + i - 1, i, "CPIAUCSL")
  model_list[[i]] <- model
  cat(i, "\n")
}

forecasts <- Reduce(cbind, lapply(model_list, function(x) head(x$forecast, nwindows)))

forecasts <- accumulate_model(forecasts)

save(forecasts, file = paste("forecasts/", model_name, ".rda", sep = ""))


plot(tail(data[, "CPIAUCSL"], nwindows), type = "l")
lines(forecasts[, 1], col = 3)
