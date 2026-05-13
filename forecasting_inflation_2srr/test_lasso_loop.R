source('00_prog/00_setup.R')
load('10_data/data.rda')
dates_col <- data$date
data$date <- NULL
rownames(data) <- as.character(dates_col)
variable <- 'CPIAUCSL'
horizon <- 1
window_size <- nrow(data) - 180
cat('Running for windows...\n')
for (i in 1:180) {
  ind <- i:(window_size + i - 1)
  res <- try(runlasso(ind, data, variable, horizon), silent=TRUE)
  if (inherits(res, 'try-error')) {
    cat('Failed at window:', i, '\n')
    cat('Error:', attr(res, 'condition')$message, '\n')
    break
  }
}
