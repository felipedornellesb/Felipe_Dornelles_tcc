source('00_prog/00_setup.R')
load('10_data/data.rda')
dates_col <- data$date
data$date <- NULL
rownames(data) <- as.character(dates_col)
variable <- 'CPIAUCSL'
nwindows <- 180
maxh <- 12
options(error = traceback)
res <- tryCatch({
  rolling_window(runlasso, data, nwindows, maxh, variable)
}, error=function(e) {
  print(e)
  traceback()
})
