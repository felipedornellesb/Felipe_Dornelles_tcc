source('00_prog/00_setup.R')
load('10_data/data.rda')
dates_col <- data$date
data$date <- NULL
rownames(data) <- as.character(dates_col)
variable <- 'CPIAUCSL'
nwindows <- 1
maxh <- 1
options(error = traceback)
res <- try(rolling_window(runlasso, data, nwindows, maxh, variable))
if (inherits(res, 'try-error')) print(geterrmessage())
