library(here)
setwd(here("forecasting_inflation"))

#### gets out of sample y and computes random walk forecasts ###
library(roll)
library(tidyverse)

load("data/data.rda")
dates = data$date
data = data%>%select(-date)%>%as.matrix()
rownames(data) = as.character(dates)

nwindows = 312

# CPIAUCSL is 100 * log_diff(CPI), i.e. monthly inflation in pp. Log diffs are
# additive, so h-period accumulated inflation is the SUM. The previous
# roll_prod(1+y,h)-1 treated pp as a multiplicative gross return, inflating the
# target exponentially (RW RMSE at h=12 was ~72 instead of ~2).
y = data[,"CPIAUCSL"]
y = cbind(y,roll_sum(y,3),roll_sum(y,6),roll_sum(y,12))
yout = tail(y,nwindows)

rw = matrix(NA,nwindows,12)
for(i in 1:12){
  aux=data[(nrow(data)-nwindows-i+1):(nrow(data)-i),"CPIAUCSL"]
  rw[,i]=aux;
}

rw3 = tail(embed(y[,2],4)[,4],nwindows)
rw6 = tail(embed(y[,3],7)[,7],nwindows)
rw12 = tail(embed(y[,4],13)[,13],nwindows)
rw = cbind(rw,rw3,rw6,rw12)
colnames(rw) = c(paste("t+",1:12,sep = ""),"acc3","acc6","acc12")

save(yout,file = "forecasts/yout.rda")
save(rw,file = "forecasts/rw.rda")



