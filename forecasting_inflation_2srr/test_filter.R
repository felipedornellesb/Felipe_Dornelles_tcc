y <- 1:10
h <- 3
y_h <- as.numeric(stats::filter(y, rep(1, h), sides = 1))
y_h[is.na(y_h)] <- y_h[which.min(is.na(y_h))]
print(y_h)
