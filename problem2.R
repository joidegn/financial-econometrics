# check if mgarch is installed and if not install from github (it is not available on CRAN)
if (! 'mgarch' %in% installed.packages()[, 'Package']) {
    require('devtools') || (install.packages('devtools') && require('devtools'))
    install_github('vst/mgarch', subdir='package/mgarch')
}

#library(rmgarch)
params <- c(1,1,1,0.5,0,0,0.5,0.5,0,0,0.5)
series <- mvBEKK.sim(2, 1000, order=c(1,1), params=params)

par(mfrow=c(2,1))
plot(series$eps[[1]], type='l', ylab='yt 1')
plot(series$eps[[2]], type='l', ylab='yt 2')
savePlot('graphs/problem2_simulated_series.png', type='png')

estimated.model <- mvBEKK.est(as.data.frame(series$eps), order=c(1,1), method='BFGS')
cat('estimated Parameters: ', estimated.model$estimation$par, '\n')
