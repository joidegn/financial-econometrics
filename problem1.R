source('DGP1.R')
require('zoo')  # to make lagging variables easier

# #######################   1)   ####################### #
make.graphs <- function() {
    for (T in c(200, 500, 1000)) {
        for (alpha in seq(0, 1, 0.1)) {
            cat('alpha: ', alpha, '\n')
            data.A <- DGP(T, 0, alpha, 1)
            png(paste0('graphs/problem1_ModelA_T', T, '_alpha', alpha), width=640, height=480)
            plot(data.A[, 1], type='l', main=paste0('AR(1) T=', T, ' a=0, b=', alpha), xlab='t', ylab='y')
            dev.off()
        }
    }  # because why the hell not?

#dev.new()
#plot(DGP(500, 0.5, 0.5, 1)[,1], type='l')  # positive drift (i.e. a>0)

    for (T in c(200, 500, 1000)) {
        for (alpha in seq(0, 1, 0.1)) {
            cat('alpha: ', alpha, '\n')
            data.A <- DGP(T, 0, alpha, c(0.1, 0.5, 0.4, 1, 1), 'ARCH2')  # ARCH 2 data with constant 0.1, alpha1=0.5, alpha2=0.4, initial values 1 and 1
            png(paste0('graphs/problem1_ModelB_T', T, '_alpha', alpha), width=640, height=480)
            plot(data.A[, 1], type='l', main=paste0('ARCH(2) T=', T, ' a=0, b=', alpha), xlab='t', ylab='y')
            points(data.A[, 2], type='l', col="green")
            dev.off()
        }
    }  # because why the hell not?

    for (T in c(200, 500, 1000)) {
        for (alpha in seq(0, 1, 0.1)) {
            cat('alpha: ', alpha, '\n')
            data.A <- DGP(T, 0, alpha, c(0.1, 0.5, 0.4, 1, 1), 'GARCH11')  # GARCH 1,1 data with constant 0.1, alpha1=0.5, alpha2=0.4, initial values 1 and 1
            png(paste0('graphs/problem1_ModelB_GARCH_T', T, '_alpha', alpha), width=640, height=480)
            plot(data.A[, 1], type='l', main=paste0('GARCH(1,1) T=', T, ' a=0, b=', alpha), xlab='t', ylab='y')
            points(data.A[, 2], type='l', col="green")
            dev.off()
        }
    }  # because why the hell not?

    for (T in c(200, 500, 1000)) {
        for (alpha in seq(0, 1, 0.1)) {
            cat('alpha: ', alpha, '\n')
            data.A <- DGP(T, 0, alpha, c(0.1, 0.5, 0.4, 0.5, 1, 1), 'GJR111')  # GJR 1,1,1 data with constant 0.1, alpha1=0.5, alpha2=0.4, alpha3=0.5, initial values 1 and 1
            png(paste0('graphs/problem1_ModelB_GJR_T', T, '_alpha', alpha), width=640, height=480)
            plot(data.A[, 1], type='l', main=paste0('GJR(1,1,1) T=', T, ' a=0, b=', alpha), xlab='t', ylab='y')
            points(data.A[, 2], type='l', col="green")
            dev.off()
        }
    }  # because why the hell not?
}
make.graphs()

# #######################   2)   ####################### #

# Engle test (ARCH-test)
# Model A:
run.tests <- function() {
    print(system.time({
        N <- 1000  # we repeat the tests a 1000 times per data series
        results <- matrix(rep(NA, N*3), ncol=3)
        colnames(results) <- c(200, 500, 1000)
        args.list <- list(  # so I can call DGP in a loop
            list(a=0, b=0.5, epsilon.params=1, epsilon.specification="gaussian", initial.y=0),
            list(a=0, b=0.5, epsilon.params=c(0.1, 0.5, 0.4, 1, 1), epsilon.specification="ARCH2", initial.y=0),
            list(a=0, b=0.5, epsilon.params=c(0.1, 0.5, 0.5, 1, 1), epsilon.specification="GARCH11", initial.y=0),
            list(a=0, b=0.5, epsilon.params=c(0.1, 0.5, 0.5, 0.5, 1, 1), epsilon.specification="GJR111", initial.y=0)
        )
        cat('Engle-tests with', N, 'repetitions for AR(1), ARCH(2), GARCH(1,1), GJR(1,1,1) models with T=(200, 500, 1000). This will take a while! \n')
        for (model.num in 1:4) {
            for (T in c(200, 500, 1000)) {
                cat('Starting tests with T=', T, '\n')
                critical.value <- qchisq(0.95, 2)
                for (n in 1:N) {
                    #data <- zoo(DGP(T, 0, 0.5, 1))
                    data <- zoo(do.call(DGP, c(list(T=T), args.list[[model.num]])))
                    # conduct Engle test with m=2 lagged squared residuals 
                    data$lagged.y <- lag(data$y, 1, na.pad=T)  # note that lag has been overwritten by the zoo package
                    reg.y <- lm(y~lagged.y, data=data)
                    # TODO: is this correct? not too many rejections?
                    residual.reg.data <- data.frame(resid.sq=reg.y$residuals[1:(T-2)]^2, resid.sq.lagged=lag(zoo(reg.y$residuals))^2)
                    reg.residuals <- lm(resid.sq~resid.sq.lagged, data=residual.reg.data)  # regress epsilon^2 on lags q lags of epsilon^2 and a constant
                    #f.stat <- summary(reg)$fstatistic[1]  # I read somewhere to use that 
                    TR.sq <- (nrow(data) - 2) * summary(reg.residuals)$r.squared
                    results[n, as.character(T)] <- TR.sq > critical.value
                }
            }
            rejections <- apply(results, 2, sum) / 1000
            cat('rejections: ', rejections, '\n')
        }
    }))
}
#run.tests()

# #######################   3)   ####################### #
T <- 1000  # discard first 500 values
# for now use an AR(1) without drift
sample <- DGP(T=1000, a=0, b=0.8, epsilon.params=c(0.3, 0.4, 0.4, 0.3, 0.1, 0.1), epsilon.specification="GJR111", initial.y=0)[500:1000, ]
plot(sample$y, type='l')
points(sample$epsilon, type='l', col="green", lty=2, xlab="time")

#loglik.garch.gjr <- function(params, data) {  # params consists of 9 parameters: 2 AR params and 6 GJR params
#    ht <- params[8]
#    lik <- 0
#    for (i in 3:(nrow(data)+1)) {
#        resid <- data$y[i-1] - (params[1] + params[2]*data$y[i-2])  # note that this is resid[i-1] so to say, i.e. residual of the previous period
#        ht <- params[3] + params[4]*resid^2 + params[5]*ht + params[6]*resid^2*(resid<0)
#        lik <- lik + (-0.5*log(2*pi*ht)-0.5*(resid/ht))
#    }
#    return(lik)
#}
#loglik.garch <- function(params, data) {  # params consists of 7 parameters: 2 AR params and 5 GJR params
#    ht <- params[7]
#    lik <- 0
#    for (i in 3:(nrow(data)+1)) {
#        resid <- data$y[i-1] - (params[1] + params[2]*data$y[i-2])  # note that this is resid[i-1] so to say, i.e. residual of the previous period
#        ht <- params[3] + params[4]*resid^2 + params[5]*ht
#        lik <- lik + (-0.5*log(2*pi*ht)-0.5*(resid/ht))
#    }
#    return(lik)
#}
#optim(c(0, 0.6, 0.4, 0.3, 0.3, 0.4, 0.2, 0.2), function(params) loglik.garch.gjr(params, sample))

library(rugarch)
ugarchfit.gjr.model <- ugarchfit(ugarchspec(mean.model=list(armaOrder=c(1,0)), variance.model=list(model="gjrGARCH", garchOrder=c(1,1,1))), data=sample)
# #######################   4)   ####################### #
ugarchfit.garch.model <- ugarchfit(ugarchspec(mean.model=list(armaOrder=c(1,0)), variance.model=list(model="fGARCH", submodel="GARCH", garchOrder=c(1,1))), data=sample)
