source('DGP1.R')
require('zoo')  # to mag lagging variables easier

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
            data.A <- DGP(T, 0, alpha, c(0.1, 0.5, 0.4, 1, 1), 'ARCH2')  # ARCH 2 data with constant 0, alpha1=0.5, alpha2=0.4, initial values 1 and 1
            png(paste0('graphs/problem1_ModelB_T', T, '_alpha', alpha), width=640, height=480)
            plot(data.A[, 1], type='l', main=paste0('ARCH(2) T=', T, ' a=0, b=', alpha), xlab='t', ylab='y')
            points(data.A[, 2], type='l', col="green")
            dev.off()
        }
    }  # because why the hell not?
}
# #######################   2)   ####################### #

# Engle test (ARCH-test)
# Model A:
system.time({
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
                data$epsilon.lag1 <- lag(data$epsilon, -1, na.pad=T)  # note that lag has been overwritten by the zoo package
                data$epsilon.lag2 <- lag(data$epsilon, -2, na.pad=T)  # note that lag has been overwritten by the zoo package
                reg <- lm(epsilon^2~epsilon.lag1^2+epsilon.lag2^2, data=data)  # regress epsilon^2 on lags q lags of epsilon^2 and a constant
                #f.stat <- summary(reg)$fstatistic[1]  # I read somewhere to use that 
                TR.sq <- (nrow(data) - 2) * summary(reg)$r.squared
                results[n, as.character(T)] <- TR.sq > critical.value
            }
        }
        rejections <- apply(results, 2, sum) / 1000
        cat('rejections: ', rejections, '\n')
    }
})