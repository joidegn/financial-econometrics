DGP <- function(T, a, b, epsilon.params, epsilon.specification="gaussian", initial.y=0) {

    generate.epsilon <- function(spec, params) {
        #cat('generating epsilons of type: ', spec, ' with params: ', params, '\n')
        if (spec=='gaussian')  # epsilon.params is supposed to contain only a value for the std.deviation of epsilon
            return(data.frame(epsilon=rnorm(T, 0, params), ht=rep(1, T)))
        if (spec=='ARCH2') {  # epsilon.params is supposed to contain (a, b, c, epsilon[0], epsilon[1]) the three params for the AR(2) process of the variance process, sigma of the error term of the variance process, as well as the initial values for the variance process. I have problems with negative values but bollerslev engle et al (82) write the constant "almost surely needs to be > 0 for positive variance"
            nu <- rnorm(T, 0, 1)  # TODO: should I take a gaussian here?
            epsilon <- c(params[4], params[5], rep(0, T-2))
            ht <- c(sqrt(params[4]), sqrt(params[5]), rep(0, T-2))  # TODO: initialize ht as sqrt of initial epsilon values?
            for (i in 3:T) {
                ht[i] <- params[1] + params[2]*epsilon[i-1]^2 + params[3]*epsilon[i-2]^2  # AR(2) specification of epsilon^2
                epsilon[i] <- sqrt(ht[i]) * nu[i]
            }
            return(data.frame(epsilon=epsilon, ht=ht))
        }
        if (spec=='GARCH11') {  # epsilon.params is supposed to contain (a, b, c, epsilon[0], h[0])
            ht <- numeric(T)
            nu <- rnorm(T, 0, 1)  #TODO: is this right?
            epsilon <- c(params[4], rep(0, T-1))
            ht[1] <- params[5]
            for (i in 3:T) {
                ht[i] <- params[1] + params[2]*epsilon[i-1]^2 + params[3]*ht[i-1]  # AR(2) specification of epsilon^2
                epsilon[i] <- sqrt(ht[i]) * nu[i]
            }
            return(data.frame(epsilon=epsilon), ht=ht)
        }
        if (spec=='GJR111') {  # epsilon.params is supposed to contain (a, b, c, d, epsilon[0], h[0])
            ht <- numeric(T)
            nu <- rnorm(T, 0, 1)
            epsilon <- c(params[5], rep(0, T-1))
            ht[1] <- params[6]
            I <- epsilon < 0
            for (i in 3:T) {
                ht[i] <- params[1] + params[2]*epsilon[i-1]^2 + params[3]*ht[i-1] + params[4]*epsilon[i-1]^2*I[i-1]
                epsilon[i] <- sqrt(ht[i]) * nu[i]
                I[i] <- epsilon[i] < 0
            }
            return(data.frame(epsilon=epsilon, ht=ht))
        }
    }

    res <- generate.epsilon(epsilon.specification, epsilon.params)  # epsilon is generated according to model specification
    epsilon <- res$epsilon
    ht <- res$ht
    y <- epsilon
    y[0] <- initial.y
    for (i in 2:T)
        y[i] <- a + b*y[i-1] + y[i]  # AR(1) process for y
    return(data.frame(y=y, epsilon=epsilon, ht=ht))  # return also conditional variance ht
}
