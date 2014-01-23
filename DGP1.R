DGP <- function(T, a, b, epsilon.params, epsilon.specification="gaussian", initial.y=0) {

    generate.epsilon <- function(spec, params) {
        #cat('generating epsilons of type: ', spec, ' with params: ', params, '\n')
        if (spec=='gaussian')  # epsilon.params is supposed to contain only a value for the std.deviation of epsilon
            return(rnorm(T, 0, params))
        if (spec=='ARCH2') {  # epsilon.params is supposed to contain (a, b, c, epsilon[0], epsilon[1]) the three params for the AR(2) process of the variance process, sigma of the error term of the variance process, as well as the initial values for the variance process. I have problems with negative values but bollerslev engle et al (82) write the constant "almost surely needs to be > 0 for positive variance"
            ht <- numeric(T)
            nu <- runif(T, 0, 1)  #TODO: uniform is probably not what she had in mind when she wrote vt~iid(0,1)
            epsilon <- c(params[4], params[5], rep(0, T-2))
            for (i in 3:T) {
                ht[i] <- params[1] + params[2]*epsilon[i-1]^2 + params[3]*epsilon[i-2]^2  # AR(2) specification of epsilon^2
                epsilon[i] <- sqrt(ht[i]) * nu[i]
            }
            return(epsilon)  # dont forget to take sqrt
        }
        if (spec=='GARCH11') {  # epsilon.params is supposed to contain (a, b, c, epsilon[0], h[0])
            ht <- numeric(T)
            nu <- runif(T, 0, 1)  #TODO: uniform is probably not what she had in mind when she wrote vt~iid(0,1)
            epsilon <- c(params[4], rep(0, T-1))
            ht[1] <- params[5]
            for (i in 3:T) {
                ht[i] <- params[1] + params[2]*epsilon[i-1]^2 + params[3]*ht[i-1]  # AR(2) specification of epsilon^2
                epsilon[i] <- sqrt(ht[i]) * nu[i]
            }
            return(epsilon)  # dont forget to take sqrt
        }
        if (spec=='GJR111') {  # epsilon.params is supposed to contain (a, b, c, d, epsilon[0], h[0])
            ht <- numeric(T)
            nu <- runif(T, 0, 1)  #TODO: uniform is probably not what she had in mind when she wrote vt~iid(0,1)
            epsilon <- c(params[5], rep(0, T-1))
            ht[1] <- params[6]
            I <- epsilon < 0
            for (i in 3:T) {
                ht[i] <- params[1] + params[2]*epsilon[i-1]^2 + params[3]*ht[i-1] + params[4]*epsilon[i-1]^2*I[i-1]
                epsilon[i] <- sqrt(ht[i]) * nu[i]
                I[i] <- epsilon[i] < 0
            }
            return(epsilon)  # dont forget to take sqrt
        }
    }

    epsilon <- generate.epsilon(epsilon.specification, epsilon.params)  # epsilon is generated according to model specification
    y <- epsilon
    y[0] <- initial.y
    for (i in 2:T)
        y[i] <- a + b*y[i-1] + y[i]  # AR(1) process for y
    return(data.frame(y=y, epsilon=epsilon))
}
