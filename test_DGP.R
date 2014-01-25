source('DGP1.R')

graphs <- function() {
    # plot should show normally distributed error terms
    plot(DGP(T=1000, a=0, b=0.5, epsilon.params=c(1), epsilon.specification="gaussian", initial.y=0)$epsilon, type='l')
    devAskNewPage(T)
    # plot should show some normal ARCH(2,2) effects
    plot(DGP(T=1000, a=0, b=0.5, epsilon.params=c(0.3, 0.7, 0.5, 0.1, 0.1), epsilon.specification="ARCH2", initial.y=0)$epsilon, type='l')
    # plot should show some normal GARCH(1,1) effects
    plot(DGP(T=1000, a=0, b=0.5, epsilon.params=c(0.3, 0.4, 0.4, 0.1, 0.1), epsilon.specification="GARCH11", initial.y=0)$epsilon, type='l')
    # plot should show some normal GJR(1,1,1) effects
    plot(DGP(T=1000, a=0, b=0.5, epsilon.params=c(0.3, 0.4, 0.4, 0.5, 0.1, 0.1), epsilon.specification="GJR111", initial.y=0)$epsilon, type='l')
    devAskNewPage(F)
}
