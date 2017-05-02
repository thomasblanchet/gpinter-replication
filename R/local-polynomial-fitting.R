# ---------------------------------------------------------------------------- #
# “Generalized Pareto Curves: Theory and Applications”, 2017
# Thomas Blanchet, Juliette Fournier, Thomas Piketty
# ---------------------------------------------------------------------------- #
# Local polynomial fitting on a function and its derivative.
# ---------------------------------------------------------------------------- #

# Local polynomial fitting on a function and its derivative
lpoly <- function(x0, x, y, dydx, alpha) {
    n <- length(x)
    t <- x - x0
    h <- sort(abs(t))[floor(alpha*n)]
    w <- dnorm(t, sd=h)/h

    Y <- c(y, dydx)

    X0 <- c(rep(1, n), rep(0, n))
    X1 <- c(t, rep(1, n))
    X2 <- c(t^2/2, t)
    X3 <- c(t^3/6, t^2/2)

    fit <- lm(Y ~ 0 + X0 + X1 + X2 + X3, weights=c(w, w))
    return(coef(fit)["X3"])
}
