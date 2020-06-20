# ---------------------------------------------------------------------------- #
# “Generalized Pareto Curves: Theory and Applications”, 2017
# Thomas Blanchet, Juliette Fournier, Thomas Piketty
# ---------------------------------------------------------------------------- #
# Functions for the different interpolation methods that exist in the
# litterature.
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# Method 1: Constant Pareto coefficient above a threshold
# ---------------------------------------------------------------------------- #

method1 <- function(p, pk, qk, mk, average) {
    k <- cut(p, breaks=c(pk, 1), include.lowest=TRUE, labels=FALSE)

    # Pareto coefficients
    b <- mk[k]/((1 - pk[k])*qk[k])
    a <- b/(b - 1)
    # Thresholds
    q <- qk[k]*((1 - p)/(1 - pk[k]))^(-1/a)
    # Averages above
    e <- b*q
    # Top shares
    s <- e*(1 - p)/average

    return(list(threshold=q, top_share=s))
}

# ---------------------------------------------------------------------------- #
# Method 2: Piecewise Pareto with threshold information only
# ---------------------------------------------------------------------------- #

method2 <- function(p, pk, qk, average) {
    # Number of thresholds
    n <- length(pk)

    # Pareto coefficient in each bracket
    ak <- -log((1 - pk[2:n])/(1 - pk[1:(n - 1)]))/log(qk[2:n]/qk[1:(n -1)])
    # Truncated average
    uk <- ak/(ak - 1) * qk[1:(n - 1)]^ak * (qk[1:(n - 1)]^(1 - ak) - qk[2:n]^(1 - ak)) * (1 - pk[1:(n-1)])
    # Beyond the last bracket
    uk <- c(uk, ak[n - 1]/(ak[n - 1] - 1)*(1 - pk[n])*qk[n])
    ak <- c(ak, ak[n - 1])
    # Lorenz curve
    mk <- rev(cumsum(rev(uk)))

    k <- cut(p, breaks=c(pk, 1), include.lowest=TRUE, labels=FALSE)

    q <- qk[k]*((1 - p)/(1 - pk[k]))^(-1/ak[k])
    u <- ak[k]/(ak[k] - 1) * qk[k]^ak[k] * (qk[k]^(1 - ak[k]) - q^(1 - ak[k])) * (1 - pk[k])

    return(list(threshold=q, top_share=(mk[k] - u)/average))
}

# ---------------------------------------------------------------------------- #
# Method 3: Mean-split histogram
# ---------------------------------------------------------------------------- #

method3 <- function(p, pk, qk, mk, average) {
    n <- length(pk)

    p0 <- pk[1:(n - 1)]
    p1 <- pk[2:n]
    q0 <- qk[1:(n - 1)]
    q1 <- qk[2:n]
    m0 <- mk[1:(n - 1)]
    m1 <- mk[2:n]

    uk <- -diff(mk)/diff(pk)

    f0 <- (p1 - p0)/(q1 - q0)*(q1 - uk)/(uk - q0)
    f1 <- (p1 - p0)/(q1 - q0)*(uk - q0)/(q1 - uk)

    qstar <- uk
    pstar <- p0 + (p1 - p0)*(q1 - uk)/(q1 - q0)
    mstar <- (p1 - pstar)*(qstar + 0.5*(p1 - pstar)/f1) + m1

    k <- cut(p, breaks=c(pk, 1), include.lowest=TRUE, labels=FALSE)

    return(list(
        threshold = ifelse(p <= pstar[k],
            q0[k] + (p - p0[k])/f0[k],
            qstar[k] + (p - pstar[k])/f1[k]
        ),
        top_share = ifelse(p <= pstar[k],
            (pstar[k] - p)*(q0[k] + 0.5*(p - p0[k] + pstar[k] - p0[k])/f0[k]) + mstar[k],
            (p1[k] - p)*(qstar[k] + 0.5*(p - pstar[k] + p1[k] - pstar[k])/f1[k]) + m1[k]
        )/average
    ))
}

# ---------------------------------------------------------------------------- #
# Method 4: Piecewise Pareto with both thresholds and mean information
# ---------------------------------------------------------------------------- #

f <- function(a, q0, q1) {
    if (a == 0) {
        return((q1 - q0)/(log(q1) - log(q0)))
    } else if (a == 1) {
        return(q0*q1*(log(q1) - log(q0))/(q1 - q0))
    } else {
        return(a/(a - 1)*(q1^(1 - a) - q0^(1 - a))/(q1^(-a) - q0^(-a)))
    }
}

method4 <- function(p, pk, qk, mk, average) {
    # Number of thresholds
    n <- length(pk)

    # Mean inside each threshold
    uk <- -diff(mk)/diff(pk)

    # Solve an equation for the Pareto coefficient
    ak <- sapply(1:(n - 1), function(i) {
        return(uniroot(function(x) f(x, qk[i], qk[i + 1]) - uk[i],
            lower = -10,
            upper = 10,
            extendInt = "downX",
            tol = sqrt(.Machine$double.eps)
        )$root)
    })
    kk <- -ak*(pk[2:n] - pk[1:(n - 1)])/(qk[2:n]^(-ak) - qk[1:(n - 1)]^(-ak))

    # Constant Pareto coefficient in the last bracket
    bn <- mk[n]/(qk[n]*(1 - pk[n]))
    an <- bn/(bn - 1)
    ak <- c(ak, an)
    kk <- c(kk, an*(1 - pk[n])/qk[n]^(-an))

    k <- cut(p, breaks=c(pk, 1), include.lowest=TRUE, labels=FALSE)

    q <- (qk[k]^(-ak[k]) - ak[k]/kk[k]*(p - pk[k]))^(-1/ak[k])
    m <- mk[k] - kk[k]/(ak[k] - 1)*(qk[k]^(1 - ak[k]) - q^(1 - ak[k]))

    return(list(threshold=q, top_share=m/average))
}

# ---------------------------------------------------------------------------- #
# Method 5: general quadratic Lorenz curve (Villasenor & Arnold, 1989)
# ---------------------------------------------------------------------------- #

method5 <- function(p, pk, mk, average) {
    # Calculate Lorenz curve values
    lk <- 1 - mk/average

    # Define the variables from Villasenor & Arnold
    tk <- lk*(1 - lk)
    uk <- pk^2 - lk
    vk <- lk*(pk - 1)
    wk <- pk - lk

    # Fit the model (constrained version)
    fit <- lm(tk - wk ~ 0 + vk)

    # Extract parameters
    a <- 0 #coef(fit)[1]
    b <- coef(fit)[1]
    d <- 1 #coef(fit)[3]

    e     <- -(a + b + d + 1)
    alpha <- b^2 - 4*a
    beta  <- 2*b*e - 4*d

    # Estimate Lorenz curve and quantile function
    l <- suppressWarnings(0.5*(-(b*p + e) - sqrt(alpha*p^2 + beta*p + e^2)))
    q <- suppressWarnings(0.5*(-b - (beta + 2*alpha*p)/(2*sqrt((-1 - a - b - d)^2 + beta*p + alpha*p^2))))

    l[p == 0] <- 0
    l[p == 1] <- 1
    q[p == 0] <- 0
    q[p == 1] <- +Inf

    return(list(threshold = q*average, top_share = 1 - l))
}

# ---------------------------------------------------------------------------- #
# Method 6: Beta method Kakwani and Podder (1976)
# ---------------------------------------------------------------------------- #

method6 <- function(p, pk, mk, average) {
    lk <- 1 - mk/average

    rk <- (pk + lk)/sqrt(2)
    yk <- (pk - lk)/sqrt(2)

    fit <- lm(log(yk) ~ log(rk) + log(sqrt(2) - rk))

    a <- exp(coef(fit)[1])
    alpha <- coef(fit)[2]
    beta <- coef(fit)[3]

    f <- function(pi) {
        eta <- a*pi^alpha*(sqrt(2) - pi)^beta

        eta[pi <= .Machine$double.eps] <- 0
        eta[pi >= sqrt(2) - .Machine$double.eps] <- 0
        return(eta)
    }

    get_l <- function(p) {
        # Value of Lorenz curve
        sol <- uniroot(
            f = function(l) (p - l)/sqrt(2) - f((p + l)/sqrt(2)),
            lower = 0,
            upper = 1,
            extendInt = "yes"
        )
        return(sol$root)
    }

    l <- sapply(p, get_l)
    q <- sapply(p, function(p) numDeriv::grad(get_l, p))

    return(list(threshold = q*average, top_share = 1 - l))
}
