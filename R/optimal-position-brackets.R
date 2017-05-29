# ---------------------------------------------------------------------------- #
# “Generalized Pareto Curves: Theory and Applications”, 2017
# Thomas Blanchet, Juliette Fournier, Thomas Piketty
# ---------------------------------------------------------------------------- #
# Determine the optimal position of brackets in generalized Pareto
# interpolation.
# ---------------------------------------------------------------------------- #

# Calculate the median profile for phi''' for France and the United States
data_pretax <- subset(dina_data, income_type_short == "pretax")
data_pretax <- subset(data_pretax, p >= 10000 & p <= 99900)

pmin <- 0.1
pmax <- 0.999
xmin <- -log(1 - pmin)
xmax <- -log(1 - pmax)
xout <- seq(xmin, xmax, length.out=1000)

phi_d3_data <- ddply(data_pretax, c("year", "country"), function(data) {
    p <- data$p/1e5
    x <- -log(1 - p)

    phi_d3 <- sapply(xout, function(x0) lpoly(x0, x, data$phi, data$dphi, 0.05))

    return(data.frame(x=xout, phi_d3=phi_d3))
})

# Calculate median profile
phi_d3_median <- ddply(phi_d3_data, "x", function(data) {
    return(data.frame(x=data$x[1], phi_d3=median(data$phi_d3)))
})

# Smooth that profile
phi_d3_smooth <- smooth.spline(phi_d3_median$x, phi_d3_median$phi_d3, df=40)
phi_d3_fun <- function(x) predict(phi_d3_smooth, x)$y

# Function that gives the maximal error of a given placement of x1, ..., xn
max_error <- function(xk) {
    n <- length(xk)
    # Make the approximation that the maximum error is at the midpoints
    midpoints <- (xk[1:(n - 1)] + xk[2:n])/2

    if (xk[n - 1] >= xk[n]) {
        return(+Inf)
    }

    return(max(abs(interpolation_value_error_bound_noncons(midpoints, xk, phi_d3_fun))))
}

# Solve an optimization problem to find the optimal position of the thresholds
results <- list()
for (n in c(4, 5, 6, 7, 8)) {
    xk <- -log(1 - seq(pmin, pmax, length.out=n)) # seq(xmin, xmax, length.out=n)
    # Solve the optimization problem using bracket size d[k] = x[k+1] - x[k],
    # which leads to simpler constraints on the parameter space.
    dk <- diff(xk)[1:(length(xk) - 2)]

    fit <- optim(
        par = log(dk),
        fn = function(theta) {
            dk <- exp(theta)
            xk <- c(cumsum(c(xmin, dk)), xmax)
            return(max_error(xk))
        },
        method = "Nelder-Mead",
        control = list(trace=1, maxit=10000)
    )

    maxerr <- fit$value
    dk_fit <- exp(fit$par)
    xk_fit <- c(cumsum(c(xmin, dk_fit)), xmax)
    # Run the algorithm repeatedly to ensure convergence is really OK
    repeat {
        maxerr_old <- maxerr
        fit <- optim(
            par = log(dk_fit),
            fn = function(theta) {
                dk <- exp(theta)
                xk <- c(cumsum(c(xmin, dk)), xmax)
                return(max_error(xk))
            },
            method = "Nelder-Mead",
            control = list(trace=1, maxit=10000)
        )
        maxerr <- fit$value
        dk_fit <- exp(fit$par)
        xk_fit <- c(cumsum(c(xmin, dk_fit)), xmax)
        if (abs((maxerr - maxerr_old)/maxerr_old) < 1e-3) {
            break
        }
    }

    results[[n]] <- list(
        pk_fit = 1 - exp(-xk_fit),
        maxerr = maxerr
    )
    cat(paste0(rep("*", 97), collapse=""))
    cat(paste0("\nOptimal position for ", n - 1, " brackets:\n"))
    cat(paste0("--> Initial value: ", paste(sprintf("%2.2f%%", 100*(1 - exp(-xk))) , collapse=", ")), "\n")
    cat(paste0("--> Solution     : ", paste(sprintf("%2.2f%%", 100*(1 - exp(-xk_fit))) , collapse=", ")), "\n")
    cat(paste0(rep("*", 97), collapse=""))
    cat("\n")
}

# Create the LaTeX table comparing the results
filename <- "output/tables/optimal-position-brackets.tex"

sink(filename)
cat("\\begin{tabular}{cccccc}\\toprule\n")
cat("& 3 brackets & 4 brackets & 5 brackets & 6 brackets & 7 brackets \\\\ \\midrule \n")
for (i in 1:8) {
    if (i == 1) {
        cat("\\multirow{8}{*}{\\begin{tabular}[c]{@{}c@{}}optimal placement\\\\ of thresholds\\end{tabular}} ")
    }
    for (n in c(4, 5, 6, 7, 8)) {
        if (!is.na(results[[n]]$pk_fit[i])) {
            cat(sprintf("& %2.1f\\%% ", 100*results[[n]]$pk_fit[i]))
        } else {
            cat("& ")
        }
    }
    cat("\\\\ \n")
}
cat("\\midrule \n")
cat("\\begin{tabular}{@{}c}maximum relative \\\\ error on top shares\\end{tabular} & ")
cat(paste0(sapply(c(4, 5, 6, 7, 8), function(n) sprintf("%2.2f\\%%", 100*results[[n]]$maxerr)), collapse=" & "))
cat(" \\\\ \\bottomrule \n\\end{tabular}\n")
sink()
