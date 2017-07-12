# ---------------------------------------------------------------------------- #
# “Generalized Pareto Curves: Theory and Applications”, 2017
# Thomas Blanchet, Juliette Fournier, Thomas Piketty
# ---------------------------------------------------------------------------- #
# Estimate and plot the sampling error for the US distribution of labor income
# in 1962 (which has finite variance).
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# Function for estimating the mean absolute sampling error
# ---------------------------------------------------------------------------- #

# Spline basis functions
h00 <- function(x) {
    return(1 - 10*x^3 + 15*x^4 - 6*x^5)
}
h01 <- function(x) {
    return(10*x^3 - 15*x^4 + 6*x^5)
}
h10 <- function(x) {
    return(x - 6*x^3 + 8*x^4 - 3*x^5)
}
h11 <- function(x) {
    return(-4*x^3 + 7*x^4 - 3*x^5)
}
h20 <- function(x) {
    return(x^2/2 - (3*x^3)/2 + (3*x^4)/2 - x^5/2)
}
h21 <- function(x) {
    return(x^3/2 - x^4 + x^5/2)
}

# First derivatives
h00d1 <- function(x) {
    return(-30*x^2 + 60*x^3 - 30*x^4)
}
h01d1 <- function(x) {
    return(30*x^2 - 60*x^3 + 30*x^4)
}
h10d1 <- function(x) {
    return(1 - 18*x^2 + 32*x^3 - 15*x^4)
}
h11d1 <- function(x) {
    return(-12*x^2 + 28*x^3 - 15*x^4)
}
h20d1 <- function(x) {
    return(x - (9*x^2)/2 + 6*x^3 - (5*x^4)/2)
}
h21d1 <- function(x) {
    return((3*x^2)/2 - 4*x^3 + (5*x^4)/2)
}

sampling_error_finite <- function(x, pk, samplesize, dist) {
    p <- 1 - exp(-x)
    k <- cut(p, pk, labels=FALSE, include.lowest=TRUE)
    if (is.na(k)) {
        return(c(NA, NA))
    }
    pk <- sort(c(pk, p))
    n <- length(pk)

    xk <- -log(1 - pk)
    mk <- (1 - pk)*top_average(dist, pk)
    qk <- as.integer(fitted_quantile(dist, pk))

    muk <- vector(mode="numeric", length=n)
    for (i in 1:(n - 1)) {
        muk[i] <- integrate(function(x) {
            return(fitted_quantile(dist, x))
        }, lower=pk[i], upper=pk[i+1])$value
    }
    muk[n] <- integrate(function(x) {
        return(fitted_quantile(dist, x))
    }, lower=pk[n], upper=1)$value

    sigmak <- vector(mode="numeric", length=n)
    for (i in 1:(n - 1)) {
        sigmak[i] <- integrate(function(x) {
            return(fitted_quantile(dist, x)^2)
        }, lower=pk[i], upper=pk[i+1])$value - muk[i]^2
    }
    sigmak[n] <- integrate(function(x) {
        return(fitted_quantile(dist, x)^2)
    }, lower=pk[n], upper=1)$value - muk[n]^2

    M01 <- matrix(0, nrow=n, ncol=n)
    for (i in 1:n) {
        for (j in 1:n) {
            if (i == j) {
                M01[i, j] <- sigmak[i]
            } else {
                M01[i, j] <- -muk[i]*muk[j]
            }
        }
    }

    M02 <- matrix(0, nrow=n, ncol=n)
    for (i in 1:n) {
        for (j in 1:n) {
            M02[i, j] <- pk[min(i, j)]*(1 - pk[max(i, j)])
        }
    }

    M03 <- matrix(nrow=n, ncol=n)
    for (i in 1:n) {
        for (j in 1:n) {
            if (i >= j) {
                M03[i, j] <- -pk[j]*muk[i]
            } else {
                M03[i, j] <- (1 - pk[j])*muk[i]
            }
        }
    }

    M04 <- diag(1/fitted_density(dist, qk))

    M05 <- matrix(0, nrow=n, ncol=n)
    for (i in 1:n) {
        for (j in 1:n) {
            if (i == j) {
                M05[i, j] <- qk[i]
            } else if (j == i + 1) {
                M05[i, j] <- -qk[j]
            }
        }
    }

    M06 <- rbind(
        cbind(diag(n), M05),
        cbind(matrix(0, nrow=n, ncol=n), -M04)
    )

    M07 <- rbind(
        cbind(M01, M03),
        cbind(t(M03), M02)
    )

    ham <- M06 %*% M07 %*% t(M06)

    M01 <- matrix(0, nrow=2*n, ncol=2*n)
    for (i in 1:(2*n)) {
        for (j in 1:(2*n)) {
            if (i <= n && j <= n) {
                if (i <= j) {
                    M01[i, j] <- 1
                }
            } else {
                if (i == j) {
                    M01[i, j] <- 1
                }
            }
        }
    }

    M02 <- matrix(0, nrow=2*n, ncol=2*n)
    for (i in 1:(2*n)) {
        for (j in 1:(2*n)) {
            if (i <= n && j <= n) {
                if (i == j) {
                    M02[i, j] <- -1/mk[i]
                }
            } else if (i > n && j <= n) {
                if (i - n == j) {
                    M02[i, j] <- -exp(-xk[j])*qk[j]/mk[j]^2
                }
            } else if (i > n && j > n) {
                if (i == j) {
                    M02[i, j] <- exp(-xk[j - n])/mk[j - n]
                }
            }
        }
    }

    M03 <- matrix(0, nrow=2*n - 3, ncol=2*n - 2)
    for (i in 1:(2*n - 3)) {
        for (j in 1:(2*n - 2)) {
            if (i <= n - 2) {
                if (i == j) {
                    M03[i, j] <- -1
                } else if (i + 1 == j) {
                    M03[i, j] <- 1
                }
            } else {
                if (i + 1 == j) {
                    M03[i, j] <- 1
                }
            }
        }
    }

    M03 <- cbind(
        M03[, 1:k],
        0,
        M03[, (k + 1):(n - 1 + k)],
        0,
        M03[, (n + k):(2*n - 2)]
    )

    M04 <- matrix(0, nrow=n-1, ncol=n-1)
    xk_ <- xk[-(k + 1)]
    for (i in 1:(n - 1)) {
        for (j in 1:(n - 1)) {
            if (i == 1 && j == 1) {
                M04[i, j] <- 9/(xk_[2] - xk_[1])
            } else if (i == n - 1 && j == n - 1) {
                M04[i, j] <- 1
            } else if (i == n - 1) {
                M04[i, j] <- 0
            } else if (i == j) {
                M04[i, j] <- 9/(xk_[i] - xk_[i - 1]) + 9/(xk_[i + 1] - xk_[i])
            } else if (i + 1 == j) {
                M04[i, j] <- -3/(xk_[i + 1] - xk_[i])
            } else if (i == j + 1) {
                M04[i, j] <- -3/(xk_[j + 1] - xk_[j])
            }
        }
    }

    M05 <- matrix(0, nrow=n-1, ncol=n-2)
    for (i in 1:(n - 1)) {
        for (j in 1:(n - 2)) {
            if (i == j) {
                M05[i, j] <- 60/(xk_[j + 1] - xk_[j])^3
            } else if (i == j + 1 && i <= n - 2) {
                M05[i, j] <- -60/(xk_[i] - xk_[i - 1])^3
            }
        }
    }

    M06 <- matrix(0, nrow=n-1, ncol=n-1)
    for (i in 1:(n-1)) {
        for (j in 1:(n-1)) {
            if (i == 1 && j == 1) {
                M06[i, j] <- -36/(xk_[2] - xk_[1])^2
            } else if (i == 1 && j == 2) {
                M06[i, j] <- -24/(xk_[2] - xk_[1])^2
            } else if (i == j && i <= n - 2) {
                M06[i, j] <- 36/(xk_[i] - xk_[i - 1])^2 - 36/(xk_[i + 1] - xk_[i])^2
            } else if (i == j + 1 && i <= n - 2) {
                M06[i, j] <- 24/(xk_[i] - xk_[i - 1])^2
            } else if (i + 1 == j && i <= n - 2) {
                M06[i, j] <- -24/(xk_[j] - xk_[j - 1])^2
            } else if (i == n - 1 && j == n - 2) {
                M06[i, j] <- -1/(xk_[n - 1] - xk_[n - 2])
            } else if (i == n - 1 && j == n - 1) {
                M06[i, j] <- 1/(xk_[n - 1] - xk_[n - 2])
            }
        }
    }

    M08 <- solve(M04) %*% cbind(M05, M06) %*% M03

    M09 <- matrix(0, nrow=8, ncol=2*n)
    M09[1, k + 0] <- 1
    M09[2, k + 1] <- 1
    M09[3, k + 2] <- 1
    M09[4, n + k + 0] <- 1
    M09[5, n + k + 1] <- 1
    M09[6, n + k + 2] <- 1
    M09[7, ] <- M08[k, ]
    M09[8, ] <- M08[k + 1, ]

    x0 <- xk_[k]
    x1 <- xk_[k+1]
    t <- (x - x0)/(x1 - x0)
    u <- (x1 - x0)
    M10 <- rbind(
        c(h00(t), -1, h01(t), u*h10(t), 0, u*h11(t), u^2*h20(t), u^2*h21(t)),
        c(h00d1(t)/u, 0, h01d1(t)/u, h10d1(t), -1, h11d1(t), u*h20d1(t), u*h21d1(t))
    )

    bread <- M10 %*% M09 %*% M02 %*% M01

    sigma <- sqrt(diag(bread %*% ham %*% t(bread))/samplesize)

    return(sqrt(2/pi)*sigma)
}

# ---------------------------------------------------------------------------- #
# Application to US labor income, 1962
# ---------------------------------------------------------------------------- #

data_us_labor <- subset(data_us_fr, year == 1970 & var_code == "pllin" & iso == "US")

# Interpolate the distribution (for estimates of the density, etc.)
data_us_labor$bracket <- cut(data_us_labor$p,
    breaks = round(1e5*c(0, 0.1, 0.5, 0.9, 0.99, 1)),
    include.lowest = TRUE,
    labels = FALSE,
    right = FALSE
)
short_tab <- ddply(data_us_labor, "bracket", function(data) {
    return(data.frame(
        p = min(data$p)/1e5,
        threshold = data$threshold[which.min(data$p)],
        top_share = data$top_share[which.min(data$p)]
    ))
})
average <- data_us_labor$average[1]
# Remove the first bracket, which includes zero or negative values
short_tab <- short_tab[-1, ]
short_tab$m <- average*short_tab$top_share

# Generalized Pareto interpolation
dist <- tabulation_fit(short_tab$p, short_tab$threshold,
    average, topshare=short_tab$top_share)

# Estimate sampling error
n <- data_us_labor$population[1]
pk <- short_tab$p
xk <- -log(1 - pk)
x_out <- seq(-log(1 - 0.1), -log(1 - 0.99), length.out=1000)
sampling_error <- pbsapply(x_out, function(x) sampling_error_finite(x, pk, n, dist))
sampling_error <- data.frame(
    x    = x_out,
    phi  = as.vector(sampling_error[1, ]),
    dphi = as.vector(sampling_error[2, ])
)

# Plot
dir.create("output/plots/sampling-error", showWarnings=FALSE, recursive=TRUE)
filename <- "output/plots/sampling-error/error-phi-finite-variance.pdf"
pdf(filename, family=plot_font, width=4.5, height=3.5)
print(ggplot(data=sampling_error) + geom_line(aes(x=x, y=phi), na.rm=TRUE) +
    geom_vline(xintercept=xk, linetype="dashed") +
    annotate("text", label="paste(italic(p) == 10, '%')", x=xk[1] + 0.16, y=2.8e-5, angle=90, parse=TRUE) +
    annotate("text", label="paste(italic(p) == 50, '%')", x=xk[2] + 0.16, y=2.8e-5, angle=90, parse=TRUE) +
    annotate("text", label="paste(italic(p) == 90, '%')", x=xk[3] + 0.16, y=2.8e-5, angle=90, parse=TRUE) +
    annotate("text", label="paste(italic(p) == 99, '%')", x=xk[4] - 0.16, y=2.8e-5, angle=90, parse=TRUE) +
    xlab(expression(paste(italic(x) == -log, "(", 1 - italic(p), ")"))) +
    scale_y_continuous(name="mean absolute error", labels=fancy_scientific) +
    ggtitle(expression(paste(phi1, '(', italic(x), ") for US labor income (1970)")),
        subtitle=expression("for a tabulation with"~italic(p)~"= 10%, 50%, 90% and 99%")) +
    theme_bw() + theme(
        plot.title = element_text(hjust=0.5),
        plot.subtitle = element_text(hjust=0.5)),
        plot.background = element_rect(fill=plot_bg, color=plot_bg),
        panel.background = element_rect(fill=plot_bg),
        legend.key = element_rect(fill=plot_bg),
        text = element_text(color=plot_text_color)
    )
dev.off()
embed_fonts(path.expand(filename))

filename <- "output/plots/sampling-error/error-deriv-phi-finite-variance.pdf"
pdf(filename, family=plot_font, width=4.5, height=3.5)
print(ggplot(data=sampling_error) + geom_line(aes(x=x, y=dphi), na.rm=TRUE) +
    geom_vline(xintercept=xk, linetype="dashed") +
    annotate("text", label="paste(italic(p) == 10, '%')", x=xk[1] + 0.16, y=8e-5, angle=90, parse=TRUE) +
    annotate("text", label="paste(italic(p) == 50, '%')", x=xk[2] + 0.16, y=8e-5, angle=90, parse=TRUE) +
    annotate("text", label="paste(italic(p) == 90, '%')", x=xk[3] + 0.16, y=8e-5, angle=90, parse=TRUE) +
    annotate("text", label="paste(italic(p) == 99, '%')", x=xk[4] - 0.16, y=8e-5, angle=90, parse=TRUE) +
    xlab(expression(paste(italic(x) == -log, "(", 1 - italic(p), ")"))) +
    scale_y_continuous(name="mean absolute error", labels=fancy_scientific, limits=c(0, 9e-5)) +
    ggtitle(expression(paste(phi1, "'(", italic(x), ") for US labor income (1970)")),
        subtitle=expression("for a tabulation with"~italic(p)~"= 10%, 50%, 90% and 99%")) +
    theme_bw() + theme(
        plot.title = element_text(hjust=0.5),
        plot.subtitle = element_text(hjust=0.5)),
        plot.background = element_rect(fill=plot_bg, color=plot_bg),
        panel.background = element_rect(fill=plot_bg),
        legend.key = element_rect(fill=plot_bg),
        text = element_text(color=plot_text_color)
    )
dev.off()
embed_fonts(path.expand(filename))

