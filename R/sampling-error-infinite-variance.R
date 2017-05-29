# ---------------------------------------------------------------------------- #
# “Generalized Pareto Curves: Theory and Applications”, 2017
# Thomas Blanchet, Juliette Fournier, Thomas Piketty
# ---------------------------------------------------------------------------- #
# Estimate and plot the sampling error for the US distribution of capital income
# in 1962 (which has infinite variance).
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# Function for estimating the mean absolute sampling error
# ---------------------------------------------------------------------------- #

# The stable distribution
dstable <- function(x, alpha) {
    param <- setParam(alpha=alpha, location=0, logscale=0, pm="M")
    return(dEstable(x, param))
}
pstable <- function(x, alpha) {
    param <- setParam(alpha=alpha, location=0, logscale=0, pm="M")
    return(pEstable(x, param))
}
qstable <- function(x, alpha) {
    param <- setParam(alpha=alpha, location=0, logscale=0, pm="M")
    return(qEstable(x, param))
}
mean_abs_error_stable <- function(alpha) {
    param <- setParam(alpha=alpha, location=0, logscale=0, pm="M")

    I1 <- integrate(function(x) {
        return(-x*dEstable(x, param))
    }, lower=-Inf, upper=0)$value

    I2 <- integrate(function(x) {
        return(x*dEstable(x, param))
    }, lower=0, upper=+Inf)$value

    return(I1 + I2)
}

sampling_error_infinite <- function(x, pk, samplesize, dist) {
    p <- 1 - exp(-x)
    k <- cut(p, pk, labels=FALSE, include.lowest=TRUE)
    if (is.na(k)) {
        return(c(NA, NA))
    }
    pk <- sort(c(pk, p))
    n <- length(pk)

    xk <- -log(1 - pk)
    qk <- fitted_quantile(dist, pk)
    mk <- threshold_share(dist, qk)*dist$average

    cons <- (1 - pk[n])*(dist$xi_top/dist$sigma_top)^(-1/dist$xi_top)
    alpha <- 1/dist$xi_top

    g <- (pi*cons/(2*gamma(alpha)*sin(alpha*pi/2)))^(1/alpha)

    S <- matrix(0, nrow=2*n, ncol=1)
    S[n, 1] <- mean_abs_error_stable(alpha)

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

    return(abs(samplesize^(1/alpha - 1) * g * M10 %*% M09 %*% M02 %*% M01 %*% S))
}

# ---------------------------------------------------------------------------- #
# Application to US capital income, 1962
# ---------------------------------------------------------------------------- #

data_us_capital <- subset(dina_data_us, year == 1962 & income_type_short == "capital")

# Interpolate the distribution (for estimates of the density, etc.)
data_us_capital$bracket <- cut(data_us_capital$p,
    breaks = round(1e5*c(0, 0.1, 0.5, 0.9, 0.99, 1)),
    include.lowest = TRUE,
    labels = FALSE,
    right = FALSE
)
short_tab <- ddply(data_us_capital, "bracket", function(data) {
    return(data.frame(
        p = min(data$p)/1e5,
        threshold = data$threshold[which.min(data$p)],
        topshare = data$topshare[which.min(data$p)]
    ))
})
average <- data_us_capital$average[1]
# Remove the first bracket, which includes zero or negative values
short_tab <- short_tab[-1, ]
short_tab$m <- average*short_tab$topshare

# Generalized Pareto interpolation
dist <- tabulation_fit(short_tab$p, short_tab$threshold,
    average, topshare=short_tab$topshare)

# Estimate sampling error
n <- data_us_labor$population[1]
pk <- short_tab$p
xk <- -log(1 - pk)
x_out <- seq(-log(1 - 0.1), -log(1 - 0.99), length.out=1000)
sampling_error <- pbsapply(x_out, function(x) sampling_error_infinite(x, pk, n, dist))
sampling_error <- data.frame(
    x    = x_out,
    phi  = as.vector(sampling_error[1, ]),
    dphi = as.vector(sampling_error[2, ])
)

# Plot
filename <- "output/plots/sampling-error/error-phi-infinite-variance.pdf"
pdf(filename, family=plot_font, width=4.5, height=3.5)
print(ggplot(data=sampling_error) + geom_line(aes(x=x, y=phi), na.rm=TRUE) +
    geom_vline(xintercept=xk, linetype="dashed") +
    annotate("text", label="paste(italic(p) == 10, '%')", x=xk[1] + 0.16, y=1.8e-5, angle=90, parse=TRUE) +
    annotate("text", label="paste(italic(p) == 50, '%')", x=xk[2] + 0.16, y=1.8e-5, angle=90, parse=TRUE) +
    annotate("text", label="paste(italic(p) == 90, '%')", x=xk[3] + 0.16, y=1.8e-5, angle=90, parse=TRUE) +
    annotate("text", label="paste(italic(p) == 99, '%')", x=xk[4] - 0.16, y=1.8e-5, angle=90, parse=TRUE) +
    xlab(expression(paste(italic(x) == -log, "(", 1 - italic(p), ")"))) +
    scale_y_continuous(name="mean absolute error", labels=fancy_scientific) +
    ggtitle(expression(paste(phi1, '(', italic(x), ") for US capital income (1962)")),
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

filename <- "output/plots/sampling-error/error-deriv-phi-infinite-variance.pdf"
pdf(filename, family=plot_font, width=4.5, height=3.5)
print(ggplot(data=sampling_error) + geom_line(aes(x=x, y=dphi), na.rm=TRUE) +
    geom_vline(xintercept=xk, linetype="dashed") +
    annotate("text", label="paste(italic(p) == 10, '%')", x=xk[1] + 0.16, y=4.5e-5, angle=90, parse=TRUE) +
    annotate("text", label="paste(italic(p) == 50, '%')", x=xk[2] + 0.16, y=4.5e-5, angle=90, parse=TRUE) +
    annotate("text", label="paste(italic(p) == 90, '%')", x=xk[3] + 0.16, y=4.5e-5, angle=90, parse=TRUE) +
    annotate("text", label="paste(italic(p) == 99, '%')", x=xk[4] - 0.16, y=4.5e-5, angle=90, parse=TRUE) +
    xlab(expression(paste(italic(x) == -log, "(", 1 - italic(p), ")"))) +
    scale_y_continuous(name="mean absolute error", labels=fancy_scientific, limits=c(0, 5e-5)) +
    ggtitle(expression(paste(phi1, "'(", italic(x), ") for US capital income (1962)")),
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


