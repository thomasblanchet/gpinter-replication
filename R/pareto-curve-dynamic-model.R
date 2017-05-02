# ---------------------------------------------------------------------------- #
# “Generalized Pareto Curves: Theory and Applications”, 2017
# Thomas Blanchet, Juliette Fournier, Thomas Piketty
# ---------------------------------------------------------------------------- #
# Generate the Pareto curve for a given profile of the variance of labor
# income shocks.
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# Function to generate the Pareto of the stationary distribution for a
# random growth with volatility profile of growth shocks given by a
# rational function
# ---------------------------------------------------------------------------- #

stationary_invpareto <- function(theta, x_out) {
    # Variance of income shocks
    sigma2 <- function(x) {
        return((theta[1]*x^2 + theta[2]*x + theta[3])/
            (theta[4]*x^2 + theta[5]*x + theta[6]))
    }
    # Derivative of that function
    dsigma2 <- function(x) {
    return((-(theta[3]*(theta[5] + 2*theta[4]*x)) + theta[1]*x*(2*theta[6] + theta[5]*x) +
        theta[2]*(theta[6] - theta[4]*x^2))/(theta[6] + x*(theta[5] + theta[4]*x))^2)
    }
    # Asymptotic value of the function
    sigma2_inf <- theta[1]/theta[4]
    # Local Pareto (Zipf) exponent of the stationary distribution
    zeta <- function(x) {
        return(1 + 2/sigma2(x) + x/sigma2(x)*dsigma2(x))
    }
    zeta_inf <- 1 + 2/sigma2_inf
    # (Non-normalized) density of the stationary distribution
    pdf <- Vectorize(function(x) {
        return(x^(-zeta_inf - 1)*exp(-log(sigma2(x)) + suppressWarnings(integrate(function(t) {
            return((zeta_inf - zeta(t))/t)
        }, lower=1, upper=x)$value)))
    })
    # Calculate the constant in that density
    cons <- tryCatch(integrate(pdf, lower=1, upper=+Inf)$value, error = function(e) {
        return(NA)
    })
    if (is.na(cons)) return(NA)
    # Cumulative distribution function
    cdf <- function(x) {
        return(integrate(function(t) pdf(t)/cons, lower=1, upper=x)$value)
    }
    # Truncated mean
    truncmean <- function(x) {
        return(integrate(function(t) t*pdf(t)/cons, lower=x, upper=+Inf)$value)
    }
    # Return the inverted Pareto coefficient
    return(tryCatch({
        sapply(x_out, function(x) {
            p <- cdf(x)
            b <- truncmean(x)/(x*(1 - p))
            return(c(p, b))
        })
    }, error = function(e) {
        return(NA)
    }))
}

# ---------------------------------------------------------------------------- #
# Calibrate the parameters of the volatility profile to match the actual
# Pareto curve
# ---------------------------------------------------------------------------- #

data_target <- subset(dina_data_us, year == 2014 & income_type_short == "labor")
data_target <- subset(data_target, p >= 10000 & p <= 99900)

theta_init <- c(
    191.239293061311,
    -941.879895132741,
    27798.5469931347,
    51.4653781154658,
    768.100007149152,
    457.513214071661
)
p_target <- c(seq(0.1, 0.9, 0.1), seq(0.91, 0.99, 0.01))
q_target <- sapply(p_target, function(p) {
    data_target$threshold[round(1e5*p) == data_target$p]
})
b_target <- sapply(p_target, function(p) {
    data_target$invpareto[round(1e5*p) == data_target$p]
})
q_target <- q_target/q_target[1]

fit <- optim(
    par = theta_init,
    fn = function(theta) {
        b <- stationary_invpareto(theta, q_target)
        if (length(b) == 1 && is.na(b)) {
            return(+Inf)
        } else {
            return(sum((b[2, ] - b_target)^2))
        }
    },
    method = "Nelder-Mead",
    control = list(trace=1, maxit=10000)
)
theta <- fit$par

# ---------------------------------------------------------------------------- #
# Plot the volatility profile
# ---------------------------------------------------------------------------- #

x <- seq(1, 1000, 0.1)
df <- data.frame(
    x = x,
    sigma = (theta[1]*x^2 + theta[2]*x + theta[3])/(theta[4]*x^2 + theta[5]*x + theta[6])
)

filename <- "output/plots/pareto-curve-dynamic-model/volatility.pdf"
pdf(filename, family=plot_font, width=4.5, height=3.5)
print(ggplot(df) + geom_line(aes(x=x, y=sigma)) +
    scale_x_log10(breaks=c(1, 2, 5, 10, 20, 50, 100, 200, 500, 1000)) +
    xlab("normalized earnings") + ylab("standard deviation") +
    ggtitle("Volatility of earnings growth",
        subtitle = "calibrated to match the distribution \nof the top 90% of US labor income in 2014") +
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

# ---------------------------------------------------------------------------- #
# Compare the Pareto curve of the model with reality
# ---------------------------------------------------------------------------- #

b_model <- stationary_invpareto(theta, seq(1, 100, 0.1))

df <- rbind(
    data.frame(
        p = 0.1 + 0.9*b_model[1, ],
        b = b_model[2, ],
        source = rep("model", ncol(b_model))
    ),
    data.frame(
        p = data_target$p/1e5,
        b = data_target$invpareto,
        source = rep("data", length(data_target$p))
    )
)

filename <- "output/plots/pareto-curve-dynamic-model/stationary-dist.pdf"
pdf(filename, family=plot_font, width=4.5, height=3.5)
print(ggplot(df) + geom_line(aes(x=p, y=b, linetype=source, color=source), na.rm=TRUE) + theme_bw() +
    ggtitle("Generalized Pareto curve", subtitle="top 90% of US labor income in 2014") +
    ylab("inverted Pareto coefficient") + xlab("fractile") +
    theme(legend.justification = c(1, 1), legend.position = c(1, 1), legend.title = element_blank(),
        legend.background = element_rect(linetype="solid", color="black", size=0.25),
        legend.box.margin = margin(10, 10, 10, 10),
        legend.direction = "horizontal",
        plot.title=element_text(hjust=0.5),
        plot.subtitle=element_text(hjust=0.5),
        plot.background = element_rect(fill=plot_bg, color=plot_bg),
        panel.background = element_rect(fill=plot_bg),
        legend.key = element_rect(fill=plot_bg),
        text = element_text(color=plot_text_color)
    ) +
    scale_x_continuous(limits = c(0.1, 1), breaks=seq(0.1, 1, 0.1)) +
    scale_color_brewer("source", type="qual", palette="Set1") +
    scale_linetype_manual("source", values=c("model"="solid", "data"="longdash"))
)
dev.off()
embed_fonts(path.expand(filename))

