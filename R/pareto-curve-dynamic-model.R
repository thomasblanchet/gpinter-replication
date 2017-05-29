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

stationary_dist <- function(theta, x) {
    c1 <- theta[1]
    c2 <- theta[2]
    c3 <- theta[3]
    c4 <- theta[4]

    # Create the rational function inside the equation
    zeta <- rationalfun(
        numer = c(-c1, 0, 2 + c2 - 2*c1*c4, 0, 3*c3 + 4*c4 + 2*c2*c4 - c1*c4^2, 0, c3*c4 + 2*c4^2 + c2*c4^2),
        denom = c(0, c1, 0, c2 + 2*c1*c4, 0, c3 + 2*c2*c4 + c1*c4^2, 0, c3*c4 + c2*c4^2)
    )
    # Integrate the rational function
    zeta_int <- tryCatch(int2fun(integral(zeta)), error=function(e) NA)
    if (suppressWarnings(is.na(zeta_int))) {
        return(NA)
    }
    # Non-normalized density
    pdf_nonnorm <- function(x) {
        return(exp(-zeta_int(x))/x)
    }
    # Calculate the constant
    cons <- tryCatch(integrate(pdf_nonnorm, lower=0, upper=+Inf)$value, error=function(e) NA)
    if (is.na(cons)) {
        return(NA)
    }
    # Normalized density
    pdf <- function(x) pdf_nonnorm(x)/cons
    # Distribution function
    cdf <- Vectorize(function(x) integrate(pdf, lower=0, upper=x)$value)
    # Non-normalized Lorenz curve
    lorenz <- Vectorize(function(x) integrate(function(t) t*pdf(t), lower=x, upper=+Inf)$value)
    # Calculate density, CDF, Pareto coefficients
    return(tryCatch({
        f <- pdf(x)
        p <- cdf(x)
        m <- lorenz(x)
        return(list(x=x, p=p, f=f, b=m/(1 - p)/x))
    }, error = function(e) NA))
}

# ---------------------------------------------------------------------------- #
# Calibrate the parameters of the volatility profile to match the actual
# Pareto curve
# ---------------------------------------------------------------------------- #

data_target <- subset(dina_data_us, year == 2014 & income_type_short == "labor")
data_target <- subset(data_target, p >= 10000 & p <= 99900)

theta_init <- rep(1, 4)

p_target <- c(seq(0.1, 0.9, 0.1), seq(0.91, 0.99, 0.01))
q_target <- sapply(p_target, function(p) {
    data_target$threshold[round(1e5*p) == data_target$p]
})
b_target <- sapply(p_target, function(p) {
    data_target$invpareto[round(1e5*p) == data_target$p]
})
q_target <- q_target/q_target[p_target == 0.5]

fit <- optim(
    par = theta_init,
    fn = function(theta) {
        dist <- suppressWarnings(stationary_dist(theta, q_target))
        if (length(dist) == 1 && is.na(dist)) {
            return(+Inf)
        } else {
            return(sum((dist$b - b_target)^2))
        }
    },
    method = "Nelder-Mead",
    control = list(trace=1, maxit=10000)
)
theta_init <- fit$par

repeat {
    oldval <- fit$value
    fit <- optim(
        par = theta_init,
        fn = function(theta) {
            dist <- suppressWarnings(stationary_dist(theta, q_target))
            if (length(dist) == 1 && is.na(dist)) {
                return(+Inf)
            } else {
                return(sum((dist$b - b_target)^2))
            }
        },
        method = "Nelder-Mead",
        control = list(trace=1, maxit=10000)
    )
    if (abs((fit$value - oldval)/oldval) < 1e-8) {
        theta <- fit$par
        break
    }
    theta_init <- fit$par
}


x <- seq(0.1, 100, length.out = 1000)
dist <- stationary_dist(theta, x)

b <- dist$b[dist$p >= 0.1]
p <- dist$p[dist$p >= 0.1]
p <- p[b > 0]
b <- b[b > 0]
plot(p, b, type='l', col='red')
points(data_target$p/1e5, data_target$invpareto, type='l')

# ---------------------------------------------------------------------------- #
# Simulate the calibrated model
# ---------------------------------------------------------------------------- #

c1 <- theta[1]
c2 <- theta[2]
c3 <- theta[3]
c4 <- theta[4]

# Create the rational function inside the equation
zeta <- rationalfun(
    numer = c(-c1, 0, 2 + c2 - 2*c1*c4, 0, 3*c3 + 4*c4 + 2*c2*c4 - c1*c4^2, 0, c3*c4 + 2*c4^2 + c2*c4^2),
    denom = c(0, c1, 0, c2 + 2*c1*c4, 0, c3 + 2*c2*c4 + c1*c4^2, 0, c3*c4 + c2*c4^2)
)
# Integrate the rational function
zeta_int <- tryCatch(int2fun(integral(zeta)), error=function(e) NA)
if (suppressWarnings(is.na(zeta_int))) {
    return(NA)
}
# Non-normalized density
f_nonnorm <- function(x) {
    return(exp(-zeta_int(x))/x)
}
# Calculate the constant
cons <- tryCatch(integrate(f_nonnorm, lower=0, upper=+Inf)$value, error=function(e) NA)
if (is.na(cons)) {
    return(NA)
}
# Normalized density
f <- function(x) f_nonnorm(x)/cons
# Variance of income shocks
sigma2 <- function(x) {
    return((c1 + c2*x^2)/x^2 + (c3*x^2)/(1 + c4*x^2))
}

set.seed(19920902)

n <- 5e4
t <- 5e3

dt <- 1e-3

x <- rep(1, n)
pb <- txtProgressBar(min=0, max=t, initial=0, style=3)
for (i in 1:t) {
    sigma <- sqrt(sigma2(x))

    x_new <- x - dt*x + x*sqrt(dt)*sigma*rnorm(n)
    x_new[x_new <= 0] <- x[x_new <= 0]
    x <- x_new

    setTxtProgressBar(pb, i)
}
close(pb)

hist(x, breaks=5000, probability=TRUE, xlim=c(0, 10))
curve(f(x), col="red", add=TRUE)

y <- rev(sort(x))[1:(n/10)]
plot(sort(log(y)), log(1 - (1:length(y))/(length(y) + 1)))
lm(y ~ 0 + log(1 - (1:length(y))/(length(y) + 1)))

# ---------------------------------------------------------------------------- #
# Plot the volatility profile
# ---------------------------------------------------------------------------- #

x <- seq(1, 200, 0.1)
df <- data.frame(
    x = x,
    sigma = sqrt(sigma2(x))
)

filename <- "output/plots/pareto-curve-dynamic-model/volatility.pdf"
pdf(filename, family=plot_font, width=4.5, height=3.5)
print(ggplot(df) + geom_line(aes(x=x, y=sigma)) +
    scale_x_log10(breaks=c(1, 2, 5, 10, 20, 50, 100, 200)) +
    xlab("multiple of median earnings") +
    ylab(expression(paste("coefficient of variation ", sigma, "(", italic(x), ")/", mu, "(", italic(x), ")"))) +
    ggtitle("Volatility of earnings growth",
        subtitle = "calibrated to match the distribution \nof of US labor income in 2014") +
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

x <- seq(0.1, 100, length.out = 10000)
dist <- stationary_dist(theta, x)

df <- rbind(
    data.frame(
        p = dist$p,
        b = dist$b,
        source = rep("model", length(dist$p))
    ),
    data.frame(
        p = data_target$p/1e5,
        b = data_target$invpareto,
        source = rep("data", length(data_target$p))
    )
)
df <- df[df$p >= 0.1 & df$p <= 0.999, ]

filename <- "output/plots/pareto-curve-dynamic-model/stationary-dist.pdf"
pdf(filename, family=plot_font, width=4.5, height=3.5)
print(ggplot(df) + geom_line(aes(x=p, y=b, linetype=source, color=source), na.rm=TRUE) + theme_bw() +
    ggtitle("Generalized Pareto curve", subtitle="United States labor income in 2014") +
    ylab(expression(paste("inverted Pareto coefficient ", italic(b), "(", italic(p), ")"))) +
    xlab(expression(paste("rank ", italic(p)))) +
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

# ---------------------------------------------------------------------------- #
# Discrete-time model
# ---------------------------------------------------------------------------- #

model1_invpareto <- function(p, mu1, sigma1, mu2, sigma2, n, t) {
    set.seed(19920902)

    # Check that the mean of multiplicative shocks is lower than one
    if (exp(mu1 + sigma1^2/2) >= 1) {
        return(NA)
    }

    # Simulate the model
    x <- rlnorm(n, mu2, sigma2)
    for (i in 1:t) {
        x <- rlnorm(n, mu1, sigma1)*x + rlnorm(n, mu2, sigma2)
    }

    x <- sort(x)
    b <- rev(cumsum(rev(x)))/x/(n:1)

    return(spline((1:n)/(n + 1), b, xout=p)$y)
}

# Fit the standard Kesten process
n <- 1e4
t <- 1e2

theta_init <- c(-1, 1, 1, 1)
p_target <- seq(0.1, 0.90, 0.01)
b_target <- sapply(p_target, function(p) {
    data_target$invpareto[round(1e5*p) == data_target$p]
})

library(nloptr)

fit <- nloptr(theta_init,
    eval_f = function(theta) {
        b <- model1_invpareto(p_target, theta[1], theta[2], theta[3], theta[4], n, t)
        if (length(b) == 1 && is.na(b)) {
            return(+Inf)
        } else {
            return(sum((b - b_target)^2))
        }
    },
    lb = c(-100, 0, -100, 0),
    ub = c(100, 100, 100, 100),
    opts = list(algorithm="NLOPT_GN_DIRECT", maxeval=3e3, print_level=1)
)

theta_init <- fit$solution
fit <- optim(
    par = theta_init,
    fn = function(theta) {
        b <- model1_invpareto(p_target, theta[1], theta[2], theta[3], theta[4], n, t)
        if (length(b) == 1 && is.na(b)) {
            return(+Inf)
        } else {
            return(sum((b - b_target)^2))
        }
    },
    method = "Nelder-Mead",
    control = list(trace=10, maxit=10000)
)
theta_init <- fit$par


mu1 <- theta_init[1]
sigma1 <- theta_init[2]

mu2 <- theta_init[3]
sigma2 <- theta_init[4]

n <- 1e5
t <- 1e2

x <- rlnorm(n, mu2, sigma2)

for (i in 1:t) {
    x <- rlnorm(n, mu1, sigma1)*x + rlnorm(n, mu2, sigma2)
}

x <- sort(x)
p <- (1:n)/(n+1)
b <- rev(cumsum(rev(x)))/x/(n:1)

p_out <- seq(0.1, 0.99, 0.01)

df <- data.frame(
    p = p_out,
    b_model = spline(p, b, xout=p_out)$y,
    b_data = sapply(p_out, function(p) data_target$invpareto[data_target$p == round(1e5*p)])
)
df <- df[df$p <= 0.99, ]
df <- df[df$p >= 0.10, ]

ggplot(df) +
    geom_line(aes(x=p, y=b_model), color="blue") +
    geom_line(aes(x=p, y=b_data), color="red")



