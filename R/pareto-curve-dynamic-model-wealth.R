# ---------------------------------------------------------------------------- #
# “Generalized Pareto Curves: Theory and Applications”, 2017
# Thomas Blanchet, Juliette Fournier, Thomas Piketty
# ---------------------------------------------------------------------------- #
# Generate the Pareto curve for the US distribution of wealth in 2010
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# Plot the input data
# ---------------------------------------------------------------------------- #

data_target <- subset(data_us_fr, year == 2010 & var_code == "hweal" & iso == "US")
data_target <- subset(data_target, p >= 50000 & p <= 99900)

# Generalized Pareto curve for the top tail
p <- data_target$p[data_target$p >= 50000]/1e5
b <- data_target$b[data_target$p >= 50000]

plot(p, b, type='l')

# ---------------------------------------------------------------------------- #
# Calibrate the model using variance
# ---------------------------------------------------------------------------- #

p_target <- data_target$p[data_target$p >= 80000]
q_target <- data_target$threshold[data_target$p >= 80000]
b_target <- data_target$b[data_target$p >= 80000]

q_target <- q_target/data_target$average[1]

p0_top <- 0.8
q0_top <- q_target[p_target == round(1e5*p0_top)]

stationary_dist_variance <- function(theta, x, p0_top, q0_top) {
    c1 <- theta[1]
    c2 <- theta[2]
    c3 <- theta[3]
    c4 <- theta[4]

    # Create the rational function inside the equation
    zeta <- rationalfun(
        numer = c(0, 2 + 2*c2, 0, 4*c3 + 4*c4 + 4*c2*c4, 0, 2*c3*c4 + 2*c4^2 + 2*c2*c4^2),
        denom = c(c1, 0, c2 + 2*c1*c4, 0, c3 + 2*c2*c4 + c1*c4^2, 0, c3*c4 + c2*c4^2)
    )
    # Integrate the rational function
    zeta_int <- tryCatch(int2fun(integral(zeta)), error=function(e) NA)
    if (suppressWarnings(is.na(zeta_int))) {
        return(NA)
    }
    # Non-normalized density
    pdf_nonnorm <- function(x) {
        return(exp(-zeta_int(x)))
    }
    # Calculate the constant
    cons <- tryCatch(integrate(pdf_nonnorm, lower=q0_top, upper=+Inf)$value, error=function(e) NA)
    if (is.na(cons)) {
        return(NA)
    }
    # Normalized density
    pdf <- function(x) pdf_nonnorm(x)/cons*(1 - p0_top)
    # Distribution function
    cdf <- Vectorize(function(x) p0_top + integrate(pdf, lower=q0_top, upper=x)$value)
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

theta_init <- c(0.9945929428, 5.0189544187, 0.0004407915, 0.0004620301)
fit <- optim(
    par = theta_init,
    fn = function(theta) {
        dist <- suppressWarnings(stationary_dist_variance(theta, q_target, p0_top, q0_top))
        if (length(dist) == 1 && is.na(dist)) {
            return(+Inf)
        } else {
            return(sum((dist$b - b_target)^2))
        }
    },
    method = "Nelder-Mead",
    control = list(trace=1, maxit=1000, REPORT=10)
)
theta_init <- fit$par

repeat {
    oldval <- fit$value
    fit <- optim(
        par = theta_init,
        fn = function(theta) {
            dist <- suppressWarnings(stationary_dist_variance(theta, q_target, p0_top, q0_top))
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

dist <- stationary_dist_variance(theta, q_target, p0_top, q0_top)

plot(dist$p, dist$b, type='l', col='red')
points(p_target/1e5, b_target, type='l')

# ---------------------------------------------------------------------------- #
# Simulate the calibrated model
# ---------------------------------------------------------------------------- #

c1 <- theta[1]
c2 <- theta[2]
c3 <- theta[3]
c4 <- theta[4]

# Create the rational function inside the equation
zeta <- rationalfun(
    numer = c(0, 2 + 2*c2, 0, 4*c3 + 4*c4 + 4*c2*c4, 0, 2*c3*c4 + 2*c4^2 + 2*c2*c4^2),
    denom = c(c1, 0, c2 + 2*c1*c4, 0, c3 + 2*c2*c4 + c1*c4^2, 0, c3*c4 + c2*c4^2)
)
# Integrate the rational function
zeta_int <- tryCatch(int2fun(integral(zeta)), error=function(e) NA)
if (suppressWarnings(is.na(zeta_int))) {
    return(NA)
}
# Non-normalized density
f_nonnorm <- function(x) {
    return(exp(-zeta_int(x)))
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

hist(x, breaks=5000, probability=TRUE, xlim=c(0, 5))
curve(f(x), col="red", add=TRUE)

y <- rev(sort(x))[1:(n/10)]
plot(sort(log(y)), log(1 - (1:length(y))/(length(y) + 1)))
lm(y ~ 0 + log(1 - (1:length(y))/(length(y) + 1)))

# ---------------------------------------------------------------------------- #
# Plot the volatility profile
# ---------------------------------------------------------------------------- #

x <- seq(2, 500, 0.1)
df <- data.frame(
    x = x,
    sigma = sqrt(sigma2(x))
)

dir.create("output/plots/pareto-curve-dynamic-model", showWarnings=FALSE, recursive=TRUE)
filename <- "output/plots/pareto-curve-dynamic-model/wealth-volatility.pdf"
pdf(filename, family=plot_font, width=4.5, height=3.5)
print(ggplot(df) + geom_line(aes(x=x, y=sigma)) +
    scale_x_log10(breaks=c(0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500)) +
    xlab("multiple of mean wealth") +
    ylab(expression(paste("coefficient of variation |", sigma, "(", italic(x), ")/", mu, "(", italic(x), ")|"))) +
    ggtitle("Volatility of wealth growth",
        subtitle = "calibrated to match the distribution \nof of US personal wealth in 2010") +
    theme_bw() + theme(
        legend.background = element_rect(linetype="solid", color="black", size=0.25, fill=plot_bg),
        legend.box.margin = margin(10, 10, 10, 10),
        legend.direction = "horizontal",
        plot.title=element_text(hjust=0.5),
        plot.subtitle=element_text(hjust=0.5),
        plot.background = element_rect(fill=plot_bg, color=plot_bg),
        panel.background = element_rect(fill=plot_bg),
        legend.key = element_rect(fill=plot_bg),
        text = element_text(color=plot_text_color)
    ))
dev.off()
embed_fonts(path.expand(filename))

# ---------------------------------------------------------------------------- #
# Compare the Pareto curve of the model with reality
# ---------------------------------------------------------------------------- #

dist <- stationary_dist_variance(theta, q_target, p0_top, q0_top)

df <- rbind(
    data.frame(
        p = dist$p,
        b = dist$b,
        source = rep("model", length(dist$p))
    ),
    data.frame(
        p = p_target/1e5,
        b = b_target,
        source = rep("data", length(p_target))
    )
)
df <- df[df$p <= 0.999, ]

filename <- "output/plots/pareto-curve-dynamic-model/wealth-var-stationary-dist.pdf"
pdf(filename, family=plot_font, width=4.5, height=3.5)
print(ggplot(df) + geom_line(aes(x=p, y=b, linetype=source, color=source), na.rm=TRUE) + theme_bw() +
    ggtitle("Generalized Pareto curve", subtitle="United States personal wealth in 2010") +
    ylab(expression(paste("inverted Pareto coefficient ", italic(b), "(", italic(p), ")"))) +
    xlab(expression(paste("rank ", italic(p)))) +
    theme(legend.justification = c(1, 1), legend.position = c(1, 1), legend.title = element_blank(),
        legend.background = element_rect(linetype="solid", color="black", size=0.25, fill=plot_bg),
        legend.box.margin = margin(10, 10, 10, 10),
        legend.direction = "horizontal",
        plot.title=element_text(hjust=0.5),
        plot.subtitle=element_text(hjust=0.5),
        plot.background = element_rect(fill=plot_bg, color=plot_bg),
        panel.background = element_rect(fill=plot_bg),
        legend.key = element_rect(fill=plot_bg),
        text = element_text(color=plot_text_color)
    ) +
    scale_x_continuous(limits = c(p0_top, 1), breaks=seq(0.1, 1, 0.05)) +
    scale_color_brewer("source", type="qual", palette="Set1") +
    scale_linetype_manual("source", values=c("model"="solid", "data"="longdash"))
)
dev.off()
embed_fonts(path.expand(filename))

# ---------------------------------------------------------------------------- #
# Calibrate the model using the mean
# ---------------------------------------------------------------------------- #

p_target <- data_target$p[data_target$p >= 80000]
q_target <- data_target$threshold[data_target$p >= 80000]
b_target <- data_target$b[data_target$p >= 80000]

q_target <- q_target/data_target$average[1]

p0_top <- 0.8
q0_top <- q_target[p_target == round(1e5*p0_top)]

stationary_dist_mean <- function(theta, x, p0_top, q0_top) {
    c1 <- theta[1]
    c2 <- theta[2]
    c3 <- theta[3]
    c4 <- theta[4]

    # Create the rational function inside the equation
    zeta <- rationalfun(
        numer = c(0, 2 + 2*c1, 0, -2*c2 + 2*c3 + 2*c1*c3),
        denom = c(c4, 0, 1 + c3*c4, 0, c3)
    )
    # Integrate the rational function
    zeta_int <- tryCatch(int2fun(integral(zeta)), error=function(e) NA)
    if (suppressWarnings(is.na(zeta_int))) {
        return(NA)
    }
    # Non-normalized density
    pdf_nonnorm <- function(x) {
        return(exp(-zeta_int(x)))
    }
    # Calculate the constant
    cons <- tryCatch(integrate(pdf_nonnorm, lower=q0_top, upper=+Inf)$value, error=function(e) NA)
    if (is.na(cons)) {
        return(NA)
    }
    # Normalized density
    pdf <- function(x) pdf_nonnorm(x)/cons*(1 - p0_top)
    # Distribution function
    cdf <- Vectorize(function(x) p0_top + integrate(pdf, lower=q0_top, upper=x)$value)
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

theta_init <- c(1.28887964, 0.03269871, 0.03401005, 2.57372409)
fit <- optim(
    par = theta_init,
    fn = function(theta) {
        dist <- suppressWarnings(stationary_dist_mean(theta, q_target, p0_top, q0_top))
        if (length(dist) == 1 && is.na(dist)) {
            return(+Inf)
        } else {
            return(sum((dist$b - b_target)^2))
        }
    },
    method = "Nelder-Mead",
    control = list(trace=1, maxit=1000, REPORT=10)
)
theta_init <- fit$par

repeat {
    oldval <- fit$value
    fit <- optim(
        par = theta_init,
        fn = function(theta) {
            dist <- suppressWarnings(stationary_dist_mean(theta, q_target, p0_top, q0_top))
            if (length(dist) == 1 && is.na(dist)) {
                return(+Inf)
            } else {
                return(sum((dist$b - b_target)^2))
            }
        },
        method = "Nelder-Mead",
        control = list(trace=1, maxit=10000)
    )
    if (abs((fit$value - oldval)/oldval) < 1e-10) {
        theta <- fit$par
        break
    }
    theta_init <- fit$par
}

dist <- stationary_dist_mean(theta, q_target, p0_top, q0_top)

plot(dist$p, dist$b, type='l', col='red')
points(p_target/1e5, b_target, type='l')

# ---------------------------------------------------------------------------- #
# Simulate the calibrated model
# ---------------------------------------------------------------------------- #

c1 <- theta[1]
c2 <- theta[2]
c3 <- theta[3]
c4 <- theta[4]

# Create the rational function inside the equation
zeta <- rationalfun(
    numer = c(0, 2 + 2*c1, 0, -2*c2 + 2*c3 + 2*c1*c3),
    denom = c(c4, 0, 1 + c3*c4, 0, c3)
)
# Integrate the rational function
zeta_int <- tryCatch(int2fun(integral(zeta)), error=function(e) NA)
if (suppressWarnings(is.na(zeta_int))) {
    return(NA)
}
# Non-normalized density
f_nonnorm <- function(x) {
    return(exp(-zeta_int(x)))
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
    return((c4 + x^2)/x^2)
}
mu <- function(x) {
    return(-c1 + c2*x^2/(1 + c3*x^2))
}

set.seed(19920902)

n <- 5e4
t <- 5e3

dt <- 1e-3

x <- rep(1, n)
pb <- txtProgressBar(min=0, max=t, initial=0, style=3)
for (i in 1:t) {
    sigma <- sqrt(sigma2(x))

    x_new <- x + mu(x)*dt*x + x*sqrt(dt)*sigma*rnorm(n)
    x_new[x_new <= 0] <- x[x_new <= 0]
    x <- x_new

    setTxtProgressBar(pb, i)
}
close(pb)

hist(x, breaks=5000, probability=TRUE, xlim=c(0, 5))
curve(f(x), col="red", add=TRUE)

y <- rev(sort(x))[1:(n/10)]
plot(sort(log(y)), log(1 - (1:length(y))/(length(y) + 1)))
lm(y ~ 0 + log(1 - (1:length(y))/(length(y) + 1)))

# ---------------------------------------------------------------------------- #
# Plot the mean profile
# ---------------------------------------------------------------------------- #

x <- seq(2, 500, 0.1)
df <- data.frame(
    x = x,
    mu = mu(x)
)

dir.create("output/plots/pareto-curve-dynamic-model", showWarnings=FALSE, recursive=TRUE)
filename <- "output/plots/pareto-curve-dynamic-model/wealth-mean.pdf"
pdf(filename, family=plot_font, width=4.5, height=3.5)
print(ggplot(df) + geom_line(aes(x=x, y=mu)) +
    scale_x_log10(breaks=c(0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500)) +
    xlab("multiple of mean wealth") +
    ylab(expression(paste("normalized mean ", mu, "(", italic(x), ")"))) +
    ggtitle("Mean of normalized wealth growth",
        subtitle = "calibrated to match the distribution \nof of US personal wealth in 2010") +
    theme_bw() + theme(
        legend.background = element_rect(linetype="solid", color="black", size=0.25, fill=plot_bg),
        legend.box.margin = margin(10, 10, 10, 10),
        legend.direction = "horizontal",
        plot.title=element_text(hjust=0.5),
        plot.subtitle=element_text(hjust=0.5),
        plot.background = element_rect(fill=plot_bg, color=plot_bg),
        panel.background = element_rect(fill=plot_bg),
        legend.key = element_rect(fill=plot_bg),
        text = element_text(color=plot_text_color)
    ))
dev.off()
embed_fonts(path.expand(filename))

# ---------------------------------------------------------------------------- #
# Compare the Pareto curve of the model with reality
# ---------------------------------------------------------------------------- #

dist <- stationary_dist_mean(theta, q_target, p0_top, q0_top)

df <- rbind(
    data.frame(
        p = dist$p,
        b = dist$b,
        source = rep("model", length(dist$p))
    ),
    data.frame(
        p = p_target/1e5,
        b = b_target,
        source = rep("data", length(p_target))
    )
)
df <- df[df$p <= 0.999, ]

filename <- "output/plots/pareto-curve-dynamic-model/wealth-mean-stationary-dist.pdf"
pdf(filename, family=plot_font, width=4.5, height=3.5)
print(ggplot(df) + geom_line(aes(x=p, y=b, linetype=source, color=source), na.rm=TRUE) + theme_bw() +
    ggtitle("Generalized Pareto curve", subtitle="United States personal wealth in 2010") +
    ylab(expression(paste("inverted Pareto coefficient ", italic(b), "(", italic(p), ")"))) +
    xlab(expression(paste("rank ", italic(p)))) +
    theme(legend.justification = c(1, 1), legend.position = c(1, 1), legend.title = element_blank(),
        legend.background = element_rect(linetype="solid", color="black", size=0.25, fill=plot_bg),
        legend.box.margin = margin(10, 10, 10, 10),
        legend.direction = "horizontal",
        plot.title=element_text(hjust=0.5),
        plot.subtitle=element_text(hjust=0.5),
        plot.background = element_rect(fill=plot_bg, color=plot_bg),
        panel.background = element_rect(fill=plot_bg),
        legend.key = element_rect(fill=plot_bg),
        text = element_text(color=plot_text_color)
    ) +
    scale_x_continuous(limits = c(p0_top, 1), breaks=seq(p0_top, 1, 0.05)) +
    scale_color_brewer("source", type="qual", palette="Set1") +
    scale_linetype_manual("source", values=c("model"="solid", "data"="longdash"))
)
dev.off()
embed_fonts(path.expand(filename))
