# ---------------------------------------------------------------------------- #
# “Generalized Pareto Curves: Theory and Applications”, 2017
# Thomas Blanchet, Juliette Fournier, Thomas Piketty
# ---------------------------------------------------------------------------- #
# Estimate and plot phi''' in France and the United States.
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

# Estimate phi'''
phi_d3_estim <- ddply(dina_data, c("iso", "country", "year", "income_type", "income_type_short"), function(data) {
    data <- data[data$p >= 10000 & data$p <= 99000, ]

    p <- data$p/1e5
    x <- -log(1 - data$p/1e5)

    x_out <- seq(-log(1 - min(p)), -log(1 - max(p)), length.out=1000)
    d3_phi <- sapply(x_out, function(x0) lpoly(x0, x, data$phi, data$dphi, 0.05))

    return(data.frame(x=x_out, d3_phi=d3_phi))
})

# Summarize the estimates
phi_d3_summ <- ddply(phi_d3_estim, c("iso", "country", "income_type", "income_type_short", "x"), function(data) {
    p1 <- quantile(data$d3_phi, 0.1)
    p9 <- quantile(data$d3_phi, 0.9)
    mean <- mean(data$d3_phi)
    median <- median(data$d3_phi)
    year_min <- min(data$year)
    year_max <- max(data$year)
    return(data.frame(p1=p1, p9=p9, mean=mean, median=median, year_min=year_min, year_max=year_max))
})

# Plot the function phi'''
d_ply(phi_d3_summ, c("country", "income_type"), function(data) {
    country           <- data$country[1]
    income_type       <- data$income_type[1]
    year_min          <- data$year_min[1]
    year_max          <- data$year_max[1]
    iso               <- data$iso[1]
    income_type_short <- data$income_type_short[1]

    filename <- paste0("output/plots/estim-deriv3-phi/estim-phi-d3-", iso, "-", income_type_short, ".pdf")
    pdf(filename, family="CM Roman", width=4.5, height=3.5)
    print(ggplot(data) +
        geom_line(aes(x=x, y=p1), color="grey", na.rm=TRUE) +
        geom_line(aes(x=x, y=p9), color="grey", na.rm=TRUE) +
        geom_line(aes(x=x, y=median), na.rm=TRUE) +
        geom_vline(xintercept=c(-log(1 - 0.5), -log(1 - 0.9), -log(1 - 0.99)), linetype="dashed") +
        annotate("text", label="paste(italic(p) == 50, '%')", x=-log(1 - 0.50) + 0.16, y=-5.2, angle=90, parse=TRUE) +
        annotate("text", label="paste(italic(p) == 90, '%')", x=-log(1 - 0.90) + 0.16, y=-5.2, angle=90, parse=TRUE) +
        annotate("text", label="paste(italic(p) == 99, '%')", x=-log(1 - 0.99) - 0.16, y=-5.2, angle=90, parse=TRUE) +
        ggtitle(paste0(country, ", ", year_min, "-", year_max)) +
        xlab(expression(paste(italic(x) == -log, "(", 1 - italic(p), ")"))) +
        ylab(expression(paste(phi1, "'''", "(", italic(x), ")"))) +
        ylim(c(-6, 0.5)) + theme_bw() +
        theme(plot.title = element_text(hjust = 0.5)))
    dev.off()
    embed_fonts(path.expand(filename))
})

