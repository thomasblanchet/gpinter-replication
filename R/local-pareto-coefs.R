# ---------------------------------------------------------------------------- #
# “Generalized Pareto Curves: Theory and Applications”, 2017
# Thomas Blanchet, Juliette Fournier, Thomas Piketty
# ---------------------------------------------------------------------------- #
# Compare the different concepts of local Pareto coefficients
# ---------------------------------------------------------------------------- #

data_2010 <- subset(data_us_fr, year == 2010 & var_code == "ptinc" & (iso == "FR" | iso == "US"))
data_2010 <- subset(data_2010, p >= 30000)
data_2010 <- subset(data_2010, p %% 1000 == 0)
data_2010$p <- data_2010$p/1e5

# Fit model for the US
fit_us <- lm(dphi ~ p + I(p^2/2) + I(p^3/6) + I(p^4/24) + I(p^5/120), data=subset(data_2010, iso == "US"))

p0 <- 0.3
q0 <- subset(data_2010, p == p0 & iso == "US")$threshold

dphi_fit <- function(p, fit) {
    a0 <- coef(fit)[1]
    a1 <- coef(fit)[2]
    a2 <- coef(fit)[3]
    a3 <- coef(fit)[4]
    a4 <- coef(fit)[5]
    a5 <- coef(fit)[6]

    return(a0 + p*a1 + p^2*a2/2 + p^3*a3/6 + p^4*a4/24 + p^5*a5/120)
}

dphi_deriv_fit <- function(p, fit) {
    a1 <- coef(fit)[2]
    a2 <- coef(fit)[3]
    a3 <- coef(fit)[4]
    a4 <- coef(fit)[5]
    a5 <- coef(fit)[6]

    return(a1 + p*a2 + p^2*a3/2 + p^3*a4/6 + p^4*a5/24)
}

q_fit <- function(p, fit, p0, q0) {
    return(sapply(p, function(p) {
        v1 <- integrate(function(t) {
            (1 - dphi_fit(t, fit))/(1 - t)
        }, lower=p0, upper=p)$value
        v2 <- dphi_fit(p, fit)/dphi_fit(p0, fit)
        return(q0*v2*exp(v1))
    }))
}

q_deriv_fit <- function(p, fit, p0, q0) {
    return(sapply(p, function(p) {
        v1 <- integrate(function(t) {
            (1 - dphi_fit(t, fit))/(1 - t)
        }, lower=p0, upper=p)$value
        v1_deriv <- (1 - dphi_fit(p, fit))/(1 - p)

        v2 <- dphi_fit(p, fit)/dphi_fit(p0, fit)
        v2_deriv <- dphi_deriv_fit(p, fit)/dphi_fit(p0, fit)

        return(q0*(v2_deriv*exp(v1) + v2*v1_deriv*exp(v1)))
    }))
}

plot(p, q_deriv_fit(p, fit_us, p0, q0), log='y')
q_fit(0.6, fit_us, p0, q0)
q_deriv_fit(0.6, fit_us, p0, q0)

alpha0_fit <- function(p, fit) {
    return(1/(1 - dphi_fit(p, fit)))
}

alpha0_deriv_fit <- function(p, fit, p0, q0) {
    return(dphi_deriv_fit(p, fit)/(1 - dphi_fit(p, fit))^2/q_deriv_fit(p, fit, p0, q0))
}

alpha1_fit <- function(p, fit, p0, q0) {
    return(alpha0_fit(p, fit) - q_fit(p, fit, p0, q0)*alpha0_deriv_fit(p, fit, p0, q0)/(alpha0_fit(p, fit) - 1))
}

alpha1_deriv_fit <- function(p, fit, p0, q0) {
    return(grad(function(t) alpha1_fit(t, fit, p0, q0), p)/q_deriv_fit(p, fit, p0, q0))
}

alpha2_fit <- function(p, fit, p0, q0) {
    return(alpha1_fit(p, fit, p0, q0) - q_fit(p, fit, p0, q0)*alpha1_deriv_fit(p, fit, p0, q0)/alpha1_fit(p, fit, p0, q0))
}

b1_fit <- function(p, fit, p0, q0) {
    return(alpha1_fit(p, fit, p0, q0)/(alpha1_fit(p, fit, p0, q0) - 1))
}

p <- c(seq(0.5, 0.99, 0.01), seq(0.991, 0.999, 0.001), seq(0.9991, 0.9999, 0.0001))
df_plot <- data.frame(
    p = p,
    alpha0 = alpha0_fit(p, fit_us),
    alpha1 = alpha1_fit(p, fit_us, p0, q0),
    alpha2 = alpha2_fit(p, fit_us, p0, q0)
)

dir.create(paste0("output/plots/local-pareto-exponent"),
    showWarnings=FALSE, recursive=TRUE)
filename <- paste0("output/plots/local-pareto-exponent/alpha0.pdf")
pdf(filename, family=plot_font, width=3, height=3.5)
print(ggplot(df_plot) +
    geom_line(aes(x=p, y=alpha0)) +
    xlab(expression(paste("rank ", italic(p)))) +
    ylab(expression(paste("local Pareto exponent ", alpha[0]))) +
    ggtitle(expression(alpha[0])) +
    theme_bw() + theme(
        legend.justification = c(1, 1),
        legend.position = c(1, 1),
        legend.background = element_rect(
            linetype = "solid",
            color    = plot_text_color,
            size     = 0.3,
            fill     = plot_bg
        ),
        legend.title = element_blank(),
        legend.box.margin = margin(10, 10, 10, 10),
        legend.direction = "horizontal",
        plot.title=element_text(hjust=0.5),
        plot.subtitle=element_text(hjust=0.5),
        plot.background = element_rect(fill=plot_bg, color=plot_bg),
        panel.background = element_rect(fill=plot_bg),
        legend.key = element_rect(fill=plot_bg),
        text = element_text(color=plot_text_color)
    )
)
dev.off()
embed_fonts(path.expand(filename))


dir.create(paste0("output/plots/local-pareto-exponent"),
    showWarnings=FALSE, recursive=TRUE)
filename <- paste0("output/plots/local-pareto-exponent/alpha1.pdf")
pdf(filename, family=plot_font, width=3, height=3.5)
print(ggplot(df_plot) +
    geom_line(aes(x=p, y=alpha1)) +
    xlab(expression(paste("rank ", italic(p)))) +
    ylab(expression(paste("local Pareto exponent ", alpha[1]))) +
    ggtitle(expression(alpha[1])) +
    theme_bw() + theme(
        legend.justification = c(1, 1),
        legend.position = c(1, 1),
        legend.background = element_rect(
            linetype = "solid",
            color    = plot_text_color,
            size     = 0.3,
            fill     = plot_bg
        ),
        legend.title = element_blank(),
        legend.box.margin = margin(10, 10, 10, 10),
        legend.direction = "horizontal",
        plot.title=element_text(hjust=0.5),
        plot.subtitle=element_text(hjust=0.5),
        plot.background = element_rect(fill=plot_bg, color=plot_bg),
        panel.background = element_rect(fill=plot_bg),
        legend.key = element_rect(fill=plot_bg),
        text = element_text(color=plot_text_color)
    )
)
dev.off()
embed_fonts(path.expand(filename))


dir.create(paste0("output/plots/local-pareto-exponent"),
    showWarnings=FALSE, recursive=TRUE)
filename <- paste0("output/plots/local-pareto-exponent/alpha2.pdf")
pdf(filename, family=plot_font, width=3, height=3.5)
print(ggplot(df_plot) +
    geom_line(aes(x=p, y=alpha2)) +
    xlab(expression(paste("rank ", italic(p)))) +
    ylab(expression(paste("local Pareto exponent ", alpha[2]))) +
    ggtitle(expression(alpha[2])) +
    theme_bw() + theme(
        legend.justification = c(1, 1),
        legend.position = c(1, 1),
        legend.background = element_rect(
            linetype = "solid",
            color    = plot_text_color,
            size     = 0.3,
            fill     = plot_bg
        ),
        legend.title = element_blank(),
        legend.box.margin = margin(10, 10, 10, 10),
        legend.direction = "horizontal",
        plot.title=element_text(hjust=0.5),
        plot.subtitle=element_text(hjust=0.5),
        plot.background = element_rect(fill=plot_bg, color=plot_bg),
        panel.background = element_rect(fill=plot_bg),
        legend.key = element_rect(fill=plot_bg),
        text = element_text(color=plot_text_color)
    )
)
dev.off()
embed_fonts(path.expand(filename))
