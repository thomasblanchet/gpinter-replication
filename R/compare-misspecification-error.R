# ---------------------------------------------------------------------------- #
# “Generalized Pareto Curves: Theory and Applications”, 2017
# Thomas Blanchet, Juliette Fournier, Thomas Piketty
# ---------------------------------------------------------------------------- #
# Compare the misspecification error with the actual error for the US
# distribution of labor income in 1962.
# ---------------------------------------------------------------------------- #

data_us_labor <- subset(dina_data_us, income_type_short == "labor" & year == 1962)

# Interpolate the distribution
pk <- c(0.1, 0.5, 0.9, 0.99)
xk <- -log(1 - pk)
average <- data_us_labor$average[1]
data_us_labor$bracket <- cut(data_us_labor$p,
    breaks = 1e5*c(0, pk, 1),
    include.lowest = TRUE,
    labels = FALSE,
    right = FALSE
)
short_tab <- ddply(data_us_labor, "bracket", function(data) {
    return(data.frame(
        p = min(data$p)/1e5,
        threshold = data$threshold[which.min(data$p)],
        topshare = data$topshare[which.min(data$p)]
    ))
})
short_tab <- short_tab[-1, ]
short_tab$m <- average*short_tab$topshare
dist <- tabulation_fit(short_tab$p, short_tab$threshold, average, topshare=short_tab$topshare)

# Calculate the function phi'''
data_us_labor <- subset(data_us_labor, p >= 10000 & p <= 99000)
x_out <- seq(min(xk), max(xk), length.out=1000)
data_us_labor$x <- -log(1 - data_us_labor$p/1e5)
phi_d3_data <- sapply(x_out, function(x) {
    return(lpoly(x, data_us_labor$x, data_us_labor$phi, data_us_labor$dphi, 0.05))
})
phi_d3_fun <- splinefun(x_out, phi_d3_data, method="monoH.FC")

# Calculate the error on phi and phi', using both methods
x_out2 <- seq(min(xk), max(xk), length.out=200)
df <- rbind(
    data.frame(
        x = data_us_labor$x,
        err_phi = abs(phi(dist, data_us_labor$x) - data_us_labor$phi),
        err_dphi = abs(deriv_phi(dist, data_us_labor$x) - data_us_labor$dphi),
        type = rep("observed", nrow(data_us_labor))
    ),
    data.frame(
        x = x_out2,
        err_phi = abs(interpolation_value_error_noncons(x_out2, xk, phi_d3_fun)),
        err_dphi = abs(interpolation_deriv_error_noncons(x_out2, xk, phi_d3_fun)),
        type = rep("estimated", length(x_out2))
    )
)

df$err_dphi[df$type == "observed"] <- ksmooth(
    x = df$x[df$type == "observed"],
    y = df$err_dphi[df$type == "observed"],
    x.points = df$x[df$type == "observed"],
    bandwidth = 0.25
)$y

filename <- "output/plots/compare-misspecification-error/compare-error-phi.pdf"
pdf(filename, family="CM Roman", width=4.5, height=3.5)
print(ggplot(df) +
    geom_line(aes(x=x, y=err_phi, color=type, linetype=type)) +
    geom_vline(xintercept=xk, linetype="longdash") +
    scale_color_brewer("type", type="qual", palette="Set1") +
    scale_linetype_manual("type", values=c("observed"="longdash", "estimated"="solid")) +
    xlab(expression(paste(italic(x)==-log, "(", 1-italic(p), ")"))) +
    ylab("absolute error") +
    ggtitle(expression(paste("Error on ", phi1, "(", italic(x), ")")),
        subtitle=expression("US labor income, 1962")) +
    theme_bw() + theme(
        legend.justification = c(1, 1),
        legend.position = c(1, 1),
        legend.background = element_rect(linetype="solid", color="black", size=0.3),
        legend.title = element_blank(),
        legend.box.margin = margin(10, 10, 10, 10),
        legend.direction = "horizontal",
        plot.title=element_text(hjust=0.5),
        plot.subtitle=element_text(hjust=0.5)
    )
)
dev.off()
embed_fonts(path.expand(filename))

filename <- "output/plots/compare-misspecification-error/compare-error-deriv-phi.pdf"
pdf(filename, family="CM Roman", width=4.5, height=3.5)
print(ggplot(df) +
    geom_line(aes(x=x, y=err_dphi, color=type, linetype=type)) +
    geom_vline(xintercept=xk, linetype="longdash") +
    scale_color_brewer("type", type="qual", palette="Set1") +
    scale_linetype_manual("type", values=c("observed"="longdash", "estimated"="solid")) +
    xlab(expression(paste(italic(x)==-log, "(", 1-italic(p), ")"))) +
    ylab("absolute error") +
    ggtitle(expression(paste("Error on ", phi1, "'(", italic(x), ")")),
        subtitle=expression("US labor income, 1962")) +
    theme_bw() + theme(
        legend.justification = c(1, 1),
        legend.position = c(1, 1),
        legend.background = element_rect(linetype="solid", color="black", size=0.3),
        legend.title = element_blank(),
        legend.box.margin = margin(10, 10, 10, 10),
        legend.direction = "horizontal",
        plot.title=element_text(hjust=0.5),
        plot.subtitle=element_text(hjust=0.5)
    )
)
dev.off()
embed_fonts(path.expand(filename))


