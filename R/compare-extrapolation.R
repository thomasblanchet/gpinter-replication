# ---------------------------------------------------------------------------- #
# “Generalized Pareto Curves: Theory and Applications”, 2017
# Thomas Blanchet, Juliette Fournier, Thomas Piketty
# ---------------------------------------------------------------------------- #
# Compare the different extrapolation methods in tables and graphs.
# ---------------------------------------------------------------------------- #

# Perform the extrapolation with different methods
comparisons <- ddply(dina_data, c("iso", "country", "year", "income_type", "income_type_short"), function(data) {
    average           <- data$average[1]
    income_type_short <- data$income_type_short[1]
    year              <- data$year[1]
    country           <- data$country[1]
    iso               <- data$iso[1]

    # Select the desired test tabulation here
    p_in <- c(0, 0.5, 0.90, 0.99, 1)
    p_out <- 0.999
    #p_in <- c(0, 0.5, 0.8, 0.9, 1)
    #p_out <- 0.99

    # Create the short tabulation (to be used in interpolation)
    data$bracket <- cut(data$p,
        breaks = p_in*1e5,
        include.lowest = TRUE,
        labels = FALSE,
        right = FALSE
    )
    short_tab <- ddply(data, "bracket", function(data) {
        return(data.frame(
            p = min(data$p)/1e5,
            threshold = data$threshold[which.min(data$p)],
            topshare = data$topshare[which.min(data$p)]
        ))
    })

    # Remove the first bracket
    short_tab <- short_tab[-1, ]
    short_tab$m <- average*short_tab$topshare

    # Use the different interpolation methods
    m1 <- method1(p_out, short_tab$p, short_tab$threshold, short_tab$m, average)
    m2 <- method2(p_out, short_tab$p, short_tab$threshold, average)

    # Generalized Pareto interpolation
    dist <- tabulation_fit(short_tab$p, short_tab$threshold, average, topshare=short_tab$topshare)
    m0 <- list(
        threshold = fitted_quantile(dist, p_out),
        topshare = top_share(dist, p_out)
    )

    # Plot the top of the Pareto curve
    p_plot <- seq(0.9, 0.9999, 0.000001)
    df_points <- data.frame(p=data$p/1e5, b=data$invpareto)
    df_points$type <- sapply(data$p, function(p) ifelse(p %in% round(1e5*p_in), "included", "excluded"))
    df_points <- df_points[df_points$p >= 0.9 & df_points$p < 0.9991, ]
    df_line <- data.frame(p=p_plot, b=invpareto(dist, p_plot))
    df_line[df_line$p <= p_in[length(p_in) - 1], "type"] <- "interpolation"
    df_line[df_line$p > p_in[length(p_in) - 1], "type"] <- "extrapolation"


    cat(paste0("Plotting: Pareto curves (extrapolation) - ", country, " - ", income_type_short, " - ", year, "\n"))
    filename <- paste0("output/plots/pareto-curves-extrapolation/pareto-curve-", iso, "-", income_type_short, "-", year, ".pdf")
    pdf(filename, family="CM Roman", width=4.5, height=3.5)
    print(ggplot() +
        geom_line(data=df_line, aes(x=p, y=b, linetype=type)) +
        geom_point(data=df_points, aes(x=p, y=b, shape=type)) +
        scale_linetype_manual("estimation", values=c("interpolation"="solid", "extrapolation"="dashed")) +
        scale_shape_manual("data", values=c("included"=19, "excluded"=1)) +
        xlab(expression(paste("rank ", italic(p)))) +
        ylab(expression(paste("inverted Pareto coefficient ", italic(b), "(", italic(p), ")"))) +
        scale_x_continuous(breaks=seq(0.9, 1, 0.02)) +
        guides(linetype = guide_legend(title.position = "top", order=1),
            shape = guide_legend(title.position = "top", order=2)) +
        theme_bw() + theme(legend.justification = c(0, 1), legend.position = c(0, 1),
            legend.background = element_rect(linetype="solid", color="black", size=0.25),
            legend.box.margin = margin(10, 10, 10, 10),
            legend.direction = "horizontal", plot.title=element_text(hjust=0.5),
            plot.subtitle=element_text(hjust=0.5)) + ggtitle(paste0(country, ", ", year))
    )
    dev.off()
    embed_fonts(path.expand(filename))

    return(data.frame(
        p = p_out,

        threshold_actual = sapply(p_out, function(p) data[data$p == p*1e5, "threshold"])/average,
        threshold_m0 = m0$threshold/average,
        threshold_m1 = m1$threshold/average,
        threshold_m2 = m2$threshold/average,

        topshare_actual = sapply(p_out, function(p) data[data$p == p*1e5, "topshare"]),
        topshare_m0 = m0$topshare,
        topshare_m1 = m1$topshare,
        topshare_m2 = m2$topshare
    ))
})

# Calculate relative error
relerr <- function(a, b) abs(100*(b - a)/a)

extrapolation_re <- data.frame(
    country     = comparisons$country,
    year        = comparisons$year,
    income_type = comparisons$income_type,
    p           = comparisons$p,

    threshold_m0 = relerr(comparisons$threshold_actual, comparisons$threshold_m0),
    threshold_m1 = relerr(comparisons$threshold_actual, comparisons$threshold_m1),
    threshold_m2 = relerr(comparisons$threshold_actual, comparisons$threshold_m2),

    topshare_m0 = relerr(comparisons$topshare_actual, comparisons$topshare_m0),
    topshare_m1 = relerr(comparisons$topshare_actual, comparisons$topshare_m1),
    topshare_m2 = relerr(comparisons$topshare_actual, comparisons$topshare_m2)
)

# Calculate the mean of the error
extrapolation_mre <- ddply(extrapolation_re, c("country", "income_type", "p"), function(data) {
    return(data.frame(
        threshold_m0 = mean(data$threshold_m0),
        threshold_m1 = mean(data$threshold_m1),
        threshold_m2 = mean(data$threshold_m2),

        topshare_m0 = mean(data$topshare_m0),
        topshare_m1 = mean(data$topshare_m1),
        topshare_m2 = mean(data$topshare_m2)
    ))
})



