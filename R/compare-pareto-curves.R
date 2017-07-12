# ---------------------------------------------------------------------------- #
# “Generalized Pareto Curves: Theory and Applications”, 2017
# Thomas Blanchet, Juliette Fournier, Thomas Piketty
# ---------------------------------------------------------------------------- #
# Shape of generalized Pareto curves with different interpolation methods.
# ---------------------------------------------------------------------------- #

d_ply(data_us_fr, .(iso, country, widcode, var_code, var_name, year), function(data) {
    iso <- data$iso[1]
    year <- data$year[1]
    average <- data$average[1]
    country <- data$country[1]
    var_code <- data$var_code[1]
    var_name <- data$var_name[1]

    if (var_code %in% c("fiinc", "hweal", "pkkin", "pllin")) {
        p_in <- c(0, 0.4, 0.7, 0.9, 0.99, 1)
    } else if (var_code == "ptinc") {
        p_in <- c(0, 0.1, 0.5, 0.9, 0.99, 1)
    }

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
            top_share = data$top_share[which.min(data$p)]
        ))
    })
    # Remove the first bracket, which includes zero or negative values
    short_tab <- short_tab[-1, ]
    short_tab$m <- average*short_tab$top_share

    # Tight grid in output
    p_out <- seq(0, 0.999, 0.0001)

    # Use the different interpolation methods
    m1 <- method1(p_out, short_tab$p, short_tab$threshold, short_tab$m, average)
    m2 <- method2(p_out, short_tab$p, short_tab$threshold, average)
    m3 <- method3(p_out, short_tab$p, short_tab$threshold, short_tab$m, average)

    m1$b <- average*m1$top_share/(1 - p_out)/m1$threshold
    m2$b <- average*m2$top_share/(1 - p_out)/m2$threshold
    m3$b <- average*m3$top_share/(1 - p_out)/m3$threshold

    # Generalized Pareto interpolation
    dist <- tabulation_fit(short_tab$p, short_tab$threshold, average, topshare=short_tab$top_share)
    m0 <- list(b=invpareto(dist, p_out))

    # Data frame with the curve for the different interpolation methods
    df <- data.frame(
        p = p_out,
        m0 = m0$b,
        m1 = m1$b,
        m2 = m2$b,
        m3 = m3$b
    )
    df <- melt(df, id.vars="p")
    df$variable <- toupper(df$variable)
    df <- df[df$p >= 0.7, ]

    knots <- data.frame(
        p=p_in[2:5], b=unlist(sapply(p_in[2:5], function(p) data[data$p == round(p*1e5), "b"]))
    )
    knots <- knots[knots$p >= 0.7, ]

    actual <- data.frame(
        p = data$p/1e5,
        b = data$b
    )
    actual <- actual[actual$p >= 0.7 & actual$p <= 0.999, ]

    cat(paste0("* plotting: comparison Pareto curves - ", country, " - ", var_name, " - ", year, "\n"))

    dir.create("output/plots/comparison-pareto-curves", showWarnings=FALSE, recursive=TRUE)
    filename <- paste0("output/plots/comparison-pareto-curves/", iso, "-", var_code, "-", year, ".pdf")
    pdf(filename, family=plot_font, width=6, height=3.5)
    print(ggplot(df) +
        geom_line(aes(x=p, y=value, linetype="estimated", color="estimated"), na.rm=TRUE) +
        geom_line(data=actual, aes(x=p, y=b, linetype="data", color="data")) +
        geom_point(data=knots, aes(x=p, y=b)) +
        facet_wrap("variable", nrow=2, ncol=2) +
        scale_color_brewer(name="source", type="qual", palette="Set1") +
        scale_linetype_manual(name="source", values=c("data"="dashed", "estimated"="solid")) +
        xlab(expression(paste("rank ", italic(p)))) +
        ylab(expression(paste("inverted Pareto coefficient ", italic(b), "(", italic(p), ")"))) +
        theme_bw() + theme(
            legend.position = "right",
            legend.title = element_blank(),
            plot.background = element_rect(fill=plot_bg, color=plot_bg),
            panel.background = element_rect(fill=plot_bg),
            legend.background = element_rect(fill=plot_bg),
            legend.key = element_rect(fill=plot_bg),
            text = element_text(color=plot_text_color)
        )
    )
    dev.off()
    embed_fonts(path.expand(filename))
})

