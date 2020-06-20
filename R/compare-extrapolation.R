# ---------------------------------------------------------------------------- #
# “Generalized Pareto Curves: Theory and Applications”, 2017
# Thomas Blanchet, Juliette Fournier, Thomas Piketty
# ---------------------------------------------------------------------------- #
# Compare the different extrapolation methods in tables and graphs.
# ---------------------------------------------------------------------------- #

# Perform the extrapolation with different methods
comparisons <- ddply(data_micro, .(iso, country, widcode, var_code, var_name, year), function(data) {
    iso <- data$iso[1]
    year <- data$year[1]
    average <- data$average[1]
    country <- data$country[1]
    var_code <- data$var_code[1]
    var_name <- data$var_name[1]

    # Test tabulation here
    #p_in <- c(0, 0.5, 0.9, 0.95, 1)
    #p_out <- 0.99
    p_in <- c(0, 0.5, 0.9, 0.99, 1)
    p_out <- 0.999

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

    # Remove the first bracket
    short_tab <- short_tab[-1, ]
    short_tab$m <- average*short_tab$top_share

    # Use the different interpolation methods
    m1 <- method1(p_out, short_tab$p, short_tab$threshold, short_tab$m, average)
    m2 <- method2(p_out, short_tab$p, short_tab$threshold, average)
    m5 <- method5(p_out, short_tab$p, short_tab$m, average)
    m6 <- method6(p_out, short_tab$p, short_tab$m, average)

    # Generalized Pareto interpolation
    dist <- tabulation_fit(short_tab$p, short_tab$threshold, average, topshare=short_tab$top_share)
    m0 <- list(
        threshold = fitted_quantile(dist, p_out),
        top_share = top_share(dist, p_out)
    )

    # Plot the top of the Pareto curve
    # p_plot <- seq(0.9, 0.9999, 0.000001)
    # df_points <- data.frame(p=data$p/1e5, b=data$b)
    # df_points$type <- sapply(data$p, function(p) ifelse(p %in% round(1e5*p_in), "included", "excluded"))
    # df_points <- df_points[df_points$p >= 0.9 & df_points$p < 0.9991, ]
    # df_line <- data.frame(p=p_plot, b=invpareto(dist, p_plot))
    # df_line[df_line$p <= p_in[length(p_in) - 1], "type"] <- "interpolation"
    # df_line[df_line$p > p_in[length(p_in) - 1], "type"] <- "extrapolation"
    #
    # cat(paste0("* plotting: Pareto curves (extrapolation) - ", country, " - ", var_name, " - ", year, "\n"))
    # dir.create("output/plots/pareto-curves-extrapolation", showWarnings=FALSE, recursive=TRUE)
    # filename <- paste0("output/plots/pareto-curves-extrapolation/", iso, "-", var_code, "-", year, ".pdf")
    # pdf(filename, family=plot_font, width=4.5, height=3.5)
    # print(ggplot() +
    #     geom_line(data=df_line, aes(x=p, y=b, linetype=type)) +
    #     geom_point(data=df_points, aes(x=p, y=b, shape=type)) +
    #     scale_linetype_manual("estimation", values=c("interpolation"="solid", "extrapolation"="dashed")) +
    #     scale_shape_manual("data", values=c("included"=19, "excluded"=1)) +
    #     xlab(expression(paste("rank ", italic(p)))) +
    #     ylab(expression(paste("inverted Pareto coefficient ", italic(b), "(", italic(p), ")"))) +
    #     scale_x_continuous(breaks=seq(0.9, 1, 0.02)) +
    #     guides(linetype = guide_legend(title.position = "top", order=1),
    #         shape = guide_legend(title.position = "top", order=2)) +
    #     theme_bw() + theme(legend.justification = c(0, 1), legend.position = c(0, 1),
    #         legend.background = element_rect(linetype="solid", color="black", size=0.25, fill=plot_bg),
    #         legend.box.margin = margin(10, 10, 10, 10),
    #         legend.direction = "horizontal",
    #         plot.title = element_text(hjust=0.5),
    #         plot.subtitle = element_text(hjust=0.5),
    #         plot.background = element_rect(fill=plot_bg, color=plot_bg),
    #         panel.background = element_rect(fill=plot_bg),
    #         legend.key = element_rect(fill=plot_bg),
    #         text = element_text(color=plot_text_color)
    #     ) + ggtitle(paste0(country, ", ", year))
    # )
    # dev.off()
    # embed_fonts(path.expand(filename))

    return(data.frame(
        p = p_out,

        threshold_actual = sapply(p_out, function(p) data[data$p == p*1e5, "threshold"])/average,
        threshold_m0 = m0$threshold/average,
        threshold_m1 = m1$threshold/average,
        threshold_m2 = m2$threshold/average,
        threshold_m5 = m5$threshold/average,
        threshold_m6 = m6$threshold/average,

        top_share_actual = sapply(p_out, function(p) data[data$p == p*1e5, "top_share"]),
        top_share_m0 = m0$top_share,
        top_share_m1 = m1$top_share,
        top_share_m2 = m2$top_share,
        top_share_m5 = m5$top_share,
        top_share_m6 = m6$top_share
    ))
})

# Plot the time series
d_ply(comparisons, .(iso, country, widcode, var_code, var_name, p), function(data) {
    p <- data$p[1]
    iso <- data$iso[1]
    year <- data$year[1]
    average <- data$average[1]
    country <- data$country[1]
    var_code <- data$var_code[1]
    var_name <- data$var_name[1]

    data$p <- NULL
    data <- melt(data, id.vars=c("iso", "country", "widcode", "var_code", "var_name", "year"))

    cat(paste0("* plotting: comparison time series - ", country, " - ", var_name, " - ", 100*p, "\n"))

    # Create plot for thresholds
    data_thresholds <- data[grepl("threshold", data$variable), ]
    data_thresholds$method[data_thresholds$variable == "threshold_actual"] <- "data"
    data_thresholds$method[data_thresholds$variable == "threshold_m0"] <- "M0"
    data_thresholds$method[data_thresholds$variable == "threshold_m1"] <- "M1"
    data_thresholds$method[data_thresholds$variable == "threshold_m2"] <- "M2"

    p1 <- ggplot(data_thresholds) +
        geom_line(aes(x=year, y=value, color=method), na.rm=TRUE) +
        geom_point(aes(x=year, y=value, color=method, shape=method), na.rm=TRUE) +
        ylab(paste0("P", 100*p, "/average")) +
        scale_color_brewer(type="qual", palette="Set1") +
        scale_shape_manual(values=c(19, 0, 2, 4)) +
        theme_bw() + theme(
            legend.title = element_blank(),
            legend.position = "bottom",
            plot.background = element_rect(fill=plot_bg, color=plot_bg),
            panel.background = element_rect(fill=plot_bg),
            legend.key = element_rect(fill=plot_bg),
            legend.background = element_rect(fill=plot_bg),
            text = element_text(color=plot_text_color)
        )

    # Create plot for top shares
    data_topshare <- data[grepl("top_share", data$variable), ]
    data_topshare$method[data_topshare$variable == "top_share_actual"] <- "data"
    data_topshare$method[data_topshare$variable == "top_share_m0"] <- "M0"
    data_topshare$method[data_topshare$variable == "top_share_m1"] <- "M1"
    data_topshare$method[data_topshare$variable == "top_share_m2"] <- "M2"

    p2 <- ggplot(data_topshare) +
        geom_line(aes(x=year, y=value, color=method), na.rm=TRUE) +
        geom_point(aes(x=year, y=value, color=method, shape=method), na.rm=TRUE) +
        ylab(paste0("top ", 100*(1 - p), "% share")) +
        scale_y_continuous(labels=percent) +
        scale_color_brewer(type="qual", palette="Set1") +
        scale_shape_manual(values=c(19, 0, 2, 4)) +
        theme_bw() + theme(
            legend.position = "none",
            plot.background = element_rect(fill=plot_bg, color=plot_bg),
            panel.background = element_rect(fill=plot_bg),
            legend.key = element_rect(fill=plot_bg),
            legend.background = element_rect(fill=plot_bg),
            text = element_text(color=plot_text_color)
        )

    # Extract the legend
    legend <- g_legend(p1)

    dir.create("output/plots/comparison-time-series/extrapolation", showWarnings=FALSE, recursive=TRUE)
    filename <- paste0("output/plots/comparison-time-series/extrapolation/", iso, "-", var_code, "-", 100*p, ".pdf")
    pdf(filename, family=plot_font, width=9, height=4)
    grid.arrange(
        arrangeGrob(
            p1 + theme(legend.position="none"),
            p2 + theme(legend.position="none"),
            nrow = 1
        ),
        legend, nrow=2, heights=c(10, 1)
    )
    dev.off()
    embed_fonts(path.expand(filename))
})

# Calculate relative error
relerr <- function(a, b) abs(100*(b - a)/a)

extrapolation_re <- data.frame(
    p        = comparisons$p,
    iso      = comparisons$iso,
    year     = comparisons$year,
    country  = comparisons$country,
    var_code = comparisons$var_code,
    var_name = comparisons$var_name,

    threshold_m0 = relerr(comparisons$threshold_actual, comparisons$threshold_m0),
    threshold_m1 = relerr(comparisons$threshold_actual, comparisons$threshold_m1),
    threshold_m2 = relerr(comparisons$threshold_actual, comparisons$threshold_m2),
    threshold_m5 = relerr(comparisons$threshold_actual, comparisons$threshold_m5),
    threshold_m6 = relerr(comparisons$threshold_actual, comparisons$threshold_m6),

    top_share_m0 = relerr(comparisons$top_share_actual, comparisons$top_share_m0),
    top_share_m1 = relerr(comparisons$top_share_actual, comparisons$top_share_m1),
    top_share_m2 = relerr(comparisons$top_share_actual, comparisons$top_share_m2),
    top_share_m5 = relerr(comparisons$top_share_actual, comparisons$top_share_m5),
    top_share_m6 = relerr(comparisons$top_share_actual, comparisons$top_share_m6)
)

# Calculate the mean of the error
extrapolation_mre <- ddply(extrapolation_re, .(country, iso, var_code, var_name, p), function(data) {
    return(data.frame(
        threshold_m0 = mean(data$threshold_m0),
        threshold_m1 = mean(data$threshold_m1),
        threshold_m2 = mean(data$threshold_m2),
        threshold_m5 = mean(data$threshold_m5),
        threshold_m6 = mean(data$threshold_m6),

        top_share_m0 = mean(data$top_share_m0),
        top_share_m1 = mean(data$top_share_m1),
        top_share_m2 = mean(data$top_share_m2),
        top_share_m5 = mean(data$top_share_m5),
        top_share_m6 = mean(data$top_share_m6)
    ))
})

# Generate LaTeX tables comparing the results
d_ply(extrapolation_mre, "var_code", function(data) {
    dir.create("output/tables/compare-extrapolation", showWarnings=FALSE, recursive=TRUE)
    filename <- paste0("output/tables/compare-extrapolation/compare-", data$var_code[1], ".tex")

    # Number of methods to compare
    cols <- c(0, 1, 2, 5, 6)
    n <- length(cols)
    # Column for thresholds
    threshold_cols <- paste0("threshold_m", cols)
    top_share_cols <- paste0("top_share_m", cols)

    sink(filename)

    cat(paste0(c("\\begin{tabular}{cc", rep("P@{}", n), "} \\toprule\n"), collapse=""))
    cat(paste0("& & \\multicolumn{", n, "}{>{\\centering\\arraybackslash}p{", 1.5*n, "cm}}"))
    cat("{mean percentage gap between estimated and observed values}")
    cat(paste0("\\\\ \\cmidrule(l){3-", 3 + n - 1, "}\n"))
    cat("& & ")
    cat(paste0(paste0("M", cols, collapse=" & "), " \\\\ \\midrule\n"))

    # United States
    data_us <- subset(data, iso == "US")
    cat("\\multirow{4}{3cm}{\\centering United States \\\\ (1962--2014)} & ")
    cat(paste0("\\multirow{2}{*}{Top ", 100*(1 - data_us$p), "\\% share} & "))
    cat(paste0(format(data_us[1, top_share_cols], digits=2, scientific=FALSE), collapse="\\% & "))
    cat("\\% \\\\ \n")
    cat("& & \\footnotesize (ref.)")
    for (j in cols[-1]) {
        cat(" & \\footnotesize ($\\times ")
        cat(format(data_us[1, paste0("top_share_m", j)]/data_us[1, "top_share_m0"], digits=2, scientific=FALSE))
        cat("$)")
    }
    cat(paste0("\\\\ \\cmidrule(l){2-", 2 + n, "}\n"))

    cat("& ")
    cat(paste0("\\multirow{2}{*}{P", 100*data_us$p,"/average} & "))
    cat(paste0(format(data_us[1, threshold_cols], digits=2, scientific=FALSE), collapse="\\% & "))
    cat("\\% \\\\ \n")
    cat("& & \\footnotesize (ref.)")
    for (j in cols[-1]) {
        cat(" & \\footnotesize ($\\times ")
        cat(format(data_us[1, paste0("threshold_m", j)]/data_us[1, "threshold_m0"], digits=2, scientific=FALSE))
        cat("$)")
    }
    cat(paste0("\\\\ \\midrule\n"))

    # France
    data_fr <- subset(data, iso == "FR")
    cat("\\multirow{4}{3cm}{\\centering France \\\\ (1994--2012)} & ")
    cat(paste0("\\multirow{2}{*}{Top ", 100*(1 - data_fr$p), "\\% share} & "))
    cat(paste0(format(data_fr[1, top_share_cols], digits=2, scientific=FALSE), collapse="\\% & "))
    cat("\\% \\\\ \n")
    cat("& & \\footnotesize (ref.)")
    for (j in cols[-1]) {
        cat(" & \\footnotesize ($\\times ")
        cat(format(data_fr[1, paste0("top_share_m", j)]/data_fr[1, "top_share_m0"], digits=2, scientific=FALSE))
        cat("$)")
    }
    cat(paste0("\\\\ \\cmidrule(l){2-", 2 + n, "}\n"))

    cat("& ")
    cat(paste0("\\multirow{2}{*}{P", 100*data_fr$p,"/average} & "))
    cat(paste0(format(data_fr[1, threshold_cols], digits=2, scientific=FALSE), collapse="\\% & "))
    cat("\\% \\\\ \n")
    cat("& & \\footnotesize (ref.)")
    for (j in cols[-1]) {
        cat(" & \\footnotesize ($\\times ")
        cat(format(data_fr[1, paste0("threshold_m", j)]/data_fr[1, "threshold_m0"], digits=2, scientific=FALSE))
        cat("$)")
    }
    cat(paste0("\\\\ \\bottomrule\n"))
    cat("\\end{tabular}\n")

    sink()
})



