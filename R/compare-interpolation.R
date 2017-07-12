# ---------------------------------------------------------------------------- #
# “Generalized Pareto Curves: Theory and Applications”, 2017
# Thomas Blanchet, Juliette Fournier, Thomas Piketty
# ---------------------------------------------------------------------------- #
# Compare the different interpolation methods in tables and graphs.
# ---------------------------------------------------------------------------- #

# Perform the interpolation with different methods
comparisons <- ddply(data_micro, .(iso, country, widcode, var_code, var_name, year), function(data) {
    average <- data$average[1]
    var_code <- data$var_code[1]
    var_name <- data$var_name[1]

    if (var_code == "fiinc") {
        p_in <- c(0, 0.4, 0.7, 0.9, 0.99, 1)
        p_out <- c(0.5, 0.8, 0.95)
    } else {
        p_in <- c(0, 0.1, 0.5, 0.9, 0.99, 1)
        p_out <- c(0.30, 0.75, 0.95)
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

    # Use the different interpolation methods
    m1 <- method1(p_out, short_tab$p, short_tab$threshold, short_tab$m, average)
    m2 <- method2(p_out, short_tab$p, short_tab$threshold, average)
    m3 <- method3(p_out, short_tab$p, short_tab$threshold, short_tab$m, average)
    m4 <- method4(p_out, short_tab$p, short_tab$threshold, short_tab$m, average)

    # Generalized Pareto interpolation
    dist <- tabulation_fit(short_tab$p, short_tab$threshold, average, topshare=short_tab$top_share)
    m0 <- list(
        threshold = fitted_quantile(dist, p_out),
        top_share = top_share(dist, p_out)
    )

    return(data.frame(
        p = p_out,

        threshold_actual = sapply(p_out, function(p) data[data$p == round(p*1e5), "threshold"])/average,
        threshold_m0 = m0$threshold/average,
        threshold_m1 = m1$threshold/average,
        threshold_m2 = m2$threshold/average,
        threshold_m3 = m3$threshold/average,
        threshold_m4 = m4$threshold/average,

        top_share_actual = sapply(p_out, function(p) data[data$p == round(p*1e5), "top_share"]),
        top_share_m0 = m0$top_share,
        top_share_m1 = m1$top_share,
        top_share_m2 = m2$top_share,
        top_share_m3 = m3$top_share,
        top_share_m4 = m4$top_share
    ))
})

# Plot the different estimates
d_ply(comparisons, .(iso, country, widcode, var_code, var_name, p), function(data) {
    p <- data$p[1]
    iso <- data$iso[1]
    year <- data$year[1]
    average <- data$average[1]
    country <- data$country[1]
    var_code <- data$var_code[1]
    var_name <- data$var_name[1]

    data$p            <- NULL
    data$threshold_m4 <- NULL
    data$top_share_m4 <- NULL

    data <- melt(data, id.vars=c("iso", "country", "widcode", "var_code", "var_name", "year"))

    cat(paste0("* plotting: comparison time series - ", country, " - ", var_name, " - ", 100*p, "\n"))

    # Create plot for thresholds
    data_thresholds <- data[grepl("threshold", data$variable), ]
    data_thresholds$method[data_thresholds$variable == "threshold_actual"] <- "data"
    data_thresholds$method[data_thresholds$variable == "threshold_m0"] <- "M0"
    data_thresholds$method[data_thresholds$variable == "threshold_m1"] <- "M1"
    data_thresholds$method[data_thresholds$variable == "threshold_m2"] <- "M2"
    data_thresholds$method[data_thresholds$variable == "threshold_m3"] <- "M3"

    data_thresholds$shape[data_thresholds$variable == "threshold_actual"] <- 19
    data_thresholds$shape[data_thresholds$variable == "threshold_m0"] <- 0
    data_thresholds$shape[data_thresholds$variable == "threshold_m1"] <- 2
    data_thresholds$shape[data_thresholds$variable == "threshold_m2"] <- 4
    data_thresholds$shape[data_thresholds$variable == "threshold_m3"] <- 5

    p1 <- ggplot(data_thresholds) +
        geom_line(aes(x=year, y=value, color=method), na.rm=TRUE) +
        geom_point(aes(x=year, y=value, color=method, shape=method), na.rm=TRUE) +
        ylab(paste0("P", 100*p, "/average")) +
        scale_x_continuous(breaks=pretty_breaks()) +
        scale_color_brewer(type="qual", palette="Set1") +
        scale_shape_manual(values=c(19, 0, 2, 4, 5)) +
        theme_bw() + theme(
            legend.title = element_blank(),
            legend.position = "bottom",
            plot.background = element_rect(fill=plot_bg, color=plot_bg),
            panel.background = element_rect(fill=plot_bg),
            legend.key = element_rect(fill=plot_bg),
            legend.background = element_rect(fill=plot_bg),
            text = element_text(color=plot_text_color)
        )

    dir.create("output/plots/comparison-time-series/interpolation", showWarnings=FALSE, recursive=TRUE)
    filename <- paste0("output/plots/comparison-time-series/interpolation/thr-", iso, "-", var_code, "-", 100*p, ".pdf")
    pdf(filename, family=plot_font, width=4.5, height=3.5)
    print(p1)
    dev.off()
    embed_fonts(path.expand(filename))

    # Create plot for top shares
    data_top_share <- data[grepl("top_share", data$variable), ]
    data_top_share$method[data_top_share$variable == "top_share_actual"] <- "data"
    data_top_share$method[data_top_share$variable == "top_share_m0"] <- "M0"
    data_top_share$method[data_top_share$variable == "top_share_m1"] <- "M1"
    data_top_share$method[data_top_share$variable == "top_share_m2"] <- "M2"
    data_top_share$method[data_top_share$variable == "top_share_m3"] <- "M3"

    data_top_share$shape[data_thresholds$variable == "top_share_actual"] <- 19
    data_top_share$shape[data_thresholds$variable == "top_share_m0"] <- 0
    data_top_share$shape[data_thresholds$variable == "top_share_m1"] <- 2
    data_top_share$shape[data_thresholds$variable == "top_share_m2"] <- 4
    data_top_share$shape[data_thresholds$variable == "top_share_m3"] <- 5

    p2 <- ggplot(data_top_share) +
        geom_line(aes(x=year, y=value, color=method), na.rm=TRUE) +
        geom_point(aes(x=year, y=value, color=method, shape=method), na.rm=TRUE) +
        ylab(paste0("top ", 100*(1 - p), "% share")) +
        scale_x_continuous(breaks=pretty_breaks()) +
        scale_y_continuous(labels=percent) +
        scale_color_brewer(type="qual", palette="Set1") +
        scale_shape_manual(values=c(19, 0, 2, 4, 5)) +
        theme_bw() + theme(
            legend.title = element_blank(),
            legend.position = "bottom",
            plot.background = element_rect(fill=plot_bg, color=plot_bg),
            panel.background = element_rect(fill=plot_bg),
            legend.key = element_rect(fill=plot_bg),
            legend.background = element_rect(fill=plot_bg),
            text = element_text(color=plot_text_color)
        )

    dir.create("output/plots/comparison-time-series/interpolation", showWarnings=FALSE, recursive=TRUE)
    filename <- paste0("output/plots/comparison-time-series/interpolation/topshare-", iso, "-", var_code, "-", 100*p, ".pdf")
    pdf(filename, family=plot_font, width=4.5, height=3.5)
    print(p2)
    dev.off()
    embed_fonts(path.expand(filename))

    # Extract the legend
    legend <- g_legend(p1)

    dir.create("output/plots/comparison-time-series/interpolation", showWarnings=FALSE, recursive=TRUE)
    filename <- paste0("output/plots/comparison-time-series/interpolation/", iso, "-", var_code, "-", 100*p, ".pdf")
    pdf(filename, family=plot_font, width=4, height=9)
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

    # Make a plot zooming on the most recent period
    data_thresholds$value[data_thresholds$method == "M1"] <- NA
    data_thresholds$value[data_thresholds$method == "M2"] <- NA
    data_thresholds <- data_thresholds[data_thresholds$year >= 2000, ]

    data_top_share$value[data_top_share$method == "M1"] <- NA
    data_top_share$value[data_top_share$method == "M2"] <- NA
    data_top_share <- data_top_share[data_top_share$year >= 2000, ]

    # Create plot for thresholds
    p1 <- ggplot(data_thresholds) +
        geom_line(aes(x=year, y=value, color=method), na.rm=TRUE) +
        geom_point(aes(x=year, y=value, color=method, shape=method), na.rm=TRUE) +
        ylab(paste0("P", 100*p, "/average")) +
        scale_x_continuous(breaks=pretty_breaks()) +
        scale_color_brewer(type="qual", palette="Set1") +
        scale_shape_manual(values=c(19, 0, 2, 4, 5)) +
        theme_bw() + theme(
            legend.title = element_blank(),
            legend.position = "bottom",
            plot.background = element_rect(fill=plot_bg, color=plot_bg),
            panel.background = element_rect(fill=plot_bg),
            legend.key = element_rect(fill=plot_bg),
            legend.background = element_rect(fill=plot_bg),
            text = element_text(color=plot_text_color)
        )

    dir.create("output/plots/comparison-time-series/interpolation", showWarnings=FALSE, recursive=TRUE)
    filename <- paste0("output/plots/comparison-time-series/interpolation/recent-thr-", iso, "-", var_code, "-", 100*p, ".pdf")
    pdf(filename, family=plot_font, width=4.5, height=3.5)
    print(p1)
    dev.off()
    embed_fonts(path.expand(filename))

    # Create plot for top shares
    p2 <- ggplot(data_top_share) +
        geom_line(aes(x=year, y=value, color=method), na.rm=TRUE) +
        geom_point(aes(x=year, y=value, color=method, shape=method), na.rm=TRUE) +
        ylab(paste0("top ", 100*(1 - p), "% share")) +
        scale_x_continuous(breaks=pretty_breaks()) +
        scale_y_continuous(labels=percent) +
        scale_color_brewer(type="qual", palette="Set1") +
        scale_shape_manual(values=c(19, 0, 2, 4, 5)) +
        theme_bw() + theme(
            legend.title = element_blank(),
            legend.position = "bottom",
            plot.background = element_rect(fill=plot_bg, color=plot_bg),
            panel.background = element_rect(fill=plot_bg),
            legend.key = element_rect(fill=plot_bg),
            legend.background = element_rect(fill=plot_bg),
            text = element_text(color=plot_text_color)
        )

    dir.create("output/plots/comparison-time-series/interpolation", showWarnings=FALSE, recursive=TRUE)
    filename <- paste0("output/plots/comparison-time-series/interpolation/recent-topshare-", iso, "-", var_code, "-", 100*p, ".pdf")
    pdf(filename, family=plot_font, width=4.5, height=3.5)
    print(p2)
    dev.off()
    embed_fonts(path.expand(filename))

    # Extract the legend
    legend <- g_legend(p1)

    dir.create("output/plots/comparison-time-series/interpolation", showWarnings=FALSE, recursive=TRUE)
    filename <- paste0("output/plots/comparison-time-series/interpolation/recent-", iso, "-", var_code, "-", 100*p, ".pdf")
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

interpolation_re <- data.frame(
    p        = comparisons$p,
    iso      = comparisons$iso,
    year     = comparisons$year,
    country  = comparisons$country,
    var_code = comparisons$var_code,
    var_name = comparisons$var_name,

    threshold_m0 = relerr(comparisons$threshold_actual, comparisons$threshold_m0),
    threshold_m1 = relerr(comparisons$threshold_actual, comparisons$threshold_m1),
    threshold_m2 = relerr(comparisons$threshold_actual, comparisons$threshold_m2),
    threshold_m3 = relerr(comparisons$threshold_actual, comparisons$threshold_m3),
    threshold_m4 = relerr(comparisons$threshold_actual, comparisons$threshold_m4),

    top_share_m0 = relerr(comparisons$top_share_actual, comparisons$top_share_m0),
    top_share_m1 = relerr(comparisons$top_share_actual, comparisons$top_share_m1),
    top_share_m2 = relerr(comparisons$top_share_actual, comparisons$top_share_m2),
    top_share_m3 = relerr(comparisons$top_share_actual, comparisons$top_share_m3),
    top_share_m4 = relerr(comparisons$top_share_actual, comparisons$top_share_m4)
)

# Calculate the mean of the error
interpolation_mre <- ddply(interpolation_re, .(country, iso, var_code, var_name, p), function(data) {
    return(data.frame(
        threshold_m0 = mean(data$threshold_m0),
        threshold_m1 = mean(data$threshold_m1),
        threshold_m2 = mean(data$threshold_m2),
        threshold_m3 = mean(data$threshold_m3),
        threshold_m4 = mean(data$threshold_m4),

        top_share_m0 = mean(data$top_share_m0),
        top_share_m1 = mean(data$top_share_m1),
        top_share_m2 = mean(data$top_share_m2),
        top_share_m3 = mean(data$top_share_m3),
        top_share_m4 = mean(data$top_share_m4)
    ))
})

# Generate LaTeX tables comparing the results
d_ply(interpolation_mre, "var_code", function(data) {
    dir.create("output/tables/compare-interpolation", showWarnings=FALSE, recursive=TRUE)
    filename <- paste0("output/tables/compare-interpolation/compare-", data$var_code[1], ".tex")

    # Remove the M4 method from the comparison (comment to keep)
    data$threshold_m4 <- NULL
    data$top_share_m4 <- NULL
    # Number of methods left to compare
    n <- (ncol(data) - 5)/2
    # Column for thresholds
    threshold_cols <- paste0("threshold_m", 0:(n - 1))
    topshare_cols <- paste0("top_share_m", 0:(n - 1))

    sink(filename)

    cat(paste0(c("\\begin{tabular}{cc", rep("P@{}", n), "} \\toprule\n"), collapse=""))
    cat(paste0("& & \\multicolumn{", n, "}{>{\\centering\\arraybackslash}p{", 2*n, "cm}}"))
    cat("{mean percentage gap between estimated and observed values}")
    cat(paste0("\\\\ \\cmidrule(l){3-", 3 + n - 1, "}\n"))
    cat("& & ")
    cat(paste0(paste0("M", 0:(n - 1), collapse=" & "), " \\\\ \\midrule\n"))

    # United States
    data_us <- subset(data, iso == "US")
    for (i in 1:3) {
        if (i == 1) {
            cat("\\multirow{14}{3cm}{\\centering United States \\\\ (1962--2014)} & ")
        } else {
            cat("& ")
        }
        cat(sprintf("\\multirow{2}{*}{Top %.0f\\%% share} & ", 100*(1 - data_us$p[i])))
        cat(paste0(format(data_us[i, topshare_cols], digits=2, scientific=FALSE), collapse="\\% & "))
        cat("\\% \\\\ \n")
        cat("& & \\footnotesize (ref.)")
        for (j in 1:(n - 1)) {
            cat(" & \\footnotesize ($\\times ")
            cat(format(data_us[i, paste0("top_share_m", j)]/data_us[i, "top_share_m0"], digits=2, scientific=FALSE))
            cat("$)")
        }
        cat(paste0("\\\\ \\cmidrule(l){2-", 2 + n, "}\n"))
    }
    for (i in 1:3) {
        cat("& ")
        cat(sprintf("\\multirow{2}{*}{P%.0f/average} & ", 100*data_us$p[i]))
        cat(paste0(format(data_us[i, threshold_cols], digits=2, scientific=FALSE), collapse="\\% & "))
        cat("\\% \\\\ \n")
        cat("& & \\footnotesize (ref.)")
        for (j in 1:(n - 1)) {
            cat(" & \\footnotesize ($\\times ")
            cat(format(data_us[i, paste0("threshold_m", j)]/data_us[i, "threshold_m0"], digits=2, scientific=FALSE))
            cat("$)")
        }
        if (i == 3) {
            cat(paste0("\\\\ \\midrule\n"))
        } else {
            cat(paste0("\\\\ \\cmidrule(l){2-", 2 + n, "}\n"))
        }
    }

    # France
    data_fr <- subset(data, iso == "FR")
    for (i in 1:3) {
        if (i == 1) {
            cat("\\multirow{14}{3cm}{\\centering France \\\\ (1994--2012)} & ")
        } else {
            cat("& ")
        }
        cat(sprintf("\\multirow{2}{*}{Top %.0f\\%% share} & ", 100*(1 - data_fr$p[i])))
        cat(paste0(format(data_fr[i, topshare_cols], digits=2, scientific=FALSE), collapse="\\% & "))
        cat("\\% \\\\ \n")
        cat("& & \\footnotesize (ref.)")
        for (j in 1:(n - 1)) {
            cat(" & \\footnotesize ($\\times ")
            cat(format(data_fr[i, paste0("top_share_m", j)]/data_fr[i, "top_share_m0"], digits=2, scientific=FALSE))
            cat("$)")
        }
        cat(paste0("\\\\ \\cmidrule(l){2-", 2 + n, "}\n"))
    }
    for (i in 1:3) {
        cat("& ")
        cat(sprintf("\\multirow{2}{*}{P%.0f/average} & ", 100*data_fr$p[i]))
        cat(paste0(format(data_fr[i, threshold_cols], digits=2, scientific=FALSE), collapse="\\% & "))
        cat("\\% \\\\ \n")
        cat("& & \\footnotesize (ref.)")
        for (j in 1:(n - 1)) {
            cat(" & \\footnotesize ($\\times ")
            cat(format(data_fr[i, paste0("threshold_m", j)]/data_fr[i, "threshold_m0"], digits=2, scientific=FALSE))
            cat("$)")
        }
        if (i == 3) {
            cat(paste0("\\\\ \\bottomrule\n"))
        } else {
            cat(paste0("\\\\ \\cmidrule(l){2-", 2 + n, "}\n"))
        }
    }

    cat("\\end{tabular}\n")

    sink()
})

