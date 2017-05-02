# ---------------------------------------------------------------------------- #
# “Generalized Pareto Curves: Theory and Applications”, 2017
# Thomas Blanchet, Juliette Fournier, Thomas Piketty
# ---------------------------------------------------------------------------- #
# Compare the different interpolation methods in tables and graphs.
# ---------------------------------------------------------------------------- #

# Perform the interpolation with different methods
comparisons <- ddply(dina_data, c("iso", "country", "year", "income_type", "income_type_short"), function(data) {
    average <- data$average[1]
    income_type_short <- data$income_type_short[1]

    if (income_type_short == "fiscal") {
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
            topshare = data$topshare[which.min(data$p)]
        ))
    })
    # Remove the first bracket, which includes zero or negative values
    short_tab <- short_tab[-1, ]
    short_tab$m <- average*short_tab$topshare

    # Use the different interpolation methods
    m1 <- method1(p_out, short_tab$p, short_tab$threshold, short_tab$m, average)
    m2 <- method2(p_out, short_tab$p, short_tab$threshold, average)
    m3 <- method3(p_out, short_tab$p, short_tab$threshold, short_tab$m, average)
    m4 <- method4(p_out, short_tab$p, short_tab$threshold, short_tab$m, average)

    # Generalized Pareto interpolation
    dist <- tabulation_fit(short_tab$p, short_tab$threshold, average, topshare=short_tab$topshare)
    m0 <- list(
        threshold = fitted_quantile(dist, p_out),
        topshare = top_share(dist, p_out)
    )

    return(data.frame(
        p = p_out,

        threshold_actual = sapply(p_out, function(p) data[data$p == p*1e5, "threshold"])/average,
        threshold_m0 = m0$threshold/average,
        threshold_m1 = m1$threshold/average,
        threshold_m2 = m2$threshold/average,
        threshold_m3 = m3$threshold/average,
        threshold_m4 = m4$threshold/average,

        topshare_actual = sapply(p_out, function(p) data[data$p == p*1e5, "topshare"]),
        topshare_m0 = m0$topshare,
        topshare_m1 = m1$topshare,
        topshare_m2 = m2$topshare,
        topshare_m3 = m3$topshare,
        topshare_m4 = m4$topshare
    ))
})

# Plot the different estimates
d_ply(comparisons, c("country", "income_type", "income_type_short", "p"), function(data) {
    country           <- data$country[1]
    iso               <- data$iso[1]
    income_type_short <- data$income_type_short[1]
    p                 <- data$p[1]

    data$p            <- NULL
    data$threshold_m4 <- NULL
    data$topshare_m4  <- NULL

    data <- melt(data, id.vars=c("iso", "country", "year", "income_type", "income_type_short"))

    cat(paste0("Plotting: comparison time series - ", country, " - ", income_type_short, " - ", 100*p, "\n"))

    # Create plot for thresholds
    data_thresholds <- data[grepl("threshold", data$variable), ]
    data_thresholds$method[data_thresholds$variable == "threshold_actual"] <- "data"
    data_thresholds$method[data_thresholds$variable == "threshold_m0"] <- "M0"
    data_thresholds$method[data_thresholds$variable == "threshold_m1"] <- "M1"
    data_thresholds$method[data_thresholds$variable == "threshold_m2"] <- "M2"
    data_thresholds$method[data_thresholds$variable == "threshold_m3"] <- "M3"

    p1 <- ggplot(data_thresholds) +
        geom_line(aes(x=year, y=value, color=method), na.rm=TRUE) +
        geom_point(aes(x=year, y=value, color=method, shape=method), na.rm=TRUE) +
        ylab(paste0("P", 100*p, "/average")) +
        scale_color_brewer(type="qual", palette="Set1") +
        theme_bw() + theme(
            legend.title = element_blank(),
            legend.position = "bottom",
            plot.background = element_rect(fill=plot_bg, color=plot_bg),
            panel.background = element_rect(fill=plot_bg),
            legend.key = element_rect(fill=plot_bg),
            text = element_text(color=plot_text_color)
        )

    # Create plot for top shares
    data_topshare <- data[grepl("topshare", data$variable), ]
    data_topshare$method[data_topshare$variable == "topshare_actual"] <- "data"
    data_topshare$method[data_topshare$variable == "topshare_m0"] <- "M0"
    data_topshare$method[data_topshare$variable == "topshare_m1"] <- "M1"
    data_topshare$method[data_topshare$variable == "topshare_m2"] <- "M2"
    data_topshare$method[data_topshare$variable == "topshare_m3"] <- "M3"

    p2 <- ggplot(data_topshare) +
        geom_line(aes(x=year, y=value, color=method), na.rm=TRUE) +
        geom_point(aes(x=year, y=value, color=method, shape=method), na.rm=TRUE) +
        ylab(paste0("top ", 100*(1 - p), "% share")) +
        scale_y_continuous(labels=percent) +
        scale_color_brewer(type="qual", palette="Set1") +
        theme_bw() + theme(
            legend.position = "none",
            plot.background = element_rect(fill=plot_bg, color=plot_bg),
            panel.background = element_rect(fill=plot_bg),
            legend.key = element_rect(fill=plot_bg),
            text = element_text(color=plot_text_color)
        )

    # Extract the legend
    legend <- g_legend(p1)

    filename <- paste0("output/plots/comparison-time-series/time-series-", iso, "-", income_type_short, "-", 100*p, ".pdf")
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

    # Make a plot zooming on the most recent period
    data_thresholds$value[data_thresholds$method == "M1"] <- NA
    data_thresholds$value[data_thresholds$method == "M2"] <- NA
    data_thresholds <- data_thresholds[data_thresholds$year >= 2000, ]

    data_topshare$value[data_topshare$method == "M1"] <- NA
    data_topshare$value[data_topshare$method == "M2"] <- NA
    data_topshare <- data_topshare[data_topshare$year >= 2000, ]

    # Create plot for thresholds
    p1 <- ggplot(data_thresholds) +
        geom_line(aes(x=year, y=value, color=method), na.rm=TRUE) +
        geom_point(aes(x=year, y=value, color=method, shape=method), na.rm=TRUE) +
        ylab(paste0("P", 100*p, "/average")) +
        scale_color_brewer(type="qual", palette="Set1") +
        theme_bw() + theme(
            legend.title = element_blank(),
            legend.position = "bottom",
            plot.background = element_rect(fill=plot_bg, color=plot_bg),
            panel.background = element_rect(fill=plot_bg),
            legend.key = element_rect(fill=plot_bg),
            text = element_text(color=plot_text_color)
        )

    # Create plot for top shares
    p2 <- ggplot(data_topshare) +
        geom_line(aes(x=year, y=value, color=method), na.rm=TRUE) +
        geom_point(aes(x=year, y=value, color=method, shape=method), na.rm=TRUE) +
        ylab(paste0("top ", 100*(1 - p), "% share")) +
        scale_y_continuous(labels=percent) +
        scale_color_brewer(type="qual", palette="Set1") +
        theme_bw() + theme(
            legend.position = "none",
            plot.background = element_rect(fill=plot_bg, color=plot_bg),
            panel.background = element_rect(fill=plot_bg),
            legend.key = element_rect(fill=plot_bg),
            text = element_text(color=plot_text_color)
        )

    # Extract the legend
    legend <- g_legend(p1)

    filename <- paste0("output/plots/comparison-time-series/time-series-recent-", iso, "-", income_type_short, "-", 100*p, ".pdf")
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
    country           = comparisons$country,
    iso               = comparisons$iso,
    year              = comparisons$year,
    income_type       = comparisons$income_type,
    income_type_short = comparisons$income_type_short,
    p                 = comparisons$p,

    threshold_m0 = relerr(comparisons$threshold_actual, comparisons$threshold_m0),
    threshold_m1 = relerr(comparisons$threshold_actual, comparisons$threshold_m1),
    threshold_m2 = relerr(comparisons$threshold_actual, comparisons$threshold_m2),
    threshold_m3 = relerr(comparisons$threshold_actual, comparisons$threshold_m3),
    threshold_m4 = relerr(comparisons$threshold_actual, comparisons$threshold_m4),

    topshare_m0 = relerr(comparisons$topshare_actual, comparisons$topshare_m0),
    topshare_m1 = relerr(comparisons$topshare_actual, comparisons$topshare_m1),
    topshare_m2 = relerr(comparisons$topshare_actual, comparisons$topshare_m2),
    topshare_m3 = relerr(comparisons$topshare_actual, comparisons$topshare_m3),
    topshare_m4 = relerr(comparisons$topshare_actual, comparisons$topshare_m4)
)

# Calculate the mean of the error
interpolation_mre <- ddply(interpolation_re, c("country", "iso", "income_type", "income_type_short", "p"), function(data) {
    return(data.frame(
        threshold_m0 = mean(data$threshold_m0),
        threshold_m1 = mean(data$threshold_m1),
        threshold_m2 = mean(data$threshold_m2),
        threshold_m3 = mean(data$threshold_m3),
        threshold_m4 = mean(data$threshold_m4),

        topshare_m0 = mean(data$topshare_m0),
        topshare_m1 = mean(data$topshare_m1),
        topshare_m2 = mean(data$topshare_m2),
        topshare_m3 = mean(data$topshare_m3),
        topshare_m4 = mean(data$topshare_m4)
    ))
})

# Generate LaTeX tables comparing the results
d_ply(interpolation_mre, "income_type", function(data) {
    filename <- paste0("output/tables/compare-interpolation/compare-", data$income_type_short[1], ".tex")

    # Remove the M4 method from the comparison (comment to keep)
    data$threshold_m4 <- NULL
    data$topshare_m4 <- NULL
    # Number of methods left to compare
    n <- (ncol(data) - 5)/2
    # Column for thresholds
    threshold_cols <- paste0("threshold_m", 0:(n - 1))
    topshare_cols <- paste0("topshare_m", 0:(n - 1))

    #print(data)
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
        cat("& & (ref.)")
        for (j in 1:(n - 1)) {
            cat(" & ($\\times ")
            cat(format(data_us[i, paste0("topshare_m", j)]/data_us[i, "topshare_m0"], digits=2, scientific=FALSE))
            cat("$)")
        }
        cat(paste0("\\\\ \\cmidrule(l){2-", 2 + n, "}\n"))
    }
    for (i in 1:3) {
        cat("& ")
        cat(sprintf("\\multirow{2}{*}{P%.0f/average} & ", 100*data_us$p[i]))
        cat(paste0(format(data_us[i, threshold_cols], digits=2, scientific=FALSE), collapse="\\% & "))
        cat("\\% \\\\ \n")
        cat("& & (ref.)")
        for (j in 1:(n - 1)) {
            cat(" & ($\\times ")
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
        cat("& & (ref.)")
        for (j in 1:(n - 1)) {
            cat(" & ($\\times ")
            cat(format(data_fr[i, paste0("topshare_m", j)]/data_fr[i, "topshare_m0"], digits=2, scientific=FALSE))
            cat("$)")
        }
        cat(paste0("\\\\ \\cmidrule(l){2-", 2 + n, "}\n"))
    }
    for (i in 1:3) {
        cat("& ")
        cat(sprintf("\\multirow{2}{*}{P%.0f/average} & ", 100*data_fr$p[i]))
        cat(paste0(format(data_fr[i, threshold_cols], digits=2, scientific=FALSE), collapse="\\% & "))
        cat("\\% \\\\ \n")
        cat("& & (ref.)")
        for (j in 1:(n - 1)) {
            cat(" & ($\\times ")
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

