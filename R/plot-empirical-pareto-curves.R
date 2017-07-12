# ---------------------------------------------------------------------------- #
# “Generalized Pareto Curves: Theory and Applications”, 2017
# Thomas Blanchet, Juliette Fournier, Thomas Piketty
# ---------------------------------------------------------------------------- #
# Plot empirical Pareto curves using data from WID.world
# ---------------------------------------------------------------------------- #

# Plot the evolution of France and the United States on the same plot -------- #
d_ply(data_us_fr, .(widcode, year), function(data) {
    year <- data$year[1]
    var_code <- data$var_code[1]
    var_name <- data$var_name[1]
    pop_name <- data$pop_name[1]

    # Skip if we don't have both US & FR
    if (!("US" %in% data$iso && "FR" %in% data$iso)) {
        return(NULL)
    }

    # Remove the bottom of the distribution and the very top
    data <- data[data$p >= 0.3e5 & data$p <= 0.999e5, ]

    # Rescale percentiles within [0, 1]
    data$p <- data$p/1e5

    # Range of the y-axis in the graph
    if (var_code %in% c("ptinc", "pllin", "fiinc")) {
        limits <- c(1.5, 4.5)
    } else if (var_code %in% c("pkkin", "hweal")) {
        limits <- c(1.5, 20)
    }

    cat(paste0("* plotting: US & FR empirical Pareto curves - ", var_name, " - ", year, "\n"))

    # Create plot
    dir.create(paste0("output/plots/empirical-pareto-curves/us-fr/", var_code),
        showWarnings=FALSE, recursive=TRUE)
    filename <- paste0("output/plots/empirical-pareto-curves/us-fr/", var_code, "/", year, ".pdf")
    pdf(filename, family=plot_font, width=4.5, height=3.5)
    print(ggplot(data) +
        geom_line(aes(x=p, y=b, linetype=country, color=country), na.rm=TRUE) +
        scale_color_brewer(type="qual", palette="Set1") +
        scale_x_continuous(limits = c(0.3, 1), breaks=seq(0.3, 1, 0.1)) +
        scale_linetype_manual(values=c("United States"="solid", "France"="longdash")) +
        ylim(limits) +
        xlab(expression(paste("rank ", italic(p)))) +
        ylab(expression(paste("inverted Pareto coefficient ", italic(b), "(", italic(p), ")"))) +
        ggtitle(paste("Year", year)) +
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
})

# Plot Pareto curves for all countries in WID.world -------------------------- #

# NOTE: This is done for exploratory purposes. Most of the data currently in
# the World Wealth and Income Database was not constructed in a way that lead
# to correct generalized Pareto curves.

d_ply(data_wid, .(widcode, iso, year), function(data) {
    year <- data$year[1]
    country <- data$country[1]
    average <- data$average[1]
    var_code <- data$var_code[1]
    var_name <- data$var_name[1]
    pop_name <- data$pop_name[1]

    # Rescale percentiles within [0, 1]
    data$p <- data$p/1e5

    if (length(data$p) == 127) {
        # We have access to detailed g-percentile data: no need to interpolate
        df_plot <- data.frame(
            p = data$p[(data$p >= 0.3) & (data$p <= 0.999)],
            b = data$b[(data$p >= 0.3) & (data$p <= 0.999)]
        )

        df_points <- data.frame(p=numeric(0), b=numeric(0))
    } else if (length(data$p) >= 3) {
        # We only have a few data points: we interpolate
        dist <- tryCatch(tabulation_fit(data$p, data$threshold, average,
            topshare = data$top_share), error = function(e) {
                cat(paste0("* inconsistent input data: ", country, " - ", var_name, " - ", year, "\n"))
                return(NULL)
            }
        )
        # Skip if inconsistent data
        if (is.null(dist)) {
            return(NULL)
        }
        p_min <- max(min(data$p), 0.3)
        p_max <- max(data$p)
        p_all <- seq(p_min, p_max, length.out=1000)

        df_plot <- data.frame(
            p = p_all,
            b = invpareto(dist, p_all)
        )
        df_points <- data.frame(
            p = data$p[data$p >= 0.3],
            b = data$b[data$p >= 0.3]
        )
    } else {
        # Not enough data to interpolate
        cat(paste0("* not enough data: ", country, " - ", var_name, " - ", year, "\n"))
        return(NULL)
    }

    cat(paste0("* plotting: Pareto curve - ", country, " - ", var_name, " - ", year, "\n"))

    # Create plot
    dir.create(paste0("output/plots/empirical-pareto-curves/all-countries/", country),
        showWarnings=FALSE, recursive=TRUE)
    filename <- paste0("output/plots/empirical-pareto-curves/all-countries/",
        country, "/", var_code, "-", year, ".pdf")
    pdf(filename, family=plot_font, width=4.5, height=3.5)
    print(ggplot(df_plot) +
        geom_line(aes(x=p, y=b), na.rm=TRUE) +
        geom_point(data=df_points, aes(x=p, y=b), na.rm=TRUE) +
        ylim(c(1, 5)) + xlim(c(0.9, 1)) +
        xlab(expression(paste("rank ", italic(p)))) +
        ylab(expression(paste("inverted Pareto coefficient ", italic(b), "(", italic(p), ")"))) +
        ggtitle(paste("Year", year)) +
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
})

