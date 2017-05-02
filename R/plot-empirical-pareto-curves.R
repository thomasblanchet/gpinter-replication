# ---------------------------------------------------------------------------- #
# “Generalized Pareto Curves: Theory and Applications”, 2017
# Thomas Blanchet, Juliette Fournier, Thomas Piketty
# ---------------------------------------------------------------------------- #
# Plot empirical Pareto curves using DINA data for France and the
# United States.
# ---------------------------------------------------------------------------- #

d_ply(dina_data, c("year", "income_type"), function(data) {
    year <- data$year[1]
    income_type <- data$income_type[1]
    income_type_short <- data$income_type_short[1]

    # Remove the bottom of the distribution and the very top
    data <- data[data$p >= 0.3e5 & data$p <= 0.999e5, ]

    # Rescale percentiles within [0, 1]
    data$p <- data$p/1e5

    # Range of the y-axis in the graph
    limits <- c(1.5, 4.5)

    cat(paste0("Plotting: empirical Pareto curves - ", income_type, " - ", year, "\n"))

    # Create plot
    filename <- paste0("output/plots/empirical-pareto-curves/pareto-curves-", income_type_short, "-", year, ".pdf")
    pdf(filename, family=plot_font, width=4.5, height=3.5)
    print(ggplot(data) +
        geom_line(aes(x=p, y=invpareto, linetype=country, color=country), na.rm=TRUE) +
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
