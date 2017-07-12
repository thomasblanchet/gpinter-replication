# ---------------------------------------------------------------------------- #
# “Generalized Pareto Curves: Theory and Applications”, 2017
# Thomas Blanchet, Juliette Fournier, Thomas Piketty
# ---------------------------------------------------------------------------- #
# Plot the misspecification error bounds for a fixed phi'''.
# ---------------------------------------------------------------------------- #

# Parameters of the tabulation in input
pk <- c(0.1, 0.5, 0.9, 0.99)
xk <- -log(1 - pk)
x <- seq(min(xk), max(xk), length.out=1e3)

# Calculate the error bounds
df <- data.frame(
    x = x,
    ey = interpolation_value_error_bound_cons(x, xk, 1),
    edy = interpolation_deriv_error_bound_cons(x, xk, 1)
)

dir.create("output/plots/misspecification-error-bound", showWarnings=FALSE, recursive=TRUE)
filename <- "output/plots/misspecification-error-bound/error-bound-phi.pdf"
pdf(filename, family=plot_font, width=4.5, height=3.5)
print(ggplot(df) + geom_line(aes(x=x, y=ey)) +
    geom_vline(xintercept=xk, linetype="dashed") +
    annotate("text", label="paste(italic(p) == 10, '%')", x=xk[1] + 0.16, y=0.125, angle=90, parse=TRUE) +
    annotate("text", label="paste(italic(p) == 50, '%')", x=xk[2] + 0.16, y=0.125, angle=90, parse=TRUE) +
    annotate("text", label="paste(italic(p) == 90, '%')", x=xk[3] + 0.16, y=0.125, angle=90, parse=TRUE) +
    annotate("text", label="paste(italic(p) == 99, '%')", x=xk[4] - 0.16, y=0.125, angle=90, parse=TRUE) +
    xlab(expression(paste(italic(x) == -log, "(", 1 - italic(p), ")"))) + ylab(expression(paste("multiple of ||", phi1, "'''||"[infinity]))) +
    ggtitle(expression(paste("Error bound on ", phi1, "(", italic(x), ")")),
        subtitle=expression("for a tabulation with"~italic(p)~"= 10%, 50%, 90% and 99%")) + theme_bw() +
    theme(
        plot.title = element_text(hjust=0.5),
        plot.subtitle = element_text(hjust=0.5),
        plot.background = element_rect(fill=plot_bg, color=plot_bg),
        panel.background = element_rect(fill=plot_bg),
        legend.key = element_rect(fill=plot_bg),
        text = element_text(color=plot_text_color)
    )
)
dev.off()
embed_fonts(path.expand(filename))

filename <- "output/plots/misspecification-error-bound/error-bound-deriv-phi.pdf"
pdf(filename, family=plot_font, width=4.5, height=3.5)
print(ggplot(df) + geom_line(aes(x=x, y=edy)) +
    geom_vline(xintercept=xk, linetype="dashed") + ylim(c(0, 0.25)) +
    annotate("text", label="paste(italic(p) == 10, '%')", x=xk[1] + 0.16, y=0.22, angle=90, parse=TRUE) +
    annotate("text", label="paste(italic(p) == 50, '%')", x=xk[2] + 0.16, y=0.22, angle=90, parse=TRUE) +
    annotate("text", label="paste(italic(p) == 90, '%')", x=xk[3] + 0.16, y=0.22, angle=90, parse=TRUE) +
    annotate("text", label="paste(italic(p) == 99, '%')", x=xk[4] - 0.16, y=0.22, angle=90, parse=TRUE) +
    xlab(expression(paste(italic(x) == -log, "(", 1 - italic(p), ")"))) + ylab(expression(paste("multiple of ||", phi1, "'''||"[infinity]))) +
    ggtitle(expression(paste("Error bound on ", phi1, "'(", italic(x), ")")),
        subtitle=expression("for a tabulation with"~italic(p)~"= 10%, 50%, 90% and 99%")) + theme_bw() +
    theme(
        plot.title = element_text(hjust=0.5),
        plot.subtitle = element_text(hjust=0.5),
        plot.background = element_rect(fill=plot_bg, color=plot_bg),
        panel.background = element_rect(fill=plot_bg),
        legend.key = element_rect(fill=plot_bg),
        text = element_text(color=plot_text_color)
    )
)
dev.off()
embed_fonts(path.expand(filename))
