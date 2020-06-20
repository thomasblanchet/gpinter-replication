# ---------------------------------------------------------------------------- #
# “Generalized Pareto Curves: Theory and Applications”, 2017
# Thomas Blanchet, Juliette Fournier, Thomas Piketty
# ---------------------------------------------------------------------------- #
# Compare the precision of the interpolation method with data.
# ---------------------------------------------------------------------------- #

data_us <- data_micro %>% filter(iso == "US" & var_code == "ptinc")

# ---------------------------------------------------------------------------- #
# Estimate top share with standard Pareto at the top
# ---------------------------------------------------------------------------- #

top_share_pareto1 <- function(x, p, modelstart) {
    n <- length(x)

    # Threshold at which the model starts
    thr <- quantile(x, modelstart)

    # Estimate Preto coef at the top via MLE
    alpha <- 1/mean(log(x[x > thr]/thr))
    beta <- alpha/(alpha - 1)

    if (p < modelstart) {
        # Mass between p and modelstart
        between <- (x >= quantile(x, p)) & (x < thr)
        mass1 <- sum(x[between])/n

        # Mass above modelstart
        mass2 <- thr*beta*(1 - modelstart)

        # Mass below p
        mass3 <- sum(x[x < quantile(x, p)])/n

        return((mass1 + mass2)/(mass1 + mass2 + mass3))
    } else if (p > modelstart) {
        # Mass between p and modelstart
        between <- (x < quantile(x, p)) & (x >= thr)
        mass1 <- sum(x[between])/n

        # Mass above p
        mass2 <- quantile(x, p)*beta*(1 - p)

        # Mass below modelstart
        mass3 <- sum(x[x < thr])/n

        return(mass2/(mass1 + mass2 + mass3))
    } else {
        # Mass above p/modelstart
        mass1 <- thr*beta*(1 - p)

        # Mass below p/modelstart
        mass2 <- sum(x[x < thr])/n

        return(mass1/(mass1 + mass2))
    }
}

# ---------------------------------------------------------------------------- #
# Estimate top share with generalized Pareto at the top
# ---------------------------------------------------------------------------- #

top_share_gpd <- function(x, p, modelstart) {
    n <- length(x)

    # Threshold at which the model starts
    mu <- quantile(x, modelstart)

    # Estimate Pareto mode at the top via MLE
    y <- x[x >= mu]

    ll <- function(xi, sigma) {
        return(sum(-(1/xi + 1)*log(1 + xi*(y - mu)/sigma) - log(sigma)))
    }
    ll_dxi <- function(xi, sigma) {
        return(sum(((-1 - 1/xi)*(-mu + y))/(sigma*(1 + (xi*(-mu + y))/sigma)) +
            log(1 + (xi*(-mu + y))/sigma)/xi^2))
    }
    ll_dsigma <- function(xi, sigma) {
        return(sum(-(1/sigma) - ((-1 - 1/xi)*xi*(-mu + y))/(sigma^2*(1 + (xi*(-mu + y))/sigma))))
    }

    a <- mean(y)
    b <- mean(y*(1 - ecdf(y)(y)))

    xi0 <- (a - 4*b + mu)/(a - 2*b)
    sigma0 <- ((a - mu)*(2*b - mu))/(a - 2*b)

    fit <- optim(
        par = c(qlogis(xi0), log(sigma0)),
        fn = function(par) -ll(plogis(par[1]), exp(par[2])),
        gr = function(par) -c(
            dlogis(par[1])*ll_dxi(plogis(par[1]), exp(par[2])),
            exp(par[2])*ll_dsigma(plogis(par[1]), exp(par[2]))
        ),
        method = "CG"
    )

    xi <- plogis(fit$par[1])
    sigma <- exp(fit$par[2])

    if (p < modelstart) {
        # Mass between p and modelstart
        between <- (x >= quantile(x, p)) & (x < mu)
        mass1 <- sum(x[between])/n

        # Mass above modelstart
        p0 <- modelstart
        p1 <- p
        mass2 <- mu - mu*p1 + ((-1 + p1)*sigma*(-1 + ((-1 + p0)/(-1 + p1))^xi + xi))/(-1 + xi)

        # Mass below p
        mass3 <- sum(x[x < quantile(x, p)])/n

        return((mass1 + mass2)/(mass1 + mass2 + mass3))
    } else if (p > modelstart) {
        p0 <- modelstart
        p1 <- p

        # Mass between p and modelstart
        mass1 <- -(mu*p0) + mu*p1 + p0*sigma - p1*sigma +
            ((-1 + p0 - ((-1 + p0)/(-1 + p1))^xi*(-1 + p1))*sigma)/(-1 + xi)

        # Mass above p
        mass2 <- mu - mu*p1 + ((-1 + p1)*sigma*(-1 + ((-1 + p0)/(-1 + p1))^xi + xi))/(-1 + xi)

        # Mass below modelstart
        mass3 <- sum(x[x < mu])/n

        return(mass2/(mass1 + mass2 + mass3))
    } else {
        # Mass above p/modelstart
        p0 <- modelstart
        p1 <- p
        mass1 <- mu - mu*p1 + ((-1 + p1)*sigma*(-1 + ((-1 + p0)/(-1 + p1))^xi + xi))/(-1 + xi)

        # Mass below p/modelstart
        mass2 <- sum(x[x < mu])/n

        return(mass1/(mass1 + mass2))
    }
}

# ---------------------------------------------------------------------------- #
# Single draw from a sample n = 1e5
# ---------------------------------------------------------------------------- #

set.seed(19920902)
cmp_top_single <- data_us %>% ddply(.(year), function(data) {
    #if (data$year[1] %% 10 != 0) {
    #    return(NULL)
    #}

    p_long <- c(0, seq(0.1, 0.9, 0.1), seq(0.91, 0.99, 0.01),
        seq(0.991, 0.999, 0.001), 1)
    p_short1 <- c(0, 0.5, 0.8, 0.9999, 1)
    p_short2 <- c(0, 0.5, 0.9, 0.95, 1)

    # Long tabulation (for simulating the distribution)
    data$bracket <- cut(data$p,
        breaks = round(p_long*1e5),
        include.lowest = TRUE,
        labels = FALSE,
        right = FALSE
    )

    long_tab <- ddply(data, .(bracket), function(data) {
        return(data.frame(
            p = min(data$p)/1e5,
            threshold = data$threshold[which.min(data$p)],
            top_share = data$top_share[which.min(data$p)]
        ))
    })
    long_tab <- long_tab[-1, ]

    # Short tabulation (for testing the interpolation method)
    data$bracket <- cut(data$p,
        breaks = round(p_short1*1e5),
        include.lowest = TRUE,
        labels = FALSE,
        right = FALSE
    )

    short_tab1 <- ddply(data, .(bracket), function(data) {
        return(data.frame(
            p = min(data$p)/1e5,
            threshold = data$threshold[which.min(data$p)],
            top_share = data$top_share[which.min(data$p)]
        ))
    })
    short_tab1 <- short_tab1[-1, ]

    data$bracket <- cut(data$p,
        breaks = round(p_short2*1e5),
        include.lowest = TRUE,
        labels = FALSE,
        right = FALSE
    )

    short_tab2 <- ddply(data, .(bracket), function(data) {
        return(data.frame(
            p = min(data$p)/1e5,
            threshold = data$threshold[which.min(data$p)],
            top_share = data$top_share[which.min(data$p)]
        ))
    })
    short_tab2 <- short_tab2[-1, ]

    average <- data$average[1]

    dist <- tabulation_fit(long_tab$p, long_tab$threshold,
        average, topshare = long_tab$top_share)
    top1 <- top_share(dist, 0.99)
    top01 <- top_share(dist, 0.999)

    dist_short1 <- tabulation_fit(short_tab1$p, short_tab1$threshold,
        average, topshare = short_tab1$top_share)
    dist_short2 <- tabulation_fit(short_tab2$p, short_tab2$threshold,
        average, topshare = short_tab2$top_share)

    # Simulate a large sample and calcualte the top share from it
    print(data$year[1])
    y <- simulate_gpinter(dist, 1e5)
    thr <- quantile(y, c(0.99, 0.999))
    top1_simul <- sum(y[y > thr[1]])/sum(y)
    top01_simul <- sum(y[y > thr[2]])/sum(y)

    top1_pareto1_top5 <- top_share_pareto1(y, 0.99, 0.95)
    top01_pareto1_top5 <- top_share_pareto1(y, 0.999, 0.95)

    top1_pareto1_top1 <- top_share_pareto1(y, 0.99, 0.99)
    top01_pareto1_top1 <- top_share_pareto1(y, 0.999, 0.99)

    top1_gpd_top5 <- top_share_gpd(y, 0.99, 0.95)
    top01_gpd_top5 <- top_share_gpd(y, 0.999, 0.95)

    top1_gpd_top1 <- top_share_gpd(y, 0.99, 0.99)
    top01_gpd_top1 <- top_share_gpd(y, 0.999, 0.99)

    top1_inter <- top_share(dist_short1, 0.99)
    top01_inter <- top_share(dist_short1, 0.999)

    top1_extra <- top_share(dist_short2, 0.99)
    top01_extra <- top_share(dist_short2, 0.999)

    return(tibble(
        p = c(99, 99.9),
        top = c(top1, top01),
        top_inter = c(top1_inter, top01_inter),
        top_extra = c(top1_extra, top01_extra),
        top_simul = c(top1_simul, top01_simul),
        top_pareto1_top5 = c(top1_pareto1_top5, top01_pareto1_top5),
        top_pareto1_top1 = c(top1_pareto1_top1, top01_pareto1_top1),
        top_gpd_top5 = c(top1_gpd_top5, top01_gpd_top5),
        top_gpd_top1 = c(top1_gpd_top1, top01_gpd_top1)
    ))
})

filename <- file.path("output", "plots", "compare-survey", "compare-survey-single-inter.pdf")
pdf(filename, family = plot_font, width = 4.5, height = 4)

cmp_top_single %>% ggplot() +
    geom_line(aes(x = year, y = top, group = p, color = "actual value")) +
    geom_point(aes(x = year, y = top, group = p, color = "actual value"), shape = 19) +
    geom_line(aes(x = year, y = top_inter, group = p, color = "interpolated value based on p = 50%, 80% and 99.99%")) +
    geom_point(aes(x = year, y = top_inter, group = p, color = "interpolated value based on p = 50%, 80% and 99.99%"), shape = 0) +
    geom_line(aes(x = year, y = top_simul, group = p, color = "large random sample (single draw, n = 100 000)")) +
    geom_point(aes(x = year, y = top_simul, group = p, color = "large random sample (single draw, n = 100 000)"), shape = 2) +
    scale_color_manual(
        values = c(
            "actual value" = "#e41a1c",
            "interpolated value based on p = 50%, 80% and 99.99%" = "#377eb8",
            "large random sample (single draw, n = 100 000)" = "#4daf4a"
        ),
        guide = guide_legend(
            override.aes = list(shape = c(19, 0, 2))
        )
    ) +
    xlab("year") + ylab("top income share") +
    scale_y_continuous(labels = percent) +
    geom_label(aes(x = 1975, y = 0.1275, label = "top 1%")) +
    geom_label(aes(x = 1975, y = 0.055, label = "top 0.1%")) +
    theme_bw() +
    theme(
        legend.background = element_rect(fill = plot_bg),
        legend.title = element_blank(),
        legend.margin = margin(t = -0.5, unit = 'cm'),
        legend.box.margin = margin(0, 0, 0, 0),
        legend.direction = "vertical",
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        plot.background = element_rect(fill = plot_bg, color = plot_bg),
        panel.background = element_rect(fill = plot_bg),
        legend.key = element_rect(fill = plot_bg),
        text = element_text(color = plot_text_color)
    )
dev.off()
embed_fonts(path.expand(filename))

filename <- file.path("output", "plots", "compare-survey", "compare-survey-single-inter-model.pdf")
pdf(filename, family = plot_font, width = 4.5, height = 4)

cmp_top_single %>% ggplot() +
    geom_line(aes(x = year, y = top, group = p, color = "actual value")) +
    geom_point(aes(x = year, y = top, group = p, color = "actual value"), shape = 19) +
    geom_line(aes(x = year, y = top_pareto1_top5, group = p, color = "large random sample + Pareto model for top 5%")) +
    geom_point(aes(x = year, y = top_pareto1_top5, group = p, color = "large random sample + Pareto model for top 5%"), shape = 2) +
    geom_line(aes(x = year, y = top_pareto1_top1, group = p, color = "large random sample + Pareto model for top 1%")) +
    geom_point(aes(x = year, y = top_pareto1_top1, group = p, color = "large random sample + Pareto model for top 1%"), shape = 2) +
    scale_color_manual(
        values = c(
            "actual value" = "#e41a1c",
            "large random sample + Pareto model for top 5%" = "#377eb8",
            "large random sample + Pareto model for top 1%" = "#4daf4a"
        ),
        guide = guide_legend(
            override.aes = list(shape = c(19, 0, 2))
        )
    ) +
    xlab("year") + ylab("top income share") +
    scale_y_continuous(labels = percent) +
    geom_label(aes(x = 1975, y = 0.1275, label = "top 1%")) +
    geom_label(aes(x = 1975, y = 0.055, label = "top 0.1%")) +
    theme_bw() +
    theme(
        legend.background = element_rect(fill = plot_bg),
        legend.title = element_blank(),
        legend.margin = margin(t = -0.5, unit = 'cm'),
        legend.box.margin = margin(0, 0, 0, 0),
        legend.direction = "vertical",
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        plot.background = element_rect(fill = plot_bg, color = plot_bg),
        panel.background = element_rect(fill = plot_bg),
        legend.key = element_rect(fill = plot_bg),
        text = element_text(color = plot_text_color)
    )
dev.off()
embed_fonts(path.expand(filename))

filename <- file.path("output", "plots", "compare-survey", "compare-survey-single-inter-gpd.pdf")
pdf(filename, family = plot_font, width = 4.5, height = 4)

cmp_top_single %>% ggplot() +
    geom_line(aes(x = year, y = top, group = p, color = "actual value")) +
    geom_point(aes(x = year, y = top, group = p, color = "actual value"), shape = 19) +
    geom_line(aes(x = year, y = top_gpd_top5, group = p, color = "large random sample + GPD model for top 5%")) +
    geom_point(aes(x = year, y = top_gpd_top5, group = p, color = "large random sample + GPD model for top 5%"), shape = 2) +
    geom_line(aes(x = year, y = top_gpd_top1, group = p, color = "large random sample + GPD model for top 1%")) +
    geom_point(aes(x = year, y = top_gpd_top1, group = p, color = "large random sample + GPD model for top 1%"), shape = 2) +
    scale_color_manual(
        values = c(
            "actual value" = "#e41a1c",
            "large random sample + GPD model for top 5%" = "#377eb8",
            "large random sample + GPD model for top 1%" = "#4daf4a"
        ),
        guide = guide_legend(
            override.aes = list(shape = c(19, 0, 2))
        )
    ) +
    xlab("year") + ylab("top income share") +
    scale_y_continuous(labels = percent) +
    geom_label(aes(x = 1975, y = 0.1275, label = "top 1%")) +
    geom_label(aes(x = 1975, y = 0.055, label = "top 0.1%")) +
    theme_bw() +
    theme(
        legend.background = element_rect(fill = plot_bg),
        legend.title = element_blank(),
        legend.margin = margin(t = -0.5, unit = 'cm'),
        legend.box.margin = margin(0, 0, 0, 0),
        legend.direction = "vertical",
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        plot.background = element_rect(fill = plot_bg, color = plot_bg),
        panel.background = element_rect(fill = plot_bg),
        legend.key = element_rect(fill = plot_bg),
        text = element_text(color = plot_text_color)
    )
dev.off()
embed_fonts(path.expand(filename))


filename <- file.path("output", "plots", "compare-survey", "compare-survey-single-extra-10m.pdf")
pdf(filename, family = plot_font, width = 4.5, height = 4)

cmp_top_single %>% ggplot() +
    geom_line(aes(x = year, y = top, group = p, color = "actual value")) +
    geom_point(aes(x = year, y = top, group = p, color = "actual value"), shape = 19) +
    geom_line(aes(x = year, y = top_extra, group = p, color = "extrapolated value based on p = 50%, 90% and 95%")) +
    geom_point(aes(x = year, y = top_extra, group = p, color = "extrapolated value based on p = 50%, 90% and 95%"), shape = 0) +
    geom_line(aes(x = year, y = top_simul, group = p, color = "large random sample (single draw, n = 10 000 000)")) +
    geom_point(aes(x = year, y = top_simul, group = p, color = "large random sample (single draw, n = 10 000 000)"), shape = 2) +
    scale_color_manual(
        values = c(
            "actual value" = "#e41a1c",
            "extrapolated value based on p = 50%, 90% and 95%" = "#377eb8",
            "large random sample (single draw, n = 10 000 000)" = "#4daf4a"
        ),
        guide = guide_legend(
            override.aes = list(shape = c(19, 0, 2))
        )
    ) +
    xlab("year") + ylab("top income share") +
    scale_y_continuous(labels = percent) +
    geom_label(aes(x = 1975, y = 0.1275, label = "top 1%")) +
    geom_label(aes(x = 1975, y = 0.055, label = "top 0.1%")) +
    theme_bw() +
    theme(
        legend.background = element_rect(fill = plot_bg),
        legend.title = element_blank(),
        legend.margin = margin(t = -0.5, unit = 'cm'),
        legend.box.margin = margin(0, 0, 0, 0),
        legend.direction = "vertical",
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        plot.background = element_rect(fill = plot_bg, color = plot_bg),
        panel.background = element_rect(fill = plot_bg),
        legend.key = element_rect(fill = plot_bg),
        text = element_text(color = plot_text_color)
    )
dev.off()
embed_fonts(path.expand(filename))

# ---------------------------------------------------------------------------- #
# Confidence interval from samples n = 1e5
# ---------------------------------------------------------------------------- #

cl <- makeForkCluster(6)
#cl <- NULL

set.seed(19920902)
cmp_top <- data_us %>% ddply(.(year), function(data) {
    # if (data$year[1] %% 10 != 0) {
    #     return(NULL)
    # }

    p_long <- c(0, seq(0.1, 0.9, 0.1), seq(0.91, 0.99, 0.01),
        seq(0.991, 0.999, 0.001), 1)
    p_short1 <- c(0, 0.5, 0.8, 0.9999, 1)
    p_short2 <- c(0, 0.5, 0.9, 0.95, 1)

    # Long tabulation (for simulating the distribution)
    data$bracket <- cut(data$p,
        breaks = round(p_long*1e5),
        include.lowest = TRUE,
        labels = FALSE,
        right = FALSE
    )

    long_tab <- ddply(data, .(bracket), function(data) {
        return(data.frame(
            p = min(data$p)/1e5,
            threshold = data$threshold[which.min(data$p)],
            top_share = data$top_share[which.min(data$p)]
        ))
    })
    long_tab <- long_tab[-1, ]

    # Short tabulation (for testing the interpolation method)
    data$bracket <- cut(data$p,
        breaks = round(p_short1*1e5),
        include.lowest = TRUE,
        labels = FALSE,
        right = FALSE
    )

    short_tab1 <- ddply(data, .(bracket), function(data) {
        return(data.frame(
            p = min(data$p)/1e5,
            threshold = data$threshold[which.min(data$p)],
            top_share = data$top_share[which.min(data$p)]
        ))
    })
    short_tab1 <- short_tab1[-1, ]

    data$bracket <- cut(data$p,
        breaks = round(p_short2*1e5),
        include.lowest = TRUE,
        labels = FALSE,
        right = FALSE
    )

    short_tab2 <- ddply(data, .(bracket), function(data) {
        return(data.frame(
            p = min(data$p)/1e5,
            threshold = data$threshold[which.min(data$p)],
            top_share = data$top_share[which.min(data$p)]
        ))
    })
    short_tab2 <- short_tab2[-1, ]

    average <- data$average[1]

    dist <- tabulation_fit(long_tab$p, long_tab$threshold,
        average, topshare = long_tab$top_share)
    top1 <- top_share(dist, 0.99)
    top01 <- top_share(dist, 0.999)

    dist_short1 <- tabulation_fit(short_tab1$p, short_tab1$threshold,
        average, topshare = short_tab1$top_share)
    dist_short2 <- tabulation_fit(short_tab2$p, short_tab2$threshold,
        average, topshare = short_tab2$top_share)

    # Simulate a large sample and calcualte the top share from it
    print(data$year[1])
    top_simul <- pbsapply(1:1e3, function(i) {
        y <- simulate_gpinter(dist, 1e5)
        thr <- quantile(y, c(0.99, 0.999))
        top1_simul <- sum(y[y > thr[1]])/sum(y)
        top01_simul <- sum(y[y > thr[2]])/sum(y)
        return(c(top1_simul, top01_simul))
    }, cl = cl)
    top1_simul_lb <- quantile(top_simul[1, ], 0.05)
    top1_simul_ub <- quantile(top_simul[1, ], 0.95)
    top01_simul_lb <- quantile(top_simul[2, ], 0.05)
    top01_simul_ub <- quantile(top_simul[2, ], 0.95)

    top1_inter <- top_share(dist_short1, 0.99)
    top01_inter <- top_share(dist_short1, 0.999)

    top1_extra <- top_share(dist_short2, 0.99)
    top01_extra <- top_share(dist_short2, 0.999)

    return(tibble(
        p = c(99, 99.9),
        top = c(top1, top01),
        top_inter = c(top1_inter, top01_inter),
        top_extra = c(top1_extra, top01_extra),
        top_simul_lb = c(top1_simul_lb, top01_simul_lb),
        top_simul_ub = c(top1_simul_ub, top01_simul_ub)
    ))
})

stopCluster(cl)

filename <- file.path("output", "plots", "compare-survey", "compare-survey-ci-inter.pdf")
pdf(filename, family = plot_font, width = 4.5, height = 4)

cmp_top %>% ggplot() +
    geom_line(aes(x = year, y = top, group = p, color = "actual value")) +
    geom_point(aes(x = year, y = top, group = p, color = "actual value"), shape = 19) +
    geom_line(aes(x = year, y = top_inter, group = p, color = "interpolated value based on p = 50%, 80% and 99.99%")) +
    geom_point(aes(x = year, y = top_inter, group = p, color = "interpolated value based on p = 50%, 80% and 99.99%"), shape = 0) +
    geom_line(aes(x = year, y = top_simul_lb, group = p, color = "large random sample (90% CI, n = 100 000)"), linetype = "longdash") +
    geom_line(aes(x = year, y = top_simul_ub, group = p, color = "large random sample (90% CI, n = 100 000)"), linetype = "longdash") +
    scale_color_manual(
        values = c(
            "actual value" = "#e41a1c",
            "interpolated value based on p = 50%, 80% and 99.99%" = "#377eb8",
            "large random sample (90% CI, n = 100 000)" = "#4daf4a"
        ),
        guide = guide_legend(
            override.aes = list(shape = c(19, 0, NA), linetype = c("solid", "solid", "longdash"))
        )
    ) +
    xlab("year") + ylab("top income share") +
    scale_y_continuous(labels = percent) +
    geom_label(aes(x = 1975, y = 0.14, label = "top 1%")) +
    geom_label(aes(x = 1975, y = 0.07, label = "top 0.1%")) +
    theme_bw() +
    theme(
        legend.background = element_rect(fill = plot_bg),
        legend.title = element_blank(),
        legend.margin = margin(t = -0.5, unit = 'cm'),
        legend.box.margin = margin(0, 0, 0, 0),
        legend.direction = "vertical",
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        plot.background = element_rect(fill = plot_bg, color = plot_bg),
        panel.background = element_rect(fill = plot_bg),
        legend.key = element_rect(fill = plot_bg),
        text = element_text(color = plot_text_color)
    )
dev.off()
embed_fonts(path.expand(filename))

filename <- file.path("output", "plots", "compare-survey", "compare-survey-ci-extra.pdf")
pdf(filename, family = plot_font, width = 4.5, height = 4)

cmp_top %>% ggplot() +
    geom_line(aes(x = year, y = top, group = p, color = "actual value")) +
    geom_point(aes(x = year, y = top, group = p, color = "actual value"), shape = 19) +
    geom_line(aes(x = year, y = top_extra, group = p, color = "extrapolated value based on p = 50%, 90% and 95%")) +
    geom_point(aes(x = year, y = top_extra, group = p, color = "extrapolated value based on p = 50%, 90% and 95%"), shape = 0) +
    geom_line(aes(x = year, y = top_simul_lb, group = p, color = "large random sample (90% CI, n = 100 000)"), linetype = "longdash") +
    geom_line(aes(x = year, y = top_simul_ub, group = p, color = "large random sample (90% CI, n = 100 000)"), linetype = "longdash") +
    scale_color_manual(
        values = c(
            "actual value" = "#e41a1c",
            "extrapolated value based on p = 50%, 90% and 95%" = "#377eb8",
            "large random sample (90% CI, n = 100 000)" = "#4daf4a"
        ),
        guide = guide_legend(
            override.aes = list(shape = c(19, 0, NA), linetype = c("solid", "solid", "longdash"))
        )
    ) +
    xlab("year") + ylab("top income share") +
    scale_y_continuous(labels = percent) +
    geom_label(aes(x = 1975, y = 0.14, label = "top 1%")) +
    geom_label(aes(x = 1975, y = 0.07, label = "top 0.1%")) +
    theme_bw() +
    theme(
        legend.background = element_rect(fill = plot_bg),
        legend.title = element_blank(),
        legend.margin = margin(t = -0.5, unit = 'cm'),
        legend.box.margin = margin(0, 0, 0, 0),
        legend.direction = "vertical",
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        plot.background = element_rect(fill = plot_bg, color = plot_bg),
        panel.background = element_rect(fill = plot_bg),
        legend.key = element_rect(fill = plot_bg),
        text = element_text(color = plot_text_color)
    )
dev.off()
embed_fonts(path.expand(filename))

# ---------------------------------------------------------------------------- #
# Calculate typical survey error in a table
# ---------------------------------------------------------------------------- #

data_us_2010 <- subset(data_micro, year == "2010" & iso == "US" & var_code == "ptinc")

dist <- tabulation_fit(
    p = data_us_2010$p/1e5,
    average = data_us_2010$average[1],
    threshold = data_us_2010$threshold,
    topshare = data_us_2010$top_share,
    fast = TRUE
)

p_in <- c(0, 0.1, 0.5, 0.9, 0.99, 1)
p_out <- c(0.30, 0.5, 0.75, 0.9, 0.95, 0.99, 0.999)

data_us_2010$bracket <- cut(data_us_2010$p,
    breaks = p_in*1e5,
    include.lowest = TRUE,
    labels = FALSE,
    right = FALSE
)

short_tab <- ddply(data_us_2010, .(bracket), function(data) {
    return(data.frame(
        p = min(data$p)/1e5,
        threshold = data$threshold[which.min(data$p)],
        top_share = data$top_share[which.min(data$p)]
    ))
})

average <- data_us_2010$average[1]
short_tab <- short_tab[-1, ]

dist <- tabulation_fit(short_tab$p, short_tab$threshold, average, topshare=short_tab$top_share)

# Simulate that distribution
set.seed(19920902)

n <- 100e6
x <- simulate_gpinter(dist, n)

x <- sort(x)
top_share_x <- rev(cumsum(rev(x)))
top_share_x <- top_share_x/top_share_x[1]

top_share_test <- top_share_x[floor(n*p_out)]
thresholds_test <- x[floor(n*p_out)]

# Simulate survey samples
k_all <- c(1e3, 1e4, 1e5, 1e6, 1e7)
replications <- pbreplicate(1000, {
    error_top_share <- matrix(nrow=length(k_all), ncol=length(p_out))
    error_thresholds <- matrix(nrow=length(k_all), ncol=length(p_out))
    for (i in 1:length(k_all)) {
        k <- k_all[i]

        y <- sample(x, k)
        y <- sort(y)

        top_share_y <- rev(cumsum(rev(y)))
        top_share_y <- top_share_y/top_share_y[1]
        top_share_y <- top_share_y[floor(k*p_out)]

        thresholds_y <- y[floor(k*p_out)]

        error_top_share[i, ] <- (top_share_y - top_share_test)/top_share_test
        error_thresholds[i, ] <- (thresholds_y - thresholds_test)/thresholds_test
    }
    return(list(error_top_share, error_thresholds))
}, simplify=FALSE)

mean_error_top_share <- matrix(data=0, nrow=length(p_out), ncol=length(k_all) + 1)
for (i in 1:length(k_all)) {
    mean_error_top_share[, i] <- 100*rowMeans(sapply(replications, function(rep) abs(rep[[1]][i, ])))
}

mean_error_thresholds <- matrix(data=0, nrow=length(p_out), ncol=length(k_all) + 1)
for (i in 1:length(k_all)) {
    mean_error_thresholds[, i] <- 100*rowMeans(sapply(replications, function(rep) abs(rep[[2]][i, ])))
}

# Generate LaTeX table
n <- length(k_all) + 1

dir.create("output/tables", showWarnings=FALSE, recursive=TRUE)
filename <- paste0("output/tables/compare-survey.tex")

sink(filename)

cat("\\begin{tabular}{rP@{}P@{}P@{}P@{}P@{}P@{}}\\toprule\n")

cat(paste0("& \\multicolumn{", n, "}{>{\\centering\\arraybackslash}p{", 2*n, "cm}}"))
cat("{mean percentage gap between estimated and observed values ")
cat("for a survey with simple random sampling and sample size $n$}")
cat(paste0("\\\\ \\cmidrule(l){2-", 2 + n - 1, "}\n"))

for (j in 1:n) {
    k <- ifelse(is.na(k_all[j]), 1e8, k_all[j])
    cat(paste0(" & $n=10^{", log10(k), "}$"))
}
cat("\\\\ \\midrule\n")

for (i in 1:length(p_out)) {
    cat(paste0("Top ", 100*(1 - p_out[i]), "\\% share"))
    for (j in 1:n) {
        cat(sprintf(" & %.2f\\%%", mean_error_top_share[i, j]))
    }
    cat("\\\\\n")
}
cat("\\midrule\n")

for (i in 1:length(p_out)) {
    cat(paste0("P", 100*p_out[i], " threshold"))
    for (j in 1:n) {
        cat(sprintf(" & %.2f\\%%", mean_error_thresholds[i, j]))
    }
    cat("\\\\\n")
}
cat("\\bottomrule\n")

cat("\\end{tabular}\n")

sink()
