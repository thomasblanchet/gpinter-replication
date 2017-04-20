# ---------------------------------------------------------------------------- #
# Code for the replication of the tables comparing the different extrapolation
# methods in Blanchet, Fournier & Piketty (2017).
# ---------------------------------------------------------------------------- #

# Perform the extrapolation with different methods
comparisons <- ddply(dina_data, c("iso", "country", "year", "income_type", "income_type_short"), function(data) {
    average <- data$average[1]
    income_type_short <- data$income_type_short[1]

    #p_in <- c(0, 0.5, 0.90, 0.99, 1)
    #p_out <- 0.999

    p_in <- c(0, 0.5, 0.8, 0.9, 1)
    p_out <- 0.99

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



