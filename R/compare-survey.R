# ---------------------------------------------------------------------------- #
# “Generalized Pareto Curves: Theory and Applications”, 2017
# Thomas Blanchet, Juliette Fournier, Thomas Piketty
# ---------------------------------------------------------------------------- #
# Compare the precision of the interpolation method with data.
# ---------------------------------------------------------------------------- #

# Estimate the distribution of pre-tax national income in the US in 2014
data_us_2010 <- subset(data_micro, year == "2010" & iso == "US" & var_code == "ptinc")

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
