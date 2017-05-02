# ---------------------------------------------------------------------------- #
# “Generalized Pareto Curves: Theory and Applications”, 2017
# Thomas Blanchet, Juliette Fournier, Thomas Piketty
# ---------------------------------------------------------------------------- #
# Provide out-of-sample estimates of error bounds.
# ---------------------------------------------------------------------------- #

# On the most recent ten years: calculate the envelope of phi'''
envelope_phi_d3 <- ddply(subset(dina_data, income_type_short == "pretax"), c("iso", "country"), function(data) {
    # Only keep the most recent ten years
    yearmax <- max(data$year)
    data <- subset(data, year > yearmax - 10)

    # Estimte phi''' in each year
    phi_d3_estim <- ddply(data, "year", function(data) {
        data <- data[data$p >= 10000 & data$p <= 99000, ]

        p <- data$p/1e5
        x <- -log(1 - data$p/1e5)

        x_out <- seq(-log(1 - min(p)), -log(1 - max(p)), length.out=1000)
        d3_phi <- sapply(x_out, function(x0) lpoly(x0, x, data$phi, data$dphi, 0.05))

        return(data.frame(x=x_out, d3_phi=d3_phi))
    })

    # Take the max to get the envelope
    return(ddply(phi_d3_estim, "x", function(data) {
        return(data.frame(x=data$x[1], max_phi_d3=max(abs(data$d3_phi))))
    }))
})

# Calculate an error bound implied by this envelope
err_bound_theory <- ddply(envelope_phi_d3, c("iso", "country"), function(data) {
    # Create a function for the envelope
    phi_d3_fun <- splinefun(data$x, data$max_phi_d3, method="monoH.FC")

    # Tabulation to consider
    p_in <- c(0.1, 0.5, 0.9, 0.99)
    p_out <- c(0.3, 0.75, 0.95)
    x_in <- -log(1 - p_in)
    x_out <- -log(1 - p_out)

    # Calculate the error on phi and phi'
    err_phi <- interpolation_value_error_bound_noncons(x_out, x_in, phi_d3_fun)
    err_dphi <- interpolation_deriv_error_bound_noncons(x_out, x_in, phi_d3_fun)

    return(data.frame(p=p_out, err_phi=err_phi, err_dphi=err_dphi))
})

# Calculate the actual error over the least recent period
err_actual <- ddply(subset(dina_data, income_type_short == "pretax"), c("iso", "country", "year"), function(data) {
    # Only keep the least recent years
    yearmax <- max(subset(dina_data, iso == data$iso[1] &
            income_type_short == data$income_type_short[1])$year)
    if (data$year[1] > yearmax - 10) {
        return(NULL)
    }
    if (data$year[1] < 1994 && data$iso[1] == "FR") {
        return(NULL)
    }

    # Calculate short tabulation
    p_in <- c(0, 0.1, 0.5, 0.9, 0.99, 1)
    p_out <- c(0.30, 0.75, 0.95)

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
    average <- data$average[1]
    short_tab <- short_tab[-1, ]
    short_tab$m <- average*short_tab$topshare

    # Generalized Pareto interpolation
    dist <- tabulation_fit(short_tab$p, short_tab$threshold, average, topshare=short_tab$topshare)

    phi_actual <- sapply(p_out, function(p) data[data$p == round(1e5*p), "phi"])
    dphi_actual <- sapply(p_out, function(p) data[data$p == round(1e5*p), "dphi"])

    return(data.frame(
        p        = p_out,
        err_phi  = abs(phi(dist, -log(1 - p_out)) - phi_actual),
        err_dphi = abs(deriv_phi(dist, -log(1 - p_out)) - dphi_actual)
    ))
})

# Get the maximum over all years
err_bound_actual <- ddply(err_actual, c("iso", "country", "p"), function(data) {
    return(data.frame(p=data$p[1], err_phi=max(data$err_phi), err_dphi=max(data$err_dphi)))
})

# Create the LaTeX table
filename <- "output/tables/out-of-sample-error-bounds.tex"

sink(filename)
cat("\\begin{tabular}{cccccccc}\\toprule\n")
cat("& & ")
cat("\\multicolumn{2}{c}{
\\begin{tabular}{@{}c}
maximum absolute \\\\
error on $\\varphi$
\\end{tabular}
} & ")
cat("\\multicolumn{2}{c}{
\\begin{tabular}{@{}c}
maximum absolute \\\\
error on $\\varphi'$
\\end{tabular}
} & ")
cat("\\multicolumn{2}{c}{
\\begin{tabular}{@{}c}
maximum relative \\\\
error on top shares
\\end{tabular}
} \\\\ \\midrule\n")
cat("& & actual & bound & actual & bound & actual & bound \\\\ \\midrule \n")

err_us_actual <- subset(err_bound_actual, iso == "US")
err_us_theory <- subset(err_bound_theory, iso == "US")
for (i in 1:3) {
    if (i == 1) {
        cat("\\multirow{3}{*}{\\begin{tabular}[c]{@{}c@{}}United States\\\\(1962--2004)\\end{tabular}} & ")
    } else {
        cat("& ")
    }
    cat(paste0("$p=", 100*err_us_actual$p[i],"\\%$ & "))

    cat(sprintf("$%.4f$ & ", err_us_actual$err_phi[i]))
    cat(sprintf("$%.4f$ & ", err_us_theory$err_phi[i]))

    cat(sprintf("$%.4f$ & ", err_us_actual$err_dphi[i]))
    cat(sprintf("$%.4f$ & ", err_us_theory$err_dphi[i]))

    cat(sprintf("$%.2f\\%%$ & ", 100*err_us_actual$err_phi[i]))
    cat(sprintf("$%.2f\\%%$ \\\\ \n", 100*err_us_theory$err_phi[i]))
}

cat("\\midrule \n")

err_fr_actual <- subset(err_bound_actual, iso == "FR")
err_fr_theory <- subset(err_bound_theory, iso == "FR")
for (i in 1:3) {
    if (i == 1) {
        cat("\\multirow{3}{*}{\\begin{tabular}[c]{@{}c@{}}France\\\\(1994--2002)\\end{tabular}} & ")
    } else {
        cat("& ")
    }
    cat(paste0("$p=", 100*err_fr_actual$p[i],"\\%$ & "))

    cat(sprintf("$%.4f$ & ", err_fr_actual$err_phi[i]))
    cat(sprintf("$%.4f$ & ", err_fr_theory$err_phi[i]))

    cat(sprintf("$%.4f$ & ", err_fr_actual$err_dphi[i]))
    cat(sprintf("$%.4f$ & ", err_fr_theory$err_dphi[i]))

    cat(sprintf("$%.2f\\%%$ & ", 100*err_fr_actual$err_phi[i]))
    cat(sprintf("$%.2f\\%%$ \\\\ \n", 100*err_fr_theory$err_phi[i]))
}

cat("\\bottomrule \n\\end{tabular}\n")
sink()

