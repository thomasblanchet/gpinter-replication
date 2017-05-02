# ---------------------------------------------------------------------------- #
# “Generalized Pareto Curves: Theory and Applications”, 2017
# Thomas Blanchet, Juliette Fournier, Thomas Piketty
# ---------------------------------------------------------------------------- #
# Main R file to reproduce the figures and tables of the article.
# ---------------------------------------------------------------------------- #

# Install the gpinter package if necessary: this package encapsulates the
# main functions required to use generalized Pareto interpolation yourself.
# It is available separately from GitHub.
#library(devtools)
#install_github("thomasblanchet/gpinter")
library(gpinter)

# Other required packages
library(plyr)
library(ggplot2)
library(scales)
library(reshape2)
library(gridExtra)
library(FMStable)
library(pbapply)

# This is required to get the same font as the rest of the LaTeX document in graphs
library(fontcm)
library(extrafont)
loadfonts()

# Graphical parameters for the graphs
# For the article
plot_font       <- "CM Roman"
plot_bg         <- "#FFFFFF"
plot_text_color <- "#000000"

# For the presentation
#plot_font       <- "CM Sans"
#plot_bg         <- "#FBFBFB"
#plot_text_color <- "#23373B"

# Function to extract legends in ggplot
g_legend <- function(gplot) {
    tmp <- ggplot_gtable(ggplot_build(gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
}

# Function to properly display the scientific notation in plots
fancy_scientific <- function(x) {
    # Turn in to character string in scientific notation
    l <- format(x, scientific = TRUE)
    # Quote the part before the exponent to keep all the digits
    l <- gsub("^(.*)e", "'\\1'e", l)
    # Turn the 'e+' into plotmath format
    l <- gsub("e", "%*%10^", l)
    # Return this as an expression
    return(ifelse(x == 0, parse(text="0"), parse(text=l)))
}

# Set the project root directory here
setwd("~/GitHub/gpinter-replication")

# Load the main data file with all the relevant DINA data from France and the
# United States that was created using the Stata programs
dina_data <- read.csv("data-clean/data-us-fr.csv", sep=";", stringsAsFactors=FALSE)

# Load another file with labor and capital income for the US only
dina_data_us <- read.csv("data-clean/data-us-labor-capital.csv", sep=";", stringsAsFactors=FALSE)

# ---------------------------------------------------------------------------- #
# Plot empirical Pareto curves for the distribution of incomes
# ---------------------------------------------------------------------------- #
source("R/plot-empirical-pareto-curves.R")

# ---------------------------------------------------------------------------- #
# Plot Pareto curves for the generalized Pareto distribution
# ---------------------------------------------------------------------------- #
source("R/plot-generalized-pareto-dist.R")

# ---------------------------------------------------------------------------- #
# Compare different interpolation methods
# ---------------------------------------------------------------------------- #

# For testing the interpolation method, we only keep data after 1993 in France
# (data before was calculated with much smaller sample sizes, or tabulations)
dina_data <- subset(dina_data, iso != "FR" | (year > 1993 & year <= 2012))

# Functions for different interpolation methods
source("R/interpolation-methods.R")

# Compare interpolation methods
source("R/compare-interpolation.R")

# Compare extrapolation methods
source("R/compare-extrapolation.R")

# Compare the Pareto curves with the different interpolation methods
source("R/compare-pareto-curves.R")

# ---------------------------------------------------------------------------- #
# Misspecification error
# ---------------------------------------------------------------------------- #

# Plot the error bounds as a multiple of phi'''
source("R/plot-error-bound.R")

# Function for local polynomial fitting on a function and its derivative
source("R/local-polynomial-fitting.R")

# Calculate and plot the estimates of phi'''
source("R/plot-deriv3-phi.R")

# Compare misspecification calculated using the Peano kernel theorem with the
# actual error
source("R/compare-misspecification-error.R")

# Estimate out-of-sample error bounds
source("R/out-of-sample-error-bounds.R")

# Estimate optimal position of brackets
source("R/optimal-position-brackets.R")

# Compare precision with (simulated) survey data
source("R/compare-survey.R")

# ---------------------------------------------------------------------------- #
# Sampling error
# ---------------------------------------------------------------------------- #

# Finite variance case
source("R/sampling-error-finite-variance.R")

# Infinite variance case
source("R/sampling-error-infinite-variance.R")

# ---------------------------------------------------------------------------- #
# Dynamic models of inequality
# ---------------------------------------------------------------------------- #

source("R/pareto-curve-dynamic-model.R")

