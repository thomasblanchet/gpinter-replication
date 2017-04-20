# ---------------------------------------------------------------------------- #
# Main R file to reproduce the tables and figures in “Generalized Pareto
# Curves: Theory and Application” (2017), Thomas Blanchet, Juliette Fournier
# and Thomas Piketty.
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

# This is required to get the same font as the rest of the LaTeX document in graphs
library(fontcm)
library(extrafont)
loadfonts()

# Function to extract legends in ggplot
g_legend <- function(gplot){
    tmp <- ggplot_gtable(ggplot_build(gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
}

# Set the project root directory here
setwd("~/GitHub/gpinter-replication")

# Load the main data file with all the relevant DINA data from France and the
# United States: that was created using the Stata programs
dina_data <- read.csv("data-clean/data-both.csv", sep=";", stringsAsFactors=FALSE)

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

# Functions for different interpolation methods
source("R/interpolation-methods.R")

# Compare interpolation methods
source("R/compare-interpolation.R")

# Compare extrapolation methods
source("R/compare-extrapolation.R")

# Compare the Pareto curves with the different methods
source("R/compare-pareto-curves.R")



