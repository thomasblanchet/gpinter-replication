# Generalized Pareto Curves: Theory and Applications

Codes to replicate our article on generalized Pareto interpolation “Generalized Pareto Curves: Theory and Applications” by Thomas Blanchet, Juliette Fournier and Thomas Piketty.

These codes replicate the tables and figures of the article. If you are interested in using generalized Pareto interpolation for your own work, please see [the `gpinter` package](https://github.com/thomasblanchet/gpinter) or our [online application](http://wid.world/gpinter/).

## Replication instructions

First run the following Stata DO files that import and clean up data from the World Wealth and Income Database.
- `stata/clean-france.do`
- `stata/clean-data.do`

Then run the file `R/main.R`.
