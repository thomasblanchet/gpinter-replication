// -------------------------------------------------------------------------- //
// Put DINA data for France and the United States in the right format
// for the R programs.
// -------------------------------------------------------------------------- //

// Set the project root directory here
cd "~/GitHub/gpinter-replication"

// Clean US DINA data
do "stata/clean-data-us.do"

// Clean France DINA data
do "stata/clean-data-france.do"

// Combine both countries into a single output
use "data-clean/data-us.dta", clear
append using "data-clean/data-fr.dta"

export delimited "data-clean/data-both.csv", replace delimiter(";")

