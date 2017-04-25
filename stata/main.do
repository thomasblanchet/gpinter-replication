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

// Combine both countries into a single output for fiscal and national income
use "data-clean/data-us.dta", clear
drop if inlist(income_type_short, "capital", "labor")
append using "data-clean/data-fr.dta"

export delimited "data-clean/data-us-fr.csv", replace delimiter(";")

// Create another file with labor and capital income data for the US only
use "data-clean/data-us.dta", clear
keep if inlist(income_type_short, "capital", "labor")
export delimited "data-clean/data-us-labor-capital.csv", replace delimiter(";")

