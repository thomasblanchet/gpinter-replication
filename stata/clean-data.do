// -------------------------------------------------------------------------- //
// Put the data from the World Wealth and Income Database in suitable format
// for use with the R code
// -------------------------------------------------------------------------- //

// Need the ADO 'kountry' to deal with country ISO codes
ssc install kountry

// Set the project root directory here
cd "~/GitHub/gpinter-replication"

// Decompress the raw data
cd "data-wid"
unzipfile "world-wealth-and-income-database.zip", replace
cd ".."

// Import the raw data
import delimited "data-wid/world-wealth-and-income-database.csv", clear

// -------------------------------------------------------------------------- //
// Basic cleaning up and formatting
// -------------------------------------------------------------------------- //

// Keep only:
//   1) pre-tax national income (ptinc)
//   2) pre-tax labor income (pllin)
//   3) pre-tax capital income (pkkin)
//   4) fiscal income (fiinc)
//   5) household wealth (hweal)
//   6) associated populations
keep iso year p ?ptinc* ?pllin* ?pkkin* ?fiinc* ?hweal* npopul999? npopul992?

// Remove m- variables (macro aggregates) and o- variables (average above,
// redundant), b- variables (to be recalculated anyway)
keep iso year p a* s* t* n*

// Reshape to put type of income/wealth/population as a variable
reshape long a s t n, i(iso p year) j(widcode) string

// Drop observations with no income/wealth data
egen to_keep = rownonmiss(a s t n)
keep if to_keep
drop to_keep

// Store the data in a temporary file for later use
save "data-temp/wid-db.dta", replace

// -------------------------------------------------------------------------- //
// Overall averages
// -------------------------------------------------------------------------- //

use "data-temp/wid-db.dta", clear

// Keep aveages only
keep iso year p widcode a
drop if missing(a)

// Keep averages for the whole distribution
keep if (p == "p0p100")
drop p

save "data-temp/averages.dta", replace

// -------------------------------------------------------------------------- //
// Top shares
// -------------------------------------------------------------------------- //

use "data-temp/wid-db.dta", clear

// Keep shares only
keep iso year p widcode s
drop if missing(s)

// Keep top share only
keep if regexm(p, "^p(.*)p100$")

// Parse percentile
generate long p_integer = round(1e3*real(regexs(1))) if regexm(p, "^p(.*)p100$")
drop p
rename p_integer p

save "data-temp/top-shares.dta", replace

// -------------------------------------------------------------------------- //
// Thresholds
// -------------------------------------------------------------------------- //

use "data-temp/wid-db.dta", clear

// Keep thresholds only
keep iso year p widcode t
drop if missing(t)

// Parse percentile
generate long p_integer = round(1e3*real(regexs(1))) if regexm(p, "^p(.*)$")
drop p
rename p_integer p

save "data-temp/thresholds.dta", replace

// -------------------------------------------------------------------------- //
// Combine the data calculated above
// -------------------------------------------------------------------------- //

use "data-temp/top-shares.dta", clear
merge 1:1 iso year widcode p using "data-temp/thresholds.dta", nogenerate keep(match)
merge n:1 iso year widcode using "data-temp/averages.dta", nogenerate keep(match)

rename a average
rename t threshold
rename s top_share

generate double top_average = average*top_share/(1 - p/1e5)
generate double b           = top_average/threshold
generate double phi         = -log(top_share*average)
generate double dphi        = threshold/top_average

save "data-temp/combined.dta", replace

// -------------------------------------------------------------------------- //
// Final cleaning up
// -------------------------------------------------------------------------- //

use "data-temp/combined.dta", clear

// Generate full country names
kountry iso, from(iso2c)
rename NAMES_STD country
replace country = "Urban China" if (iso == "CN-UR")
replace country = "Rural China" if (iso == "CN-RU")
replace country = "Zanzibar"    if (iso == "ZZ")

// Explain WID codes
generate var_code = substr(widcode, 1, 5)
generate pop_code = substr(widcode, 9, 1)

generate var_name = ""
replace var_name = "fiscal income"           if (var_code == "fiinc")
replace var_name = "personal wealth"         if (var_code == "hweal")
replace var_name = "capital income"          if (var_code == "pkkin")
replace var_name = "labor income"            if (var_code == "pllin")
replace var_name = "pre-tax national income" if (var_code == "ptinc")

keep if inlist(pop_code, "i", "j", "t")
generate pop_name = ""
replace pop_name = "individuals"             if (pop_code == "i")
replace pop_name = "equal-split individuals" if (pop_code == "j")
replace pop_name = "tax units"               if (pop_code == "t")

order iso country widcode var_code var_name pop_code pop_name year p ///
	average top_share top_average threshold b phi dphi

sort iso country widcode var_code var_name pop_code pop_name year p

// -------------------------------------------------------------------------- //
// Save the results
// -------------------------------------------------------------------------- //

export delimited "data-clean/data-wid.csv", replace

