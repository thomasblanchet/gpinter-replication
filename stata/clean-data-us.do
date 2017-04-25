// -------------------------------------------------------------------------- //
// Put United States DINA in the right format for the R programs
// -------------------------------------------------------------------------- //

import delimited "data-raw/US_income.csv", clear delimiter(";")

// Country code & country name
rename alpha2 iso
generate country = "United States"

// Parse percentiles (multiply them by 1e5 and store them as integer to avoid
// rounding issues)
generate long p = round(1000*real(substr(p2, 2, .)))
replace p = -1 if (p2 == "pall")

// Keep 4 types of income: fiscal, pre-tax total, pre-tax labor and pre-tax capital
keep iso country year p npopul992j ?fiinc992j ?ptinc992j ?ptlin992j ?ptkin992j

// Reshape to long format (put income type as a variable)
reshape long n a s t, i(year p) j(income_type) string

// Get the population size in each year
preserve
keep if income_type == "popul992j"
collapse (firstnm) n, by(iso year)
rename n population
tempfile population
save "`population'"
restore

// Get the average income in each year
drop if income_type == "popul992j"
preserve
sort year income_type p
by year income_type: replace a = a*(cond(missing(p[_n + 1]), 1e5, p[_n + 1]) - p)/1e5
collapse (sum) a, by(iso year income_type)
rename a average
tempfile average
save "`average'"
restore

// Add previously calcualted population size and average income
drop a n
drop if p == -1
merge n:1 iso year using "`population'", nogenerate assert(match)
merge n:1 iso year income_type using "`average'", nogenerate assert(match) 

// Put income type in friendlier format
replace income_type = "fiscal income"           if (income_type == "fiinc992j")
replace income_type = "pre-tax national income" if (income_type == "ptinc992j")
replace income_type = "pre-tax capital income"  if (income_type == "ptkin992j")
replace income_type = "pre-tax labor income"    if (income_type == "ptlin992j")

generate income_type_short = ""
replace income_type_short = "fiscal"  if (income_type == "fiscal income")
replace income_type_short = "pretax"  if (income_type == "pre-tax national income")
replace income_type_short = "capital" if (income_type == "pre-tax capital income")
replace income_type_short = "labor"   if (income_type == "pre-tax labor income")

generate population_type = "equal-split adults"

// Compute top shares, Pareto coefficients, etc.
rename s share
rename t threshold

gsort year income_type -p
by year income_type: generate topshare = sum(share)
sort year income_type p
// Normalize shares to get rid of a small discrepancy that prevent them from
// summing to one
recast double share topshare
by year income_type: replace share = share/topshare[1]
by year income_type: replace topshare = topshare/topshare[1]

generate double topaverage = average*topshare/(1 - p/1e5)

generate double invpareto = topaverage/threshold
generate double dphi = threshold/topaverage
generate double phi = -log(topshare*average)

// Save results
order country iso population_type income_type income_type_short year ///
	population average p threshold share topshare topaverage invpareto dphi phi
save "data-clean/data-us.dta", replace

