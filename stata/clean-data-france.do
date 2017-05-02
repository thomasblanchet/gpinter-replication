// -------------------------------------------------------------------------- //
// Put France DINA in the right format for the R programs
// -------------------------------------------------------------------------- //

// -------------------------------------------------------------------------- //
// Equal-split fiscal income
// -------------------------------------------------------------------------- //

local excel_sheets TC4 TC13

tempfile data

forvalues j = 1/2 {
	local sheet_name: word `j' of `excel_sheets'
	
	import excel "data-raw/GGP2017DINAAppendixC.xlsx", clear sheet("`sheet_name'") cellrange(A14:EF140)

	// Rename columns
	rename A p
	// List existing columns
	ds p, not
	local vars = r(varlist)
	local n: list sizeof vars
	local n = `n'/3
	forvalues i = 1/`n' {
		local year = 1969 + `i'

		local i1 = 3*(`i' - 1) + 1
		local i2 = 3*(`i' - 1) + 2
		local i3 = 3*(`i' - 1) + 3
		
		local v1: word `i1' of `vars'
		rename `v1' t`year'
		
		local v2: word `i2' of `vars'
		rename `v2' a`year'
		
		local v3: word `i3' of `vars'
		rename `v3' b`year'
	}

	// Multiply percentiles by 1e5 and store them as integer to avoid rounding issues
	replace p = round(p*1000)

	// Put years as a variable
	drop b*
	reshape long a t, i(p) j(year)

	// Get average income
	preserve
	keep if p == 0
	keep year a
	rename a average
	tempfile average
	save "`average'"
	restore

	// Add the above information on averages to the data
	merge n:1 year using "`average'", nogenerate assert(match)

	// Get top share
	rename a topaverage
	rename t threshold
	generate double topshare = topaverage*(1 - p/1e5)/average

	// Get bracket shares
	sort year p
	by year: generate double share = topshare - cond(missing(topshare[_n + 1]), 0, topshare[_n + 1])

	// Other quantities
	generate double invpareto = topaverage/threshold
	generate double dphi = threshold/topaverage
	generate double phi = -log(average*topshare)

	// Save results
	if (`j' == 1) {
		generate income_type = "fiscal income"
		generate income_type_short = "fiscal"
	}
	else if (`j' == 2) {
		generate income_type = "pre-tax national income"
		generate income_type_short = "pretax"
	}
	if (`j' != 1) {
		append using "`data'"
	}
	save "`data'", replace
}

generate country = "France"
generate iso = "FR"
generate population = .
generate population_type = "equal-split adults"

sort year income_type p
order country iso population_type income_type income_type_short year ///
	population average p threshold share topshare topaverage invpareto dphi phi

save "data-clean/data-fr.dta", replace
