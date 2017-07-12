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
	
	import excel "data-wid/GGP2017DINAAppendixC.xlsx", clear sheet("`sheet_name'") cellrange(A14:EF140)

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
	rename a top_average
	rename t threshold
	generate double top_share = top_average*(1 - p/1e5)/average

	// Other quantities
	generate double b = top_average/threshold
	generate double dphi = threshold/top_average
	generate double phi = -log(average*top_share)

	// Save results
	if ("`sheet_name'" == "TC4") {
		generate widcode = "fiinc992j"
	}
	else if ("`sheet_name'" == "TC13") {
		generate widcode = "ptinc992j"
	}
	if (`j' != 1) {
		append using "`data'"
	}
	save "`data'", replace
}

generate iso = "FR"

save "data-temp/dina-fr.dta", replace
