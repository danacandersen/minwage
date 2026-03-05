/**********************************************************************
  00_DATA_PREP.DO
  Planter Minimum Wage Project

  PURPOSE:
    Import planter_mw.csv, clean and transform variables, save
    planter_mw_prepped.dta for use by all downstream do-files.

  RUN ONCE (or whenever the raw CSV changes).

  SECTIONS:
    (0) Packages
    (1A) Import, types, IDs, day-level aggregates, outcomes, weather
    (1B) Safe string copies (avoid strL issues in collapse/reshape)
    (1C) Year-specific June 1 event time
    (1D) Contract lifecycle + season ramp + cross-season experience
    (SAVE) Write planter_mw_prepped.dta

  OUTPUT:
    logs/00_data_prep_YYYYMMDD.log
    planter_mw_prepped.dta

  REQUIREMENTS:
    reghdfe, ftools, coefplot  (installed in Section 0)
    ~/Downloads/planter_mw.csv
**********************************************************************/

clear all
set more off
version 19


/**********************************************************************
  LOG: open dated log in logs/ subdirectory
**********************************************************************/

* Create logs/ directory if it doesn't exist
cap mkdir "logs"

local logdate : display %tdCYND daily("$S_DATE","DMY")
local logfile "logs/00_data_prep_`logdate'.log"
cap log close _all
log using "`logfile'", replace text
di as text "[00] Log opened: `logfile'"
di as text "[00] Started: $S_DATE $S_TIME"


/**********************************************************************
  (0) PACKAGES — install if missing
**********************************************************************/

di as text "------------------------------------------------------------"
di as text "[0] Checking / installing packages..."
di as text "------------------------------------------------------------"

capture which reghdfe
if _rc {
    di as text "[0] Installing ftools + reghdfe..."
    ssc install ftools,   replace
    ssc install reghdfe,  replace
}
else {
    di as text "[0] reghdfe already installed."
}

capture which coefplot
if _rc {
    di as text "[0] Installing coefplot..."
    ssc install coefplot, replace
}
else {
    di as text "[0] coefplot already installed."
}


/**********************************************************************
  (1A) IMPORT + TYPES + IDs + AGGREGATES + OUTCOMES + WEATHER
**********************************************************************/

di as text "------------------------------------------------------------"
di as text "[1A] Importing planter_mw.csv..."
di as text "------------------------------------------------------------"

* --- Load raw data ---
import delimited using "~/Downloads/planter_mw.csv", clear varnames(1)

* Confirm required variables exist before proceeding
di as text "[1A] Confirming required variables..."
confirm variable date year province contract name piece prod top_up tot_hrs
confirm variable temp total_precip wind_spd experience mw
di as text "[1A] All required variables present."

di as text "[1A] N obs raw = " _N


* --- Date conversion to Stata daily format ---
* If date is stored as a string (e.g. "2018-06-15"), convert to numeric.
capture confirm numeric variable date
if _rc {
    di as text "[1A] date is string — converting with daily(date,YMD)..."
    gen double date_num = daily(date, "YMD")
    format date_num %td
    drop date
    rename date_num date
}
format date %td
di as text "[1A] date range: " %td r(min) " to " %td r(max)


* --- Destring core numeric variables if stored as strings ---
* Handles NA / . / null / NULL values gracefully.
foreach var in temp total_precip wind_spd experience tot_hrs piece prod top_up mw {
    capture confirm numeric variable `var'
    if _rc {
        di as text "[1A] destinging `var'..."
        destring `var', replace ignore("NA" "." "" "null" "NULL")
    }
}

* --- Replace sentinel / implausible weather values with missing ---
* The raw data contains large numeric sentinel codes (e.g. ~1e+14) for
* missing weather readings that were not caught by destring (they are
* already numeric). Plausible ranges: temp -50 to 50°C, precip 0-2000mm,
* wind 0-300 km/h. Anything far outside these is a data error.
* We use conservative thresholds well above any plausible value.
di as text "[1A] Replacing implausible weather sentinels with missing..."

quietly count if !missing(temp) & (temp > 500 | temp < -200)
di as text "  temp sentinels replaced: " r(N)
replace temp = . if temp > 500 | temp < -200

quietly count if !missing(total_precip) & total_precip > 5000
di as text "  total_precip sentinels replaced: " r(N)
replace total_precip = . if total_precip > 5000

quietly count if !missing(wind_spd) & wind_spd > 500
di as text "  wind_spd sentinels replaced: " r(N)
replace wind_spd = . if wind_spd > 500

di as text "[1A] Weather missingness after sentinel removal:"
foreach var in temp total_precip wind_spd {
    quietly count if missing(`var')
    di as text "  `var' missing: " r(N) " (" %5.1f (r(N)/_N*100) "%)"
}


* --- Worker / contract IDs (encoded once) ---
* encode assigns integer codes; wc_id is the worker×contract cell.
cap drop worker_id contract_id wc_id
encode name,     gen(worker_id)
encode contract, gen(contract_id)
egen   wc_id = group(worker_id contract_id)

di as text "[1A] N workers:   " r(max)   // after encode
quietly count
di as text "[1A] N obs:       " r(N)


* --- Day-level aggregates (handles multi-contract days) ---
* tot_hrs_day = total hours worked by worker on a given date across all contracts.
* n_contracts_day flags days where a worker plants in more than one contract.
cap drop tot_hrs_day n_contracts_day multi_contract_day long_day extreme_day
bysort worker_id date: egen double tot_hrs_day     = total(tot_hrs)
bysort worker_id date: egen int    n_contracts_day = count(contract_id)
gen byte multi_contract_day = (n_contracts_day > 1)
gen byte long_day           = (tot_hrs_day > 14)    // flagging unusually long days
gen byte extreme_day        = (tot_hrs_day > 18)


* --- Outcome variables ---
* ln_incentive: log of total piece-rate earnings per hour (piece × productivity).
* lprod: log trees/hour.
* ln_piece: log piece rate.
* All set to missing when inputs are non-positive (log undefined).
cap drop ln_incentive lprod ln_piece
gen double ln_incentive = ln(piece * prod) if piece > 0 & prod > 0
gen double lprod        = ln(prod)          if prod > 0
gen double ln_piece     = ln(piece)         if piece > 0

di as text "[1A] Outcome non-missing counts:"
count if !missing(ln_incentive)
count if !missing(lprod)
count if !missing(ln_piece)


* --- Weather transforms ---
* precip0 flags zero-precip days (separate intercept shift).
* precip_pos is precip conditional on > 0 (avoids square-root compression for zeros).
* Squared terms allow flexible non-linear weather-productivity relationships.
cap drop precip0 precip_pos temp2 wind2 precip2
gen byte   precip0    = (total_precip == 0)                                  if !missing(total_precip)
gen double precip_pos = cond(total_precip > 0 & !missing(total_precip), total_precip, 0)
gen double temp2      = temp^2          if !missing(temp)
gen double wind2      = wind_spd^2      if !missing(wind_spd)
gen double precip2    = precip_pos^2    if !missing(precip_pos)

* Experience squared for quadratic experience control
cap drop experience2
gen double experience2 = experience^2   if !missing(experience)


/**********************************************************************
  (1B) SAFE STRING COPIES
  Avoid strL storage type which causes errors in collapse/reshape.
  province_s / contract_s / name_s used throughout analysis.
**********************************************************************/

di as text "------------------------------------------------------------"
di as text "[1B] Creating safe string copies..."
di as text "------------------------------------------------------------"

cap drop province_s contract_s name_s
gen str2   province_s = province
gen str40  contract_s = contract
gen str40  name_s     = name

* Quick check — province values present
tab province_s, miss


/**********************************************************************
  (1C) YEAR-SPECIFIC JUNE 1 EVENT TIME
  Computes event time relative to June 1 of each observation's year.
  This lets the same variables work for all years (2016, 2017, 2018, 2019…).

  tdate_y      = June 1 of that year (Stata daily)
  post_y       = 1 if date >= June 1 of that year
  event_day_y  = calendar days since June 1
  event_week_y = event weeks (floor division; -1 = week before cutoff, 0 = cutoff week)
**********************************************************************/

di as text "------------------------------------------------------------"
di as text "[1C] Building year-specific June 1 event time..."
di as text "------------------------------------------------------------"

cap drop tdate_y post_y event_day_y event_week_y
gen double tdate_y      = mdy(6, 1, year)
format tdate_y %td
gen byte   post_y       = (date >= tdate_y)       if !missing(tdate_y)
gen int    event_day_y  = date - tdate_y           if !missing(tdate_y)
gen int    event_week_y = floor(event_day_y / 7)   if !missing(event_day_y)

di as text "[1C] event_week_y distribution:"
summ event_week_y, detail


/**********************************************************************
  (1D) CONTRACT LIFECYCLE + SEASON RAMP + CROSS-SEASON EXPERIENCE

  contract_day:        days since first day of this contract×year cell.
  contract_day2:       squared term for quadratic ramp.
  season_days_worked:  cumulative distinct work days within worker×year
                       (proxy for within-season fatigue / ramp-up).
  season_days_worked2: squared term.
  cum_seasons:         number of prior seasons observed (0 in first season).
                       Proxy for cross-season learning / accumulated skill.
**********************************************************************/

di as text "------------------------------------------------------------"
di as text "[1D] Building contract lifecycle + season ramp + cum_seasons..."
di as text "------------------------------------------------------------"

* Contract lifecycle within year
cap drop contract_start contract_day contract_day2
bysort contract_id year: egen double contract_start = min(date)
gen int contract_day  = date - contract_start
gen int contract_day2 = contract_day^2

* Within-season ramp: count distinct calendar days worked by worker within year.
* first_day flags the first row for each date within worker×year (sorted by date).
cap drop first_day season_days_worked season_days_worked2
bysort worker_id year (date): gen byte first_day          = (date != date[_n-1])
bysort worker_id year (date): gen int  season_days_worked = sum(first_day)
replace first_day = 1 if _n==1    // ensure first observation is always flagged
gen int season_days_worked2 = season_days_worked^2

* Cross-season experience: count seasons (years) observed up to and including current year.
* cum_seasons = 0 in the worker's first observed season.
cap drop yrflag cum_seasons
bysort worker_id year: gen byte yrflag    = (_n == 1)
bysort worker_id (year): gen int cum_seasons = sum(yrflag) - 1
drop yrflag

di as text "[1D] season_days_worked stats:"
summ season_days_worked
di as text "[1D] cum_seasons stats:"
summ cum_seasons


/**********************************************************************
  SAVE: write prepped dataset for downstream do-files
**********************************************************************/

di as text "------------------------------------------------------------"
di as text "[SAVE] Saving planter_mw_prepped.dta..."
di as text "------------------------------------------------------------"

save "planter_mw_prepped.dta", replace

di as text "[SAVE] Done. N=" _N " obs saved."
di as text ""
di as text "Variables in saved dataset:"
describe


/**********************************************************************
  LOG: close
**********************************************************************/

di as text ""
di as text "[00] Finished: $S_DATE $S_TIME"
log close
