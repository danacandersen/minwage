/**********************************************************************
  01_DATA_CHECKS.DO
  Planter Minimum Wage Project

  PURPOSE:
    Validate planter_mw_prepped.dta and run diagnostic checks.
    Covers coverage, missingness, outcome distributions, top-up
    validity, weather, and experience.

    KNOWN ISSUE: top_up may not align perfectly with minimum-wage
    binding (piece * prod * tot_hrs < mw * tot_hrs). This file
    quantifies the discrepancy.

  SECTIONS:
    (1) Coverage: obs/worker/contract counts by province×year
    (2) Missing value diagnostics
    (3) Outcome distributions by province and pre/post
    (4) Top-up validity check
    (5) Weather diagnostics
    (6) Experience distribution

  OUTPUT:
    logs/01_data_checks_YYYYMMDD.log

  REQUIRES: 00_data_prep.do to have been run first.
**********************************************************************/

clear all
set more off
version 19


/**********************************************************************
  LOG
**********************************************************************/

cap mkdir "logs"

local logdate : display %tdCYND daily("$S_DATE","DMY")
local logfile "logs/01_data_checks_`logdate'.log"
cap log close _all
log using "`logfile'", replace text
di as text "[01] Log opened: `logfile'"
di as text "[01] Started: $S_DATE $S_TIME"


/**********************************************************************
  LOAD
**********************************************************************/

use "planter_mw_prepped.dta", clear
di as text "[01] Loaded planter_mw_prepped.dta — N=" _N


/**********************************************************************
  (1) COVERAGE CHECKS
**********************************************************************/

di as text "============================================================"
di as text "[1] COVERAGE: obs/worker/contract by province and year"
di as text "============================================================"

* Total obs by province
di as text "--- Obs by province ---"
tab province_s, miss

* Obs by year
di as text "--- Obs by year ---"
tab year, miss

* Obs by province × year
di as text "--- Obs by province x year ---"
tab province_s year, miss

* Worker count by province × year
di as text "--- Unique workers by province x year ---"
preserve
    bysort province_s year worker_id: keep if _n==1
    tab province_s year
restore

* Contract count by province × year
di as text "--- Unique contracts by province x year ---"
preserve
    bysort province_s year contract_id: keep if _n==1
    tab province_s year
restore

* Date range by province × year
di as text "--- Date range by province x year ---"
table province_s year, statistic(min date) statistic(max date)


/**********************************************************************
  (2) MISSING VALUE DIAGNOSTICS
**********************************************************************/

di as text "============================================================"
di as text "[2] MISSINGNESS: key variables"
di as text "============================================================"

foreach var in ln_incentive lprod ln_piece tot_hrs_day ///
               piece prod top_up mw ///   // top_up is 0/1 dummy
               temp total_precip wind_spd ///
               experience cum_seasons ///
               post_y event_week_y ///
               contract_day season_days_worked {

    quietly count if missing(`var')
    local nmiss = r(N)
    quietly count
    local ntot  = r(N)
    local pct   = 100 * `nmiss' / `ntot'
    di as text "  `var':   " `nmiss' " missing of " `ntot' " (" %5.1f `pct' "%)"
}


/**********************************************************************
  (3) OUTCOME DISTRIBUTIONS
**********************************************************************/

di as text "============================================================"
di as text "[3] OUTCOME DISTRIBUTIONS"
di as text "============================================================"

* Overall
di as text "--- Overall ---"
* top_up is a 0/1 dummy; summ shows share of topped-up obs
summ ln_incentive lprod ln_piece tot_hrs_day piece prod top_up, detail

* By province
di as text "--- By province ---"
foreach prov in BC AB ON {
    di as text "  Province = `prov':"
    summ ln_incentive lprod ln_piece tot_hrs_day if province_s=="`prov'"
}

* By pre/post June 1 (using year-specific post_y)
di as text "--- By post_y (pre/post June 1, year-specific) ---"
foreach pp in 0 1 {
    local lab = cond(`pp'==0, "PRE", "POST")
    di as text "  `lab' (post_y=`pp'):"
    summ ln_incentive lprod ln_piece tot_hrs_day if post_y==`pp'
}

* By province × pre/post
di as text "--- By province x pre/post ---"
table province_s post_y, statistic(mean ln_incentive) statistic(mean lprod) ///
    statistic(mean ln_piece) nformat(%9.3f)


/**********************************************************************
  (4) TOP-UP VALIDITY CHECK

  top_up is a 0/1 dummy for whether the worker received a top-up
  on a given day. The actual dollar top-up amount is not in the data
  but can be constructed as:

      topup_amount = max(0, (mw - piece*prod) * tot_hrs)

  where:
      piece*prod  = hourly piece-rate earnings ($/hour)
      mw          = minimum wage ($/hour)
      tot_hrs     = hours worked that day

  The MW is binding (top_up should = 1) whenever piece*prod < mw.
  This block checks alignment between the top_up dummy and the
  constructed binding indicator, and quantifies the imputed top-up
  dollar amounts.
**********************************************************************/

di as text "============================================================"
di as text "[4] TOP-UP VALIDITY CHECK"
di as text "============================================================"

* Temp variables — dropped at end of section
cap drop earn_hourly mw_binding topup_amount topup_err_above topup_err_zero

* Hourly piece-rate earnings: piece ($/tree) × prod (trees/hour) = $/hour
gen double earn_hourly = piece * prod if !missing(piece) & !missing(prod)

* MW binding: hourly piece earnings fall short of the MW rate
* (comparison is hourly, so no need to multiply by tot_hrs)
gen byte mw_binding = (earn_hourly < mw) ///
    if !missing(earn_hourly) & !missing(mw)

* Imputed top-up dollar amount: daily shortfall below MW floor
* = max(0, (mw - piece*prod) * tot_hrs)
gen double topup_amount = max(0, (mw - earn_hourly) * tot_hrs) ///
    if !missing(earn_hourly) & !missing(mw) & !missing(tot_hrs)

di as text "--- top_up dummy distribution ---"
tab top_up, miss

di as text ""
di as text "--- MW binding flag: piece*prod < mw ---"
tab mw_binding, miss

di as text ""
di as text "--- Cross-tab: top_up dummy x mw_binding ---"
di as text "    (rows = top_up dummy, cols = piece*prod < mw)"
di as text "    Perfect alignment: only cells (0,0) and (1,1) should be populated."
tab top_up mw_binding, row col miss

* Error type 1: worker received top-up but piece earnings are ABOVE mw
gen byte topup_err_above = (top_up==1) & (mw_binding==0) ///
    if !missing(top_up) & !missing(mw_binding)

di as text ""
di as text "--- Error type 1: top_up=1 but piece*prod >= mw ---"
count if topup_err_above==1
di as text "    (workers flagged as topped-up but appear above MW floor)"

di as text "    Hourly earnings for these obs:"
summ earn_hourly if topup_err_above==1
di as text "    MW for these obs:"
summ mw if topup_err_above==1

* Error type 2: no top-up recorded but piece earnings are BELOW mw
gen byte topup_err_zero = (top_up==0) & (mw_binding==1) ///
    if !missing(top_up) & !missing(mw_binding)

di as text ""
di as text "--- Error type 2: top_up=0 but piece*prod < mw ---"
count if topup_err_zero==1
di as text "    (workers appear below MW floor but no top-up recorded)"

di as text "    Hourly earnings for these obs:"
summ earn_hourly if topup_err_zero==1
di as text "    MW for these obs:"
summ mw if topup_err_zero==1
di as text "    Imputed shortfall (mw - piece*prod) for these obs:"
gen double shortfall = mw - earn_hourly if topup_err_zero==1
summ shortfall if topup_err_zero==1, detail
drop shortfall

* Imputed top-up amounts for obs where binding (regardless of dummy)
di as text ""
di as text "--- Imputed top-up dollar amount when MW is binding (mw_binding==1) ---"
summ topup_amount if mw_binding==1, detail

di as text ""
di as text "--- Imputed top-up by province x year (mean, when binding) ---"
table province_s year if mw_binding==1, ///
    statistic(mean topup_amount) statistic(mean earn_hourly) statistic(mean mw) ///
    nformat(%8.3f)

di as text ""
di as text "--- Share of obs where MW is binding, by province x year ---"
table province_s year, statistic(mean mw_binding) nformat(%6.3f)

* Clean up
drop earn_hourly mw_binding topup_amount topup_err_above topup_err_zero


/**********************************************************************
  (5) WEATHER DIAGNOSTICS
**********************************************************************/

di as text "============================================================"
di as text "[5] WEATHER DIAGNOSTICS by province x year"
di as text "============================================================"

* Distribution overall
di as text "--- temp, total_precip, wind_spd overall ---"
summ temp total_precip wind_spd, detail

* By province × year
di as text "--- Mean weather by province x year ---"
table province_s year, ///
    statistic(mean temp) statistic(mean total_precip) statistic(mean wind_spd) ///
    nformat(%8.2f)

* Zero-precip share by province × year
di as text "--- Zero-precip share by province x year ---"
table province_s year if !missing(total_precip), ///
    statistic(mean precip0) nformat(%6.3f)


/**********************************************************************
  (6) EXPERIENCE DISTRIBUTION
**********************************************************************/

di as text "============================================================"
di as text "[6] EXPERIENCE DISTRIBUTION"
di as text "============================================================"

di as text "--- experience by province ---"
table province_s, statistic(mean experience) statistic(sd experience) ///
    statistic(min experience) statistic(max experience) nformat(%8.2f)

di as text "--- cum_seasons by province ---"
tab cum_seasons province_s, miss

di as text "--- season_days_worked distribution ---"
summ season_days_worked, detail


/**********************************************************************
  LOG: close
**********************************************************************/

di as text ""
di as text "[01] Finished: $S_DATE $S_TIME"
log close
