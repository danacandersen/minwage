/**********************************************************************
  02_MAIN_ANALYSIS.DO
  Planter Minimum Wage Project

  PURPOSE:
    Run the two baseline DiD specifications:
      2A) BC vs AB/ON: cross-province DiD in the treated year
      2B) Within-BC placebo years: treated year vs placebo year(s)
    Both include auto event-window selection (TAU rule), event study,
    pre-trend test, and coefplot export to output/.

  CONTROL PANEL at top — all parameters as globals, never hardcoded.

  SECTIONS:
    2.1) Control panel (EDIT THESE)
    2.2) Build controls ($X, $LIFECYCLE, $ABS, $SAMPLE_OK)
    2A)  BC vs AB/ON DiD + event study + coefplot
    2B)  Within-BC placebo years DiD + event study + coefplot

  OUTPUT:
    logs/02_main_analysis_YYYYMMDD.log
    output/event_2A_*.png
    output/event_2B_*.png

  REQUIRES: 00_data_prep.do to have been run first.
**********************************************************************/

clear all
set more off
version 19


/**********************************************************************
  LOG
**********************************************************************/

cap mkdir "logs"
cap mkdir "output"

local logdate : display %tdCYND daily("$S_DATE","DMY")
local logfile "logs/02_main_analysis_`logdate'.log"
cap log close _all
log using "`logfile'", replace text
di as text "[02] Log opened: `logfile'"
di as text "[02] Started: $S_DATE $S_TIME"


/**********************************************************************
  LOAD
**********************************************************************/

use "planter_mw_prepped.dta", clear
di as text "[02] Loaded planter_mw_prepped.dta — N=" _N


/**********************************************************************
  CLEAN RERUN STATE
  Safe drops of any variables/macros left over from a previous run.
**********************************************************************/

cap estimates clear
cap label drop ew_lbl ew_lbl2

foreach v in ///
    sample_2A treated_BC_2A post_2A event_week_2A event_week_FE ///
    in_esamp_2A bc_obs_2A ctrl_obs_2A bcN_2A ctrlN_2A okweek_2A event_week_s_2A ///
    sample_2B treated_year_2B placebo_year_2B post_y_2B event_week_y_2B event_weeky_FE ///
    in_esamp_2B t_obs_2B p_obs_2B tN_2B pN_2B okweek_2B event_week_sy_2B {
    cap drop `v'
}

foreach m in ///
    Y FE CLUSTVAR ///
    USE_WEATHER USE_EXPER USE_HOURS USE_DAYFLAGS USE_LIFECYCLE USE_TIMEFE ///
    USE_FLOORS MIN_HRS MIN_PROD MIN_EARN USE_MIN_PROD ///
    TAU CORE0 CORE1 MAXW AUTO_WIN_2A AUTO_WIN_2B ///
    WPRE_MAN_2A WPOST_MAN_2A PRE_MIN_MAN_2B POST_MAX_MAN_2B ///
    TY_2A CONTROL_2A TY_2B BCPLACEBO_2B ///
    X ABS SAMPLE_OK LIFECYCLE TIMEABS ///
    WPRE_2A WPOST_2A WINMETHOD_2A ///
    PRE_MIN_2B POST_MAX_2B WINMETHOD_2B ///
    EW_BASE_2A EW_BASE_2B {
    cap macro drop `m'
}


/**********************************************************************
  2.1) CONTROL PANEL — EDIT THESE
**********************************************************************/

di as text "============================================================"
di as text "[2.1] Control panel"
di as text "============================================================"

* ----- Outcome variable -----
* Options: ln_incentive  lprod  ln_piece  tot_hrs_day
global Y ln_incentive

* ----- Fixed effects -----
* worker_contract -> absorb(worker_id contract_id)   [BASELINE]
* wc_only         -> absorb(wc_id)
* worker_only     -> absorb(worker_id)
* contract_only   -> absorb(contract_id)
global FE worker_contract

* ----- Clustering -----
* Options: worker_id  contract_id  wc_id
global CLUSTVAR worker_id

* ----- Control toggles -----
global USE_WEATHER   1    // temp, wind, precip (quadratic)
global USE_EXPER     1    // experience quadratic
global USE_HOURS     1    // tot_hrs_day
global USE_DAYFLAGS  1    // multi_contract_day
global USE_LIFECYCLE 1    // contract ramp + season ramp + ramp×cum_seasons
global USE_TIMEFE    1    // 1 = absorb event-week FE in DiD (not in event study)

* ----- Floor filters -----
global USE_FLOORS   1
global MIN_HRS      2     // minimum hours/day
global MIN_PROD    50     // minimum productivity (trees/hour)
global MIN_EARN    30     // minimum piece-rate earnings (piece*prod)
global USE_MIN_PROD 1     // also apply MIN_PROD filter

* ----- Auto event-window (TAU rule) -----
* TAU: fraction of core-week density required at each event week.
* CORE0/CORE1: weeks defining the "core" reference window.
* MAXW: maximum weeks to search.
* AUTO_WIN_2A/2B: 1 = use TAU rule; 0 = use manual windows.
global TAU    0.25
global CORE0  0
global CORE1  1
global MAXW   30
global AUTO_WIN_2A 1
global AUTO_WIN_2B 1

* Manual windows (used if AUTO fails or is turned off)
global WPRE_MAN_2A      3
global WPOST_MAN_2A     8
global PRE_MIN_MAN_2B  -6
global POST_MAX_MAN_2B 10

* ----- 2A settings -----
global TY_2A      2018
global CONTROL_2A "AB"   // "AB", "ON", or "ABON" (auto-window disabled for ABON)

* Force manual window if ABON (two control provinces complicates TAU reshape)
if "$CONTROL_2A"=="ABON" global AUTO_WIN_2A 0

* ----- 2B settings -----
global TY_2B       2018
global BCPLACEBO_2B "2016 2017"


/**********************************************************************
  2.2) BUILD CONTROLS, LIFECYCLE, FE ABSORB, FLOOR FILTER
**********************************************************************/

di as text "------------------------------------------------------------"
di as text "[2.2] Building globals from control panel..."
di as text "------------------------------------------------------------"

* Check CLUSTVAR exists
capture confirm variable $CLUSTVAR
if _rc {
    di as error "[FATAL] CLUSTVAR=$CLUSTVAR not found. Fix global CLUSTVAR."
    exit 198
}

* Main controls
global X ""
if $USE_HOURS==1    global X "$X tot_hrs_day"
if $USE_DAYFLAGS==1 global X "$X multi_contract_day"
if $USE_EXPER==1    global X "$X c.experience##c.experience"
if $USE_WEATHER==1  global X "$X c.temp##c.temp c.wind_spd##c.wind_spd precip0 c.precip_pos##c.precip_pos"

* Lifecycle/ramp controls
global LIFECYCLE ""
if $USE_LIFECYCLE==1 {
    global LIFECYCLE "$LIFECYCLE c.contract_day##c.contract_day"
    global LIFECYCLE "$LIFECYCLE c.season_days_worked##c.season_days_worked"
    global LIFECYCLE "$LIFECYCLE c.season_days_worked##c.cum_seasons"
}

* Fixed effects absorb string
global ABS ""
if "$FE"=="worker_contract" global ABS "worker_id contract_id"
if "$FE"=="wc_only"         global ABS "wc_id"
if "$FE"=="worker_only"     global ABS "worker_id"
if "$FE"=="contract_only"   global ABS "contract_id"

* Time FE absorb (added only to baseline DiD, not event study)
global TIMEABS ""
if $USE_TIMEFE==1 global TIMEABS "event_week_FE"

* Floor / inclusion filter string
global SAMPLE_OK "1==1"
if $USE_FLOORS==1 {
    global SAMPLE_OK "$SAMPLE_OK & tot_hrs_day >= $MIN_HRS"
    global SAMPLE_OK "$SAMPLE_OK & piece > 0 & prod > 0"
    global SAMPLE_OK "$SAMPLE_OK & (piece*prod) >= $MIN_EARN"
    if $USE_MIN_PROD==1 global SAMPLE_OK "$SAMPLE_OK & prod >= $MIN_PROD"
}

* Debug output
di as text "[2.2] Y:          $Y"
di as text "[2.2] FE:         $FE  ->  ABS: $ABS"
di as text "[2.2] CLUSTVAR:   $CLUSTVAR"
di as text "[2.2] X:          $X"
di as text "[2.2] LIFECYCLE:  $LIFECYCLE"
di as text "[2.2] TIMEABS:    $TIMEABS"
di as text "[2.2] SAMPLE_OK:  $SAMPLE_OK"
di as text "[2.2] TY_2A:      $TY_2A   CONTROL_2A: $CONTROL_2A"
di as text "[2.2] TY_2B:      $TY_2B   BCPLACEBO:  $BCPLACEBO_2B"
di as text "[2.2] TAU:        $TAU   MAXW: $MAXW"


/**********************************************************************
  2A) BC vs CONTROL PROVINCE — treated year only

  Strategy: compare tree planters in BC vs AB (or ON) in the year of
  the MW reform ($TY_2A). Post is defined as date >= June 1 of that year.

  Two regressions:
    (i)  Baseline DiD: 2×2 regression with event-week FE (if USE_TIMEFE).
    (ii) Event study: interactions of event_week × treated, base=week -1.
**********************************************************************/

di as text "============================================================"
di as text "[2A] BC vs $CONTROL_2A ($TY_2A)"
di as text "============================================================"

* --- Sample + event time ---
cap drop sample_2A treated_BC_2A post_2A event_week_2A event_week_FE

gen byte sample_2A = (year == $TY_2A) & inlist(province_s, "BC", "$CONTROL_2A")
if "$CONTROL_2A"=="ABON" ///
    replace sample_2A = (year == $TY_2A) & inlist(province_s, "BC", "AB", "ON")

gen byte treated_BC_2A = .
replace treated_BC_2A = (province_s == "BC") if sample_2A==1

local tdate_2A = mdy(6, 1, $TY_2A)
gen byte post_2A       = (date >= `tdate_2A')                 if sample_2A==1
gen int  event_week_2A = floor((date - `tdate_2A') / 7)       if sample_2A==1

* Unified event-week FE for baseline DiD time effects (negatives OK since absorbed)
gen int event_week_FE  = event_week_2A if sample_2A==1

di as text "[2A] Sample counts:"
tab province_s    if sample_2A==1
tab treated_BC_2A if sample_2A==1
tab post_2A       if sample_2A==1


* --- Auto event-window selection (TAU rule) ---
* Algorithm: collapse to weekly counts by province; find contiguous sequence
* from week -1 outward where both provinces have >= TAU × core density.
cap macro drop WPRE_2A WPOST_2A WINMETHOD_2A
global WPRE_2A      .
global WPOST_2A     .
global WINMETHOD_2A "MANUAL"

if $AUTO_WIN_2A==1 {
    di as text "[2A] Running auto window selection (TAU=$TAU)..."

    local oksetup = 1
    preserve
        keep if sample_2A==1
        keep if ($SAMPLE_OK)
        keep if !missing($Y)
        keep if inrange(event_week_2A, -$MAXW, $MAXW)
        keep if inlist(province_s, "BC", "$CONTROL_2A")

        collapse (count) N=date, by(event_week_2A province_s)
        reshape wide N, i(event_week_2A) j(province_s) string

        capture confirm variable NBC
        if _rc local oksetup = 0
        local vctrl = "N$CONTROL_2A"
        capture confirm variable `vctrl'
        if _rc local oksetup = 0

        if `oksetup'==1 {
            * Core density reference
            quietly summ NBC         if inrange(event_week_2A, $CORE0, $CORE1)
            scalar coreBC = r(mean)
            quietly summ `vctrl'     if inrange(event_week_2A, $CORE0, $CORE1)
            scalar coreCT = r(mean)

            scalar thrBC = $TAU * coreBC
            scalar thrCT = $TAU * coreCT

            * Week passes if both provinces exceed threshold
            gen byte pass = (NBC >= thrBC) & (`vctrl' >= thrCT)
            sort event_week_2A

            * Expand contiguously pre: start at week -1, walk backward
            gen byte cont_pre = .
            replace cont_pre = pass if event_week_2A == -1
            forvalues k = 2/$MAXW {
                local wk = -`k'
                quietly replace cont_pre = pass & cont_pre[_n+1] if event_week_2A == `wk'
            }

            * Expand contiguously post: start at week 0, walk forward
            gen byte cont_post = .
            replace cont_post = pass if event_week_2A == 0
            forvalues k = 1/$MAXW {
                quietly replace cont_post = pass & cont_post[_n-1] if event_week_2A == `k'
            }

            quietly summ event_week_2A if cont_pre==1,  meanonly
            scalar PRE_sug  = r(min)
            quietly summ event_week_2A if cont_post==1, meanonly
            scalar POST_sug = r(max)

            quietly count if cont_pre==1
            scalar npre = r(N)
            quietly count if cont_post==1
            scalar npost = r(N)
        }
    restore

    * Apply auto result or fall back to manual
    if `oksetup'==0 {
        global WPRE_2A  $WPRE_MAN_2A
        global WPOST_2A $WPOST_MAN_2A
        global WINMETHOD_2A "MANUAL(setup_fail)"
    }
    else if missing(PRE_sug) | npre==0 {
        global WPRE_2A  $WPRE_MAN_2A
        global WPOST_2A $WPOST_MAN_2A
        global WINMETHOD_2A "MANUAL(no_pre)"
    }
    else if missing(POST_sug) | npost==0 {
        global WPRE_2A  $WPRE_MAN_2A
        global WPOST_2A $WPOST_MAN_2A
        global WINMETHOD_2A "MANUAL(no_post)"
    }
    else {
        global WPRE_2A  = -PRE_sug
        global WPOST_2A =  POST_sug
        global WINMETHOD_2A "AUTO(TAU)"
    }
}
else {
    global WPRE_2A  $WPRE_MAN_2A
    global WPOST_2A $WPOST_MAN_2A
    global WINMETHOD_2A "MANUAL"
}

di as text "[2A] Window method: $WINMETHOD_2A"
di as text "[2A] FINAL window:  [-$WPRE_2A, $WPOST_2A]"


* --- 2A Baseline DiD regression ---
di as text "[2A] Running baseline DiD..."
reghdfe $Y ///
    i.treated_BC_2A##i.post_2A ///
    $X ///
    $LIFECYCLE ///
    if sample_2A==1 ///
       & inrange(event_week_2A, -$WPRE_2A, $WPOST_2A) ///
       & ($SAMPLE_OK) ///
       & !missing($Y), ///
    absorb($ABS $TIMEABS) ///
    vce(cluster $CLUSTVAR)
estimates store did_2A
di as text "[2A] DiD: DONE"


* --- 2A Event study ---
* Restrict support to weeks with obs in BOTH BC and control province.
* event_week_s_2A = event_week_2A + 100 so negatives are valid factor values.
cap drop in_esamp_2A bc_obs_2A ctrl_obs_2A bcN_2A ctrlN_2A okweek_2A event_week_s_2A

gen byte in_esamp_2A = (sample_2A==1)
replace  in_esamp_2A = 0 if !inrange(event_week_2A, -$WPRE_2A, $WPOST_2A)
replace  in_esamp_2A = 0 if !($SAMPLE_OK)
replace  in_esamp_2A = 0 if missing($Y)
replace  in_esamp_2A = 0 if !inlist(province_s, "BC", "$CONTROL_2A")

gen byte bc_obs_2A   = (in_esamp_2A==1) & (province_s=="BC")
gen byte ctrl_obs_2A = (in_esamp_2A==1) & (province_s=="$CONTROL_2A")

bysort event_week_2A: egen int bcN_2A   = total(bc_obs_2A)   if !missing(event_week_2A)
bysort event_week_2A: egen int ctrlN_2A = total(ctrl_obs_2A) if !missing(event_week_2A)

* okweek: both sides have support
gen byte okweek_2A        = (in_esamp_2A==1) & (bcN_2A > 0) & (ctrlN_2A > 0)
gen int  event_week_s_2A  = event_week_2A + 100 if okweek_2A==1

* Require week -1 (shifted = 99) to be present
quietly count if in_esamp_2A==1 & event_week_s_2A==99
if r(N)==0 {
    di as error "[2A] FAIL: week -1 not available after support restriction. Adjust windows."
    exit 198
}
global EW_BASE_2A 99

di as text "[2A] Running event study..."
reghdfe $Y ///
    ib($EW_BASE_2A).event_week_s_2A##i.treated_BC_2A ///
    $X ///
    $LIFECYCLE ///
    if in_esamp_2A==1 & !missing(event_week_s_2A), ///
    absorb($ABS) ///
    vce(cluster $CLUSTVAR)
estimates store es_2A


* --- Pre-trend joint test ---
estimates restore es_2A

local pre_levels_2A
quietly levelsof event_week_s_2A ///
    if e(sample) & (event_week_s_2A < $EW_BASE_2A), local(pre_levels_2A)

local testlist_2A ""
foreach k of local pre_levels_2A {
    local testlist_2A `testlist_2A' `k'.event_week_s_2A#1.treated_BC_2A
}

if "`testlist_2A'"=="" {
    local pretxt_2A "Pre(joint) p=NA"
}
else {
    quietly test `testlist_2A'
    local p_pre_joint_2A = r(p)
    local pstr_2A : display %6.3f `p_pre_joint_2A'
    local pretxt_2A "Pre(joint) p=`pstr_2A'"
}
di as text "[2A] `pretxt_2A'"


* --- coefplot ---
levelsof event_week_s_2A if e(sample), local(actual_levs_2A)

local keepcoef_2A ""
local clabel_2A   ""
foreach k of local actual_levs_2A {
    local keepcoef_2A `keepcoef_2A' `k'.event_week_s_2A#1.treated_BC_2A
    local val = `k' - 100
    local clabel_2A `clabel_2A' `k'.event_week_s_2A#1.treated_BC_2A = "`val'"
}

* Position xline between week -1 and week 0 in plotted sequence
preserve
    keep if e(sample)
    keep event_week_s_2A
    duplicates drop
    sort event_week_s_2A
    gen idx = _n
    quietly summ idx if event_week_s_2A==100, meanonly
    local line_pos_2A = r(mean) - 0.5
restore

local sub2A "BC vs $CONTROL_2A | $WINMETHOD_2A | Window [-$WPRE_2A,$WPOST_2A] | TAU=$TAU | `pretxt_2A' | Cluster=$CLUSTVAR"

coefplot es_2A, ///
    keep(`keepcoef_2A') ///
    vertical omitted baselevels ///
    coeflabel(`clabel_2A') ///
    yline(0, lpattern(dash)) ///
    xline(`line_pos_2A', lwidth(vthin)) ///
    title("Event study: $Y ($TY_2A)", size(medium)) ///
    subtitle("`sub2A'", size(small)) ///
    ytitle("Effect relative to week -1") ///
    xtitle("Weeks relative to June 1 cutoff") ///
    graphregion(color(white))

graph export "output/event_2A_${Y}_${TY_2A}_${CONTROL_2A}.png", replace
di as text "[2A] Figure saved: output/event_2A_${Y}_${TY_2A}_${CONTROL_2A}.png"


/**********************************************************************
  2B) WITHIN-BC PLACEBO YEARS

  Strategy: compare BC planters in the treated year ($TY_2B) to BC
  planters in placebo years ($BCPLACEBO_2B). Post = date >= June 1
  of each year (year-specific event time from 1C).

  Two regressions (same as 2A):
    (i)  Baseline DiD with event-week FE for common time shocks.
    (ii) Event study with pre-trend test and coefplot.
**********************************************************************/

di as text "============================================================"
di as text "[2B] Within-BC: treated $TY_2B vs placebo ($BCPLACEBO_2B)"
di as text "============================================================"

* --- Sample + groups ---
cap drop treated_year_2B placebo_year_2B sample_2B post_y_2B event_week_y_2B event_weeky_FE

gen byte treated_year_2B = (year == $TY_2B)

gen byte placebo_year_2B = 0
foreach yy of numlist $BCPLACEBO_2B {
    replace placebo_year_2B = 1 if year == `yy'
}

gen byte sample_2B = (province_s == "BC") & (treated_year_2B==1 | placebo_year_2B==1)

* Use year-specific post and event-week from Section 1C
gen byte post_y_2B       = post_y       if sample_2B==1
gen int  event_week_y_2B = event_week_y if sample_2B==1

* Unified event-week FE for baseline DiD (absorbs common week-level seasonality)
gen int event_weeky_FE = event_week_y_2B if sample_2B==1

di as text "[2B] Sample counts:"
tab year      if sample_2B==1
tab post_y_2B if sample_2B==1


* --- Auto event-window for 2B ---
cap macro drop PRE_MIN_2B POST_MAX_2B WINMETHOD_2B
global PRE_MIN_2B   .
global POST_MAX_2B  .
global WINMETHOD_2B "MANUAL"

if $AUTO_WIN_2B==1 {
    di as text "[2B] Running auto window selection (TAU=$TAU)..."

    local oksetup = 1
    preserve
        keep if sample_2B==1
        keep if ($SAMPLE_OK)
        keep if !missing($Y)
        keep if inrange(event_week_y_2B, -$MAXW, $MAXW)

        collapse (count) N=date, by(event_week_y_2B treated_year_2B)
        reshape wide N, i(event_week_y_2B) j(treated_year_2B)

        capture confirm variable N0
        if _rc local oksetup = 0
        capture confirm variable N1
        if _rc local oksetup = 0

        if `oksetup'==1 {
            quietly summ N1 if inrange(event_week_y_2B, $CORE0, $CORE1)
            scalar coreT = r(mean)
            quietly summ N0 if inrange(event_week_y_2B, $CORE0, $CORE1)
            scalar coreP = r(mean)

            scalar thrT = $TAU * coreT
            scalar thrP = $TAU * coreP

            gen byte pass = (N1 >= thrT) & (N0 >= thrP)
            sort event_week_y_2B

            gen byte cont_pre = .
            replace cont_pre = pass if event_week_y_2B == -1
            forvalues k = 2/$MAXW {
                local wk = -`k'
                quietly replace cont_pre = pass & cont_pre[_n+1] if event_week_y_2B == `wk'
            }

            gen byte cont_post = .
            replace cont_post = pass if event_week_y_2B == 0
            forvalues k = 1/$MAXW {
                quietly replace cont_post = pass & cont_post[_n-1] if event_week_y_2B == `k'
            }

            quietly summ event_week_y_2B if cont_pre==1,  meanonly
            scalar PRE_sug  = r(min)
            quietly summ event_week_y_2B if cont_post==1, meanonly
            scalar POST_sug = r(max)

            quietly count if cont_pre==1
            scalar npre = r(N)
            quietly count if cont_post==1
            scalar npost = r(N)
        }
    restore

    if `oksetup'==0 {
        global PRE_MIN_2B  $PRE_MIN_MAN_2B
        global POST_MAX_2B $POST_MAX_MAN_2B
        global WINMETHOD_2B "MANUAL(setup_fail)"
    }
    else if missing(PRE_sug) | npre==0 {
        global PRE_MIN_2B  $PRE_MIN_MAN_2B
        global POST_MAX_2B $POST_MAX_MAN_2B
        global WINMETHOD_2B "MANUAL(no_pre)"
    }
    else if missing(POST_sug) | npost==0 {
        global PRE_MIN_2B  $PRE_MIN_MAN_2B
        global POST_MAX_2B $POST_MAX_MAN_2B
        global WINMETHOD_2B "MANUAL(no_post)"
    }
    else {
        global PRE_MIN_2B  = PRE_sug
        global POST_MAX_2B = POST_sug
        global WINMETHOD_2B "AUTO(TAU)"
    }
}
else {
    global PRE_MIN_2B  $PRE_MIN_MAN_2B
    global POST_MAX_2B $POST_MAX_MAN_2B
    global WINMETHOD_2B "MANUAL"
}

di as text "[2B] Window method: $WINMETHOD_2B"
di as text "[2B] FINAL window:  [$PRE_MIN_2B, $POST_MAX_2B]"


* --- 2B Baseline DiD regression ---
* Absorbs event_weeky_FE to control for common week effects across years.
di as text "[2B] Running baseline DiD..."
reghdfe $Y ///
    i.treated_year_2B##i.post_y_2B ///
    $X ///
    $LIFECYCLE ///
    if sample_2B==1 ///
       & inrange(event_week_y_2B, $PRE_MIN_2B, $POST_MAX_2B) ///
       & ($SAMPLE_OK) ///
       & !missing($Y), ///
    absorb($ABS event_weeky_FE) ///
    vce(cluster $CLUSTVAR)
estimates store did_2B
di as text "[2B] DiD: DONE"


* --- 2B Event study ---
cap drop in_esamp_2B t_obs_2B p_obs_2B tN_2B pN_2B okweek_2B event_week_sy_2B

gen byte in_esamp_2B = (sample_2B==1)
replace  in_esamp_2B = 0 if !inrange(event_week_y_2B, $PRE_MIN_2B, $POST_MAX_2B)
replace  in_esamp_2B = 0 if !($SAMPLE_OK)
replace  in_esamp_2B = 0 if missing($Y)

gen byte t_obs_2B = (in_esamp_2B==1) & (treated_year_2B==1)
gen byte p_obs_2B = (in_esamp_2B==1) & (treated_year_2B==0)

bysort event_week_y_2B: egen int tN_2B = total(t_obs_2B) if !missing(event_week_y_2B)
bysort event_week_y_2B: egen int pN_2B = total(p_obs_2B) if !missing(event_week_y_2B)

gen byte okweek_2B       = (in_esamp_2B==1) & (tN_2B > 0) & (pN_2B > 0)
gen int  event_week_sy_2B = event_week_y_2B + 100 if okweek_2B==1

quietly count if in_esamp_2B==1 & event_week_sy_2B==99
if r(N)==0 {
    di as error "[2B] FAIL: week -1 not available after support restriction. Adjust windows."
    exit 198
}
global EW_BASE_2B 99

di as text "[2B] Running event study..."
reghdfe $Y ///
    ib($EW_BASE_2B).event_week_sy_2B##i.treated_year_2B ///
    $X ///
    $LIFECYCLE ///
    if in_esamp_2B==1 & !missing(event_week_sy_2B), ///
    absorb($ABS) ///
    vce(cluster $CLUSTVAR)
estimates store es_2B


* --- Pre-trend joint test ---
estimates restore es_2B

local pre_levels_2B
quietly levelsof event_week_sy_2B ///
    if e(sample) & (event_week_sy_2B < $EW_BASE_2B), local(pre_levels_2B)

local testlist_2B ""
foreach k of local pre_levels_2B {
    local testlist_2B `testlist_2B' `k'.event_week_sy_2B#1.treated_year_2B
}

if "`testlist_2B'"=="" {
    local pretxt_2B "Pre(joint) p=NA"
}
else {
    quietly test `testlist_2B'
    local p_pre_joint_2B = r(p)
    local pstr_2B : display %6.3f `p_pre_joint_2B'
    local pretxt_2B "Pre(joint) p=`pstr_2B'"
}
di as text "[2B] `pretxt_2B'"


* --- coefplot ---
levelsof event_week_sy_2B if e(sample), local(actual_levs_2B)

local keepcoef_2B ""
local clabel_2B   ""
foreach k of local actual_levs_2B {
    local keepcoef_2B `keepcoef_2B' `k'.event_week_sy_2B#1.treated_year_2B
    local val = `k' - 100
    local clabel_2B `clabel_2B' `k'.event_week_sy_2B#1.treated_year_2B = "`val'"
}

preserve
    keep if e(sample)
    keep event_week_sy_2B
    duplicates drop
    sort event_week_sy_2B
    gen idx = _n
    quietly summ idx if event_week_sy_2B==100, meanonly
    local line_pos_2B = r(mean) - 0.5
restore

local sub2B "Treated $TY_2B vs placebo ($BCPLACEBO_2B) | $WINMETHOD_2B | Window [$PRE_MIN_2B,$POST_MAX_2B] | TAU=$TAU | `pretxt_2B' | Cluster=$CLUSTVAR"

coefplot es_2B, ///
    keep(`keepcoef_2B') ///
    vertical omitted baselevels ///
    coeflabel(`clabel_2B') ///
    yline(0, lpattern(dash)) ///
    xline(`line_pos_2B', lwidth(vthin)) ///
    title("Event study: $Y (BC only)", size(medium)) ///
    subtitle("`sub2B'", size(small)) ///
    ytitle("Effect relative to week -1") ///
    xtitle("Weeks relative to June 1 cutoff") ///
    graphregion(color(white))

graph export "output/event_2B_${Y}_${TY_2B}.png", replace
di as text "[2B] Figure saved: output/event_2B_${Y}_${TY_2B}.png"


/**********************************************************************
  SUMMARY PRINT
**********************************************************************/

di as text "============================================================"
di as text "SECTION 2: SUMMARY"
di as text "============================================================"

di as text "2A DiD result:"
estimates restore did_2A
lincom 1.treated_BC_2A#1.post_2A

di as text ""
di as text "2B DiD result:"
estimates restore did_2B
lincom 1.treated_year_2B#1.post_y_2B

di as text "============================================================"


/**********************************************************************
  LOG: close
**********************************************************************/

di as text ""
di as text "[02] Finished: $S_DATE $S_TIME"
log close
