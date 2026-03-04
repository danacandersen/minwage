
/**********************************************************************
  MASTER DO-FILE: Planter MW Analysis

  STRUCTURE:
    (0) PACKAGES
    (1) DATA PREP (run once)
        1A) Import, types, IDs, outcomes
        1B) Safe string copies (avoid strL collapse/reshape failures)
        1C) Universal June 1 event-time (year-specific)
        1D) Contract lifecycle + within-season ramp + cross-season experience
    (2) BASELINE ANALYSIS (rerunnable control panel)
        2A) BC vs AB/ON (treated-year cross-province DiD + event study)
        2B) BC placebo-years (within-BC DiD + event study)

  FEATURES:
    - Optional floors (toggle)
    - Auto event-window selection by TAU rule (2A and 2B separately)
      using support in the *actual regression sample*
    - Support-restricted event-study plots (coefplot uses only e(sample))
    - Event-study subtitle reports: AUTO vs MANUAL, window, TAU, pretrend p, cluster var

  REQUIREMENTS:
    - reghdfe (and ftools)
    - coefplot
**********************************************************************/

clear all
set more off


/**********************************************************************
  (0) PACKAGES
**********************************************************************/
capture which reghdfe
if _rc {
    ssc install ftools, replace
    ssc install reghdfe, replace
}
capture which coefplot
if _rc {
    ssc install coefplot, replace
}


/**********************************************************************
  (1) DATA PREP (RUN ONCE)
**********************************************************************/

* --- Load data ---
import delimited using "/Users/danaandersen/Downloads/planter_mw.csv", clear varnames(1)

* --- Required vars check ---
confirm variable date year province contract name piece prod top_up tot_hrs
confirm variable temp total_precip wind_spd experience

* --- Date conversion to Stata daily ---
capture confirm numeric variable date
if _rc {
    gen double date_num = daily(date, "YMD")
    format date_num %td
    drop date
    rename date_num date
}
format date %td

* --- Destring core numeric vars if needed ---
foreach var in temp total_precip wind_spd experience tot_hrs piece prod top_up {
    capture confirm numeric variable `var'
    if _rc destring `var', replace ignore("NA" "." "" "null" "NULL")
}

* --- IDs (create once) ---
cap drop worker_id contract_id wc_id
encode name,     gen(worker_id)
encode contract, gen(contract_id)
egen wc_id = group(worker_id contract_id)

* --- Day-level aggregates to handle multi-contract days ---
cap drop tot_hrs_day n_contracts_day multi_contract_day long_day extreme_day
bysort worker_id date: egen double tot_hrs_day     = total(tot_hrs)
bysort worker_id date: egen int    n_contracts_day = count(contract_id)
gen byte multi_contract_day = (n_contracts_day > 1)
gen byte long_day           = (tot_hrs_day > 14)
gen byte extreme_day        = (tot_hrs_day > 18)

* --- Outcomes (create once) ---
cap drop ln_incentive lprod ln_piece
gen double ln_incentive = ln(piece*prod) if piece>0 & prod>0
gen double lprod        = ln(prod)       if prod>0            // prod = trees/hour
gen double ln_piece     = ln(piece)      if piece>0

* --- Weather transforms (create once) ---
cap drop precip0 precip_pos temp2 wind2 precip2
gen byte   precip0    = (total_precip == 0) if !missing(total_precip)
gen double precip_pos = cond(total_precip > 0 & !missing(total_precip), total_precip, 0)
gen double temp2      = temp^2       if !missing(temp)
gen double wind2      = wind_spd^2   if !missing(wind_spd)
gen double precip2    = precip_pos^2 if !missing(precip_pos)

* --- Experience squared if not already in data ---
capture confirm variable experience2
if _rc gen double experience2 = experience^2 if !missing(experience)


/**********************************************************************
  (1B) SAFE STRING COPIES (AVOID strL ISSUES IN collapse/reshape)
**********************************************************************/
cap drop province_s contract_s name_s
gen str2   province_s = province
gen str40  contract_s = contract
gen str40  name_s     = name


/**********************************************************************
  (1C) UNIVERSAL YEAR-SPECIFIC JUNE 1 EVENT TIME (RUN ONCE)
**********************************************************************/
cap drop tdate_y post_y event_day_y event_week_y
gen double tdate_y      = mdy(6,1,year)
format tdate_y %td
gen byte   post_y       = (date >= tdate_y) if !missing(tdate_y)
gen int    event_day_y  = date - tdate_y if !missing(tdate_y)
gen int    event_week_y = floor(event_day_y/7) if !missing(event_day_y)


/**********************************************************************
  (1D) CONTRACT LIFECYCLE + SEASON RAMP + CROSS-SEASON EXPERIENCE
**********************************************************************/

* Contract lifecycle (within-year)
cap drop contract_start contract_day contract_day2
bysort contract_id year: egen double contract_start = min(date)
gen int contract_day  = date - contract_start
gen int contract_day2 = contract_day^2

* Within-season ramp: count distinct work days within worker×year
cap drop first_day season_days_worked season_days_worked2
bysort worker_id year (date): gen byte first_day = (date != date[_n-1])
bysort worker_id year (date): gen int  season_days_worked = sum(first_day)
gen int season_days_worked2 = season_days_worked^2

* Cross-season experience: number of seasons observed so far (0 in first observed season)
cap drop yrflag cum_seasons
bysort worker_id year: gen byte yrflag = (_n==1)
bysort worker_id (year): gen int cum_seasons = sum(yrflag) - 1
drop yrflag


/**********************************************************************
  (2) BASELINE ANALYSIS (RERUNNABLE CONTROL PANEL)
**********************************************************************/

* Clean rerun state (safe)
cap estimates clear
cap label drop ew_lbl ew_lbl2

foreach v in ///
    sample_2A treated_BC_2A post_2A event_week_2A ///
    in_esamp_2A bc_obs_2A ctrl_obs_2A bcN_2A ctrlN_2A okweek_2A event_week_s_2A ///
    sample_2B treated_year_2B placebo_year_2B post_y_2B event_week_y_2B ///
    in_esamp_2B t_obs_2B p_obs_2B tN_2B pN_2B okweek_2B event_week_sy_2B {
    cap drop `v'
}

foreach m in ///
    X ABS SAMPLE_OK LIFECYCLE ///
    WPRE_2A WPOST_2A PRE_MIN_2B POST_MAX_2B ///
    EW_BASE_2A EW_BASE_2B ///
    WINMETHOD_2A WINMETHOD_2B {
    cap macro drop `m'
}


/**********************************************************************
  2.1) CONTROL PANEL (EDIT THESE)
**********************************************************************/

* Outcome variable (created in Section 1)
global Y ln_incentive    // ln_incentive / lprod / ln_piece / tot_hrs_day

* Controls toggles
global USE_WEATHER     1
global USE_EXPER       1
global USE_HOURS       1
global USE_DAYFLAGS    1
global USE_LIFECYCLE   1     // contract ramp + season ramp + ramp×cum_seasons
global USE_TIMEFE      1    // 1 = include event-week FE, 0 = no time FE

* Fixed effects choice:
*   worker_contract -> absorb(worker_id contract_id)
*   wc_only         -> absorb(wc_id)
*   worker_only     -> absorb(worker_id)
*   contract_only   -> absorb(contract_id)
global FE worker_contract

* Floors (toggle)
global USE_FLOORS   1
global MIN_HRS      2
global MIN_PROD     50
global MIN_EARN     30
global USE_MIN_PROD 1

* AUTO window selection (TAU rule)
global TAU   0.25
global CORE0 0
global CORE1 1
global MAXW  30

global AUTO_WIN_2A 1
global AUTO_WIN_2B 1

* Manual windows used if AUTO is off or fails
global WPRE_MAN_2A    3
global WPOST_MAN_2A   8
global PRE_MIN_MAN_2B -6
global POST_MAX_MAN_2B 10

* --------------------------
* CLUSTERING TOGGLE (USED EVERYWHERE)
* --------------------------
* Choose ONE: contract_id, worker_id, wc_id
global CLUSTVAR worker_id

capture confirm variable $CLUSTVAR
if _rc {
    di as error "[FATAL] CLUSTVAR=$CLUSTVAR not found. Set global CLUSTVAR to an existing variable."
    exit 198
}


/**********************************************************************
  2.2) BUILD CONTROLS ($X), LIFECYCLE ($LIFECYCLE), FE ($ABS), FLOORS ($SAMPLE_OK)
**********************************************************************/

* Main controls
global X
if $USE_HOURS==1    global X "$X tot_hrs_day"
if $USE_DAYFLAGS==1 global X "$X multi_contract_day"
if $USE_EXPER==1    global X "$X c.experience##c.experience"
if $USE_WEATHER==1  global X "$X c.temp##c.temp c.wind_spd##c.wind_spd precip0 c.precip_pos##c.precip_pos"

* Lifecycle/ramp controls (built safely)
global LIFECYCLE
if $USE_LIFECYCLE==1 {
    global LIFECYCLE "$LIFECYCLE c.contract_day##c.contract_day"
    global LIFECYCLE "$LIFECYCLE c.season_days_worked##c.season_days_worked"
    global LIFECYCLE "$LIFECYCLE c.season_days_worked##c.cum_seasons"
}

* Fixed effects
global ABS
if "$FE"=="worker_contract" global ABS "worker_id contract_id"
if "$FE"=="wc_only"         global ABS "wc_id"
if "$FE"=="worker_only"     global ABS "worker_id"
if "$FE"=="contract_only"   global ABS "contract_id"

* ---- Time FE absorb list (empty unless toggled on) ----
global TIMEABS ""
if $USE_TIMEFE==1 global TIMEABS "event_week_FE"
di as text "[DEBUG] TIMEABS:     $TIMEABS"


* Floors / inclusion
global SAMPLE_OK "1==1"
if $USE_FLOORS==1 {
    global SAMPLE_OK "$SAMPLE_OK & tot_hrs_day >= $MIN_HRS"
    global SAMPLE_OK "$SAMPLE_OK & piece > 0 & prod > 0"
    global SAMPLE_OK "$SAMPLE_OK & (piece*prod) >= $MIN_EARN"
    if $USE_MIN_PROD==1 global SAMPLE_OK "$SAMPLE_OK & prod >= $MIN_PROD"
}

di as text "[DEBUG] Y:         $Y"
di as text "[DEBUG] X:         $X"
di as text "[DEBUG] LIFECYCLE: $LIFECYCLE"
di as text "[DEBUG] ABS:       $ABS"
di as text "[DEBUG] FLOORS:    $SAMPLE_OK"
di as text "[DEBUG] CLUSTER:   $CLUSTVAR"



/**********************************************************************
  2A) BC vs AB/ON — treated year only
  UPDATED:
    - Defines a unified event-week FE variable (event_week_FE) for common time effects
      so baseline DiD can use absorb(... event_week_FE) safely (negatives allowed).
    - Initializes window globals + a window-method flag (WINMETHOD_2A) cleanly.
**********************************************************************/

* --------------------------
* 2A settings
* --------------------------
global TY_2A 2018
global CONTROL_2A "AB"      // "AB" or "ON" (ABON allowed but auto-window off)

if "$CONTROL_2A"=="ABON" global AUTO_WIN_2A 0

* --------------------------
* 2A sample + event time
* --------------------------
cap drop sample_2A treated_BC_2A post_2A event_week_2A event_week_FE
gen byte sample_2A = (year==$TY_2A) & inlist(province_s,"BC","$CONTROL_2A")
if "$CONTROL_2A"=="ABON" replace sample_2A = (year==$TY_2A) & inlist(province_s,"BC","AB","ON")

gen byte treated_BC_2A = .
replace treated_BC_2A = (province_s=="BC") if sample_2A==1

local tdate_2A = mdy(6,1,$TY_2A)
gen byte post_2A       = (date >= `tdate_2A') if sample_2A==1
gen int  event_week_2A = floor((date - `tdate_2A')/7) if sample_2A==1

* Unified event-week FE variable for baseline DiD time effects (absorbed; safe with negatives)
gen int event_week_FE = event_week_2A if sample_2A==1

* sanity checks
tab province_s    if sample_2A==1
tab treated_BC_2A if sample_2A==1
tab post_2A       if sample_2A==1

* --------------------------
* 2A window globals (auto-window block will overwrite these)
* --------------------------
cap macro drop WPRE_2A WPOST_2A WINMETHOD_2A
global WPRE_2A      .
global WPOST_2A     .
global WINMETHOD_2A "MANUAL"   // set to "TAU" inside auto-window block if it succeeds

if $AUTO_WIN_2A==1 {

    local oksetup = 1
    preserve
        keep if sample_2A==1
        keep if ($SAMPLE_OK)
        keep if !missing($Y)
        keep if inrange(event_week_2A, -$MAXW, $MAXW)
        keep if inlist(province_s,"BC","$CONTROL_2A")

        collapse (count) N=date, by(event_week_2A province_s)
        reshape wide N, i(event_week_2A) j(province_s) string

        capture confirm variable NBC
        if _rc local oksetup = 0
        local vctrl = "N$CONTROL_2A"
        capture confirm variable `vctrl'
        if _rc local oksetup = 0

        if `oksetup'==1 {
            quietly summarize NBC if inrange(event_week_2A,$CORE0,$CORE1)
            scalar coreBC = r(mean)
            quietly summarize `vctrl' if inrange(event_week_2A,$CORE0,$CORE1)
            scalar coreCT = r(mean)

            scalar thrBC = $TAU*coreBC
            scalar thrCT = $TAU*coreCT

            gen byte pass = (NBC>=thrBC) & (`vctrl'>=thrCT)
            sort event_week_2A

            gen byte cont_pre = .
            replace cont_pre = pass if event_week_2A==-1
            forvalues k=2/$MAXW {
                local wk = -`k'
                quietly replace cont_pre = pass & cont_pre[_n+1] if event_week_2A==`wk'
            }

            gen byte cont_post = .
            replace cont_post = pass if event_week_2A==0
            forvalues k=1/$MAXW {
                quietly replace cont_post = pass & cont_post[_n-1] if event_week_2A==`k'
            }

            quietly summarize event_week_2A if cont_pre==1, meanonly
            scalar PRE_sug = r(min)
            quietly summarize event_week_2A if cont_post==1, meanonly
            scalar POST_sug = r(max)

            quietly count if cont_pre==1
            scalar npre = r(N)
            quietly count if cont_post==1
            scalar npost = r(N)
        }
    restore

    if (`oksetup'==0) {
        global WPRE_2A  $WPRE_MAN_2A
        global WPOST_2A $WPOST_MAN_2A
        global WINMETHOD_2A MANUAL
    }
    else if missing(PRE_sug) | npre==0 {
        global WPRE_2A  $WPRE_MAN_2A
        global WPOST_2A $WPOST_MAN_2A
        global WINMETHOD_2A MANUAL
    }
    else if missing(POST_sug) | npost==0 {
        global WPRE_2A  $WPRE_MAN_2A
        global WPOST_2A $WPOST_MAN_2A
        global WINMETHOD_2A MANUAL
    }
    else {
        global WPRE_2A  = -PRE_sug
        global WPOST_2A =  POST_sug
        global WINMETHOD_2A AUTO(TAU)
    }
}
else {
    global WPRE_2A  $WPRE_MAN_2A
    global WPOST_2A $WPOST_MAN_2A
    global WINMETHOD_2A MANUAL
}

di as text "[2A] Window method: $WINMETHOD_2A"
di as text "[2A] FINAL window: [-$WPRE_2A,$WPOST_2A]"

* 2A DiD regression
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

* 2A Event study (support-restricted) + pretrend test + plot
cap drop in_esamp_2A bc_obs_2A ctrl_obs_2A bcN_2A ctrlN_2A okweek_2A event_week_s_2A
gen byte in_esamp_2A = (sample_2A==1)
replace in_esamp_2A = 0 if !inrange(event_week_2A, -$WPRE_2A, $WPOST_2A)
replace in_esamp_2A = 0 if !($SAMPLE_OK)
replace in_esamp_2A = 0 if missing($Y)
replace in_esamp_2A = 0 if !inlist(province_s,"BC","$CONTROL_2A")

gen byte bc_obs_2A   = (in_esamp_2A==1) & (province_s=="BC")
gen byte ctrl_obs_2A = (in_esamp_2A==1) & (province_s=="$CONTROL_2A")

bysort event_week_2A: egen int bcN_2A   = total(bc_obs_2A)   if !missing(event_week_2A)
bysort event_week_2A: egen int ctrlN_2A = total(ctrl_obs_2A) if !missing(event_week_2A)

gen byte okweek_2A = (in_esamp_2A==1) & (bcN_2A>0) & (ctrlN_2A>0)
gen int  event_week_s_2A = event_week_2A + 100 if okweek_2A==1

quietly count if in_esamp_2A==1 & event_week_s_2A==99
if r(N)==0 {
    di as error "[2A] FAIL: week -1 not available after support restriction. Adjust windows."
    exit 198
}
global EW_BASE_2A 99

reghdfe $Y ///
    ib($EW_BASE_2A).event_week_s_2A##i.treated_BC_2A ///
    $X ///
    $LIFECYCLE ///
    if in_esamp_2A==1 ///
       & !missing(event_week_s_2A), ///
    absorb($ABS) ///
    vce(cluster $CLUSTVAR)
estimates store es_2A

* Pre-trend joint test (all pre-period treatment coefficients)
estimates restore es_2A

local pre_levels_2A
quietly levelsof event_week_s_2A if e(sample) & (event_week_s_2A < $EW_BASE_2A), local(pre_levels_2A)

local testlist_2A
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

* coefplot: only coefficients present in e(sample)
levelsof event_week_s_2A if e(sample), local(actual_levs_2A)

local keepcoef_2A
local clabel_2A
foreach k of local actual_levs_2A {
    local keepcoef_2A `keepcoef_2A' `k'.event_week_s_2A#1.treated_BC_2A
    local val = `k' - 100
    local clabel_2A `clabel_2A' `k'.event_week_s_2A#1.treated_BC_2A = "`val'"
}

* xline position: between -1 and 0 in the plotted sequence
preserve
    keep if e(sample)
    keep event_week_s_2A
    duplicates drop
    sort event_week_s_2A
    gen idx = _n
    quietly summarize idx if event_week_s_2A==100, meanonly
    local line_pos_2A = r(mean) - 0.5
restore

local sub2A "BC vs $CONTROL_2A | $WINMETHOD_2A | Window [-$WPRE_2A,$WPOST_2A] | TAU=$TAU | `pretxt_2A' | Cluster=$CLUSTVAR"

coefplot es_2A, ///
    keep(`keepcoef_2A') ///
    vertical ///
    omitted baselevels ///
    coeflabel(`clabel_2A') ///
    yline(0, lpattern(dash)) ///
    xline(`line_pos_2A', lwidth(vthin)) ///
    title("Event study: $Y ($TY_2A)", size(medium)) ///
    subtitle("`sub2A'", size(small)) ///
    ytitle("Effect relative to week -1") ///
    xtitle("Weeks relative to June 1 cutoff") ///
    graphregion(color(white))

graph export "event_2A_${Y}_${TY_2A}_${CONTROL_2A}.png", replace



/**********************************************************************
  2B) BC placebo-years — within BC
  UPDATED:
    - Adds a unified event-week FE variable (event_weeky_FE) for common time effects
      so baseline DiD can use absorb(... event_weeky_FE) safely (negatives allowed).
    - Initializes window globals + a window-method flag (WINMETHOD_2B) cleanly.
**********************************************************************/

* --------------------------
* 2B settings
* --------------------------
global TY_2B 2018
global BCPLACEBO_2B "2016 2017"

* --------------------------
* 2B sample + groups
* --------------------------
cap drop treated_year_2B placebo_year_2B sample_2B post_y_2B event_week_y_2B event_weeky_FE

gen byte treated_year_2B = (year==$TY_2B)

gen byte placebo_year_2B = 0
foreach yy of numlist $BCPLACEBO_2B {
    replace placebo_year_2B = 1 if year==`yy'
}

gen byte sample_2B = (province_s=="BC") & (treated_year_2B==1 | placebo_year_2B==1)

* year-specific event time already exists as post_y / event_week_y from Section 1C
gen byte post_y_2B        = post_y       if sample_2B==1
gen int  event_week_y_2B  = event_week_y if sample_2B==1

* Unified event-week FE variable for baseline DiD time effects (absorbed; safe with negatives)
gen int event_weeky_FE = event_week_y_2B if sample_2B==1

* sanity checks
tab year   if sample_2B==1
tab post_y_2B if sample_2B==1

* --------------------------
* 2B window globals (auto-window block will overwrite these)
* --------------------------
cap macro drop PRE_MIN_2B POST_MAX_2B WINMETHOD_2B
global PRE_MIN_2B      .
global POST_MAX_2B     .
global WINMETHOD_2B    "MANUAL"   // set to "TAU" inside auto-window block if it succeeds

if $AUTO_WIN_2B==1 {

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
            quietly summarize N1 if inrange(event_week_y_2B,$CORE0,$CORE1)
            scalar coreT = r(mean)
            quietly summarize N0 if inrange(event_week_y_2B,$CORE0,$CORE1)
            scalar coreP = r(mean)

            scalar thrT = $TAU*coreT
            scalar thrP = $TAU*coreP

            gen byte pass = (N1>=thrT) & (N0>=thrP)
            sort event_week_y_2B

            gen byte cont_pre = .
            replace cont_pre = pass if event_week_y_2B==-1
            forvalues k=2/$MAXW {
                local wk = -`k'
                quietly replace cont_pre = pass & cont_pre[_n+1] if event_week_y_2B==`wk'
            }

            gen byte cont_post = .
            replace cont_post = pass if event_week_y_2B==0
            forvalues k=1/$MAXW {
                quietly replace cont_post = pass & cont_post[_n-1] if event_week_y_2B==`k'
            }

            quietly summarize event_week_y_2B if cont_pre==1, meanonly
            scalar PRE_sug = r(min)
            quietly summarize event_week_y_2B if cont_post==1, meanonly
            scalar POST_sug = r(max)

            quietly count if cont_pre==1
            scalar npre = r(N)
            quietly count if cont_post==1
            scalar npost = r(N)
        }
    restore

    if (`oksetup'==0) {
        global PRE_MIN_2B  $PRE_MIN_MAN_2B
        global POST_MAX_2B $POST_MAX_MAN_2B
        global WINMETHOD_2B MANUAL
    }
    else if missing(PRE_sug) | npre==0 {
        global PRE_MIN_2B  $PRE_MIN_MAN_2B
        global POST_MAX_2B $POST_MAX_MAN_2B
        global WINMETHOD_2B MANUAL
    }
    else if missing(POST_sug) | npost==0 {
        global PRE_MIN_2B  $PRE_MIN_MAN_2B
        global POST_MAX_2B $POST_MAX_MAN_2B
        global WINMETHOD_2B MANUAL
    }
    else {
        global PRE_MIN_2B  = PRE_sug
        global POST_MAX_2B = POST_sug
        global WINMETHOD_2B AUTO(TAU)
    }
}
else {
    global PRE_MIN_2B  $PRE_MIN_MAN_2B
    global POST_MAX_2B $POST_MAX_MAN_2B
    global WINMETHOD_2B MANUAL
}

di as text "[2B] Window method: $WINMETHOD_2B"
di as text "[2B] FINAL window: [$PRE_MIN_2B,$POST_MAX_2B]"

* 2B DiD regression (UPDATED: adds common time effects via event-week FE)
* NOTE: event_weeky_FE created in the updated 2B block above.

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

* 2B Event study (support-restricted) + pretrend test + plot
cap drop in_esamp_2B t_obs_2B p_obs_2B tN_2B pN_2B okweek_2B event_week_sy_2B
gen byte in_esamp_2B = (sample_2B==1)
replace in_esamp_2B = 0 if !inrange(event_week_y_2B, $PRE_MIN_2B, $POST_MAX_2B)
replace in_esamp_2B = 0 if !($SAMPLE_OK)
replace in_esamp_2B = 0 if missing($Y)

gen byte t_obs_2B = (in_esamp_2B==1) & (treated_year_2B==1)
gen byte p_obs_2B = (in_esamp_2B==1) & (treated_year_2B==0)

bysort event_week_y_2B: egen int tN_2B = total(t_obs_2B) if !missing(event_week_y_2B)
bysort event_week_y_2B: egen int pN_2B = total(p_obs_2B) if !missing(event_week_y_2B)

gen byte okweek_2B = (in_esamp_2B==1) & (tN_2B>0) & (pN_2B>0)
gen int  event_week_sy_2B = event_week_y_2B + 100 if okweek_2B==1

quietly count if in_esamp_2B==1 & event_week_sy_2B==99
if r(N)==0 {
    di as error "[2B] FAIL: week -1 not available after support restriction. Adjust windows."
    exit 198
}
global EW_BASE_2B 99

reghdfe $Y ///
    ib($EW_BASE_2B).event_week_sy_2B##i.treated_year_2B ///
    $X ///
    $LIFECYCLE ///
    if in_esamp_2B==1 ///
       & !missing(event_week_sy_2B), ///
    absorb($ABS) ///
    vce(cluster $CLUSTVAR)
estimates store es_2B

* Pre-trend joint test (all pre-period treated-year coefficients)
estimates restore es_2B

local pre_levels_2B
quietly levelsof event_week_sy_2B if e(sample) & ///
    (event_week_sy_2B < $EW_BASE_2B), ///
    local(pre_levels_2B)

local testlist_2B
foreach k of local pre_levels_2B {
    local testlist_2B `testlist_2B' ///
        `k'.event_week_sy_2B#1.treated_year_2B
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

* coefplot: only coefficients present in e(sample)
levelsof event_week_sy_2B if e(sample), local(actual_levs_2B)

local keepcoef_2B
local clabel_2B
foreach k of local actual_levs_2B {
    local keepcoef_2B `keepcoef_2B' `k'.event_week_sy_2B#1.treated_year_2B
    local val = `k' - 100
    local clabel_2B `clabel_2B' `k'.event_week_sy_2B#1.treated_year_2B = "`val'"
}

* xline position: between -1 and 0 in the plotted sequence
preserve
    keep if e(sample)
    keep event_week_sy_2B
    duplicates drop
    sort event_week_sy_2B
    gen idx = _n
    quietly summarize idx if event_week_sy_2B==100, meanonly
    local line_pos_2B = r(mean) - 0.5
restore

local sub2B "Treated $TY_2B vs placebo ($BCPLACEBO_2B) | $WINMETHOD_2B | Window [$PRE_MIN_2B,$POST_MAX_2B] | TAU=$TAU | `pretxt_2B' | Cluster=$CLUSTVAR"

coefplot es_2B, ///
    keep(`keepcoef_2B') ///
    vertical ///
    omitted baselevels ///
    coeflabel(`clabel_2B') ///
    yline(0, lpattern(dash)) ///
    xline(`line_pos_2B', lwidth(vthin)) ///
    title("Event study: $Y (BC only)", size(medium)) ///
    subtitle("`sub2B'", size(small)) ///
    ytitle("Effect relative to week -1") ///
    xtitle("Weeks relative to June 1 cutoff") ///
    graphregion(color(white))

graph export "event_2B_${Y}_${TY_2B}.png", replace

/**********************************************************************
  (3) ADDITIONAL EMPIRICAL ANALYSIS (SELF-CONTAINED)
  - Single, unified toggle system that applies to ALL DiD/heterogeneity
    regressions in Section 3 (3.9, 3.10, 3.11, 3.12, 3.14 printing).
  - Fixes: invalid name errors from $toggles, fragile ## expansion,
    and inconsistent absorb logic across subsections.
**********************************************************************/

***********************************************************************
* 3.0) CLEAN RERUN STATE (SAFE)
***********************************************************************
cap estimates clear

* ---- Drop generated variables (safe even if they don't exist)
foreach v in ///
    sample_3 treated_BC_3 post_3 event_week_3 event_week_3s ///
    inwin_3 ///
    ln_prod abil_sample resid_preprod ability_raw ability_z ///
    wx_shock ///
    cd_sample cd_sample_ab cd_sample_bc ///
    resid_cd_ab resid_cd_bc resid_cd ///
    contract_diff_raw contract_diff_z ///
    topup_ind topup_pre phat_topup ///
    phat_mu_w phat_sd_w phat_mu_wc phat_sd_wc ///
    phat_mu_bin phat_sd_bin ///
    abil_bin_all ///
    est11 minpost11 maxpost11 spans11 any11 provBC tagc tagcs ///
{
    cap drop `v'
}

* ---- Drop globals/macros used in Section 3
foreach m in ///
    TY_3 CONTROL_3 Y_HET ///
    CLUSTER_3 VCE_3 ///
    USE_FLOORS_3 MIN_HRS_3 MIN_EARN_3 MIN_PROD_3 USE_MIN_PROD_3 ///
    USE_WEATHER_3 USE_EXPER_3 USE_HOURS_3 USE_DAYFLAGS_3 USE_LIFECYCLE_3 ///
    WPRE_3 WPOST_3 ///
    CD_POST_MODE CD_TIMEFE ///
    PHAT_VAR_LEVEL PHAT_BINS_N PHAT_SD_BINS ///
    DID_CONTRACTFE_3 DID_WEEKFE_3 DID_POST_MAIN_3 DID_POSTX_3 ///
    X_3 LIFECYCLE_3 SAMPLE_OK_3 ///
{
    cap macro drop `m'
}

***********************************************************************
* 3.1) SETTINGS (EDIT THESE)
***********************************************************************
global TY_3      2018
global CONTROL_3 "AB"
global Y_HET     "ln_incentive"

global CLUSTER_3 "worker_id"
global VCE_3     "vce(cluster $CLUSTER_3)"

* Controls toggles
global USE_WEATHER_3   1
global USE_EXPER_3     1
global USE_HOURS_3     0
global USE_DAYFLAGS_3  0
global USE_LIFECYCLE_3 1

* Floors toggle (0 = none)
global USE_FLOORS_3   0
global MIN_HRS_3      2
global MIN_EARN_3     30
global USE_MIN_PROD_3 1
global MIN_PROD_3     50

* Event window
global WPRE_3  3
global WPOST_3 8

* Contract difficulty timing:
* 0 = PRE only (AB+BC)
* 1 = AB pre+post, BC pre only   (recommended)
* 2 = POST allowed for all (AB+BC)
global CD_POST_MODE 1

* Week FE in difficulty construction
global CD_TIMEFE 0

* phat moments level
global PHAT_VAR_LEVEL "wc"

* phat bins
global PHAT_BINS_N 4
global PHAT_SD_BINS 0

* ============================================================
* GLOBAL SPEC TOGGLES (apply to ALL Section 3 DiD regressions)
* ============================================================
* Include contract FE in DiD regressions?
global DID_CONTRACTFE_3 1

* Include week FE in DiD regressions?
global DID_WEEKFE_3 1

* Include post main effect explicitly?
global DID_POST_MAIN_3 1

* Include common post×X terms for heterogeneity covariates X?
global DID_POSTX_3 1


***********************************************************************
* 3.2) BUILD CONTROLS + FLOORS STRINGS
***********************************************************************
global X_3 ""
if $USE_HOURS_3==1    global X_3 "$X_3 tot_hrs_day"
if $USE_DAYFLAGS_3==1 global X_3 "$X_3 multi_contract_day"
if $USE_EXPER_3==1    global X_3 "$X_3 c.experience##c.experience"
if $USE_WEATHER_3==1  global X_3 "$X_3 c.temp##c.temp c.wind_spd##c.wind_spd precip0 c.precip_pos##c.precip_pos"

global LIFECYCLE_3 ""
if $USE_LIFECYCLE_3==1 {
    global LIFECYCLE_3 "$LIFECYCLE_3 c.contract_day##c.contract_day"
    global LIFECYCLE_3 "$LIFECYCLE_3 c.season_days_worked##c.season_days_worked"
    global LIFECYCLE_3 "$LIFECYCLE_3 c.season_days_worked##c.cum_seasons"
}

global SAMPLE_OK_3 "1==1"
if $USE_FLOORS_3==1 {
    global SAMPLE_OK_3 "$SAMPLE_OK_3 & tot_hrs_day >= $MIN_HRS_3"
    global SAMPLE_OK_3 "$SAMPLE_OK_3 & piece > 0 & prod > 0"
    global SAMPLE_OK_3 "$SAMPLE_OK_3 & (piece*prod) >= $MIN_EARN_3"
    if $USE_MIN_PROD_3==1 global SAMPLE_OK_3 "$SAMPLE_OK_3 & prod >= $MIN_PROD_3"
}

di as text "------------------------------------------------------------"
di as text "[3] TY_3:           $TY_3   CONTROL_3: $CONTROL_3   Y_HET: $Y_HET"
di as text "[3] Window:         [-$WPRE_3,$WPOST_3]"
di as text "[3] CD_POST_MODE:   $CD_POST_MODE   CD_TIMEFE: $CD_TIMEFE"
di as text "[3] PHAT_VAR_LEVEL: $PHAT_VAR_LEVEL   PHAT_BINS_N: $PHAT_BINS_N   PHAT_SD_BINS: $PHAT_SD_BINS"
di as text "[3] DID_CONTRACTFE: $DID_CONTRACTFE_3   DID_WEEKFE: $DID_WEEKFE_3   DID_POST_MAIN: $DID_POST_MAIN_3   DID_POSTX: $DID_POSTX_3"
di as text "------------------------------------------------------------"

* ---- Safe numeric locals from globals (avoids ==1 invalid name)
local did_contractfe = real(trim("$DID_CONTRACTFE_3"))
local did_weekfe     = real(trim("$DID_WEEKFE_3"))
local did_post_main  = real(trim("$DID_POST_MAIN_3"))
local did_postx       = real(trim("$DID_POSTX_3"))
local cdmode          = real(trim("$CD_POST_MODE"))
local cd_timefe       = real(trim("$CD_TIMEFE"))

***********************************************************************
* 3.3) DEFINE SAMPLE + EVENT TIME
***********************************************************************
gen byte sample_3 = (year==$TY_3) & inlist(province_s,"BC","$CONTROL_3")

gen byte treated_BC_3 = .
replace treated_BC_3 = (province_s=="BC") if sample_3==1

local tdate_3 = mdy(6,1,$TY_3)
gen byte post_3       = (date >= `tdate_3') if sample_3==1
gen int  event_week_3 = floor((date - `tdate_3')/7) if sample_3==1
gen int  event_week_3s = event_week_3 + 100 if sample_3==1

tab province_s if sample_3==1
tab post_3     if sample_3==1

* ---- Unified absorb list for ALL DiD regressions in Section 3
local ABS_DID_3 "worker_id"
if `did_contractfe'==1 local ABS_DID_3 "`ABS_DID_3' contract_id"
if `did_weekfe'==1     local ABS_DID_3 "`ABS_DID_3' event_week_3s"

**********************************************************************
* 3.4) PRE-TREATMENT ABILITY (BC+AB, PRE ONLY)
**********************************************************************
cap drop ln_prod
gen double ln_prod = ln(prod) if prod>0

cap drop abil_sample
gen byte abil_sample = (sample_3==1) & (post_3==0) & !missing(ln_prod) & ($SAMPLE_OK_3)

cap drop resid_preprod
reghdfe ln_prod $X_3 $LIFECYCLE_3 if abil_sample==1, absorb(contract_id) resid(resid_preprod)

cap drop ability_raw ability_z
bysort worker_id: egen double ability_raw = mean(resid_preprod)

quietly summarize ability_raw if !missing(ability_raw)
scalar abil_mu = r(mean)
scalar abil_sd = r(sd)

gen double ability_z = (ability_raw - abil_mu)/abil_sd if !missing(ability_raw) & abil_sd>0

**********************************************************************
* 3.5) DAILY WEATHER SHOCK (reduced-form)
**********************************************************************
cap drop wx_shock
quietly reg ln_prod ///
    c.temp##c.temp c.wind_spd##c.wind_spd precip0 c.precip_pos##c.precip_pos ///
    if sample_3==1 & !missing(ln_prod), vce(robust)
predict double wx_shock if sample_3==1, xb

/**********************************************************************
* 3.6) CONTRACT DIFFICULTY (AB + BC, consistent CD_POST_MODE)
*   Goal: contract_diff_z defined for BOTH AB and BC contracts.
**********************************************************************/
cap drop cd_sample_ab cd_sample_bc cd_sample resid_cd_ab resid_cd_bc resid_cd contract_diff_raw contract_diff_z

gen byte cd_sample_ab = (year==$TY_3) & (province_s=="AB") & !missing(ln_prod) & !missing(ability_z) & ($SAMPLE_OK_3)
gen byte cd_sample_bc = (year==$TY_3) & (province_s=="BC") & !missing(ln_prod) & !missing(ability_z) & ($SAMPLE_OK_3)

if `cdmode'==0 {
    replace cd_sample_ab = cd_sample_ab & (post_3==0)
    replace cd_sample_bc = cd_sample_bc & (post_3==0)
}
else if `cdmode'==1 {
    * AB pre+post allowed; BC pre only
    replace cd_sample_bc = cd_sample_bc & (post_3==0)
}
else if `cdmode'==2 {
    * allow post for all: do nothing
}

gen byte cd_sample = cd_sample_ab | cd_sample_bc

* Absorb list for difficulty residualization
local ABS_CD_3 "worker_id"
if `cd_timefe'==1 local ABS_CD_3 "`ABS_CD_3' event_week_3s"

reghdfe ln_prod c.ability_z $X_3 $LIFECYCLE_3 if cd_sample_ab==1, absorb(`ABS_CD_3') resid(resid_cd_ab)
reghdfe ln_prod c.ability_z $X_3 $LIFECYCLE_3 if cd_sample_bc==1, absorb(`ABS_CD_3') resid(resid_cd_bc)

gen double resid_cd = .
replace resid_cd = resid_cd_ab if !missing(resid_cd_ab)
replace resid_cd = resid_cd_bc if !missing(resid_cd_bc)

bysort contract_id: egen double contract_diff_raw = mean(resid_cd)

quietly summarize contract_diff_raw if !missing(contract_diff_raw)
scalar cd_mu = r(mean)
scalar cd_sd = r(sd)

gen double contract_diff_z = (contract_diff_raw - cd_mu)/cd_sd if !missing(contract_diff_raw) & cd_sd>0

di as text "[3.6] Contract difficulty coverage in sample_3:"
count if sample_3==1 & province_s=="AB" & !missing(contract_diff_z)
count if sample_3==1 & province_s=="BC" & !missing(contract_diff_z)

**********************************************************************
* 3.7) PREDICTED TOP-UP RISK phat_topup (PRE only; predict all)
**********************************************************************
cap drop topup_ind topup_pre phat_topup
gen byte topup_ind = (top_up > 0) if !missing(top_up)

gen byte topup_pre = (sample_3==1) & (post_3==0) ///
    & !missing(topup_ind) & !missing(ability_z) & !missing(mw) & !missing(wx_shock) ///
    & ($SAMPLE_OK_3)

reg topup_ind c.mw c.ability_z c.wx_shock $X_3 i.contract_id if topup_pre==1, vce(cluster $CLUSTER_3)
predict double phat_topup if sample_3==1, xb
replace phat_topup = max(0, min(1, phat_topup)) if sample_3==1 & !missing(phat_topup)

**********************************************************************
* 3.8) PHAT MOMENTS + ANALYSIS WINDOW FLAG (inwin_3)
**********************************************************************
cap drop inwin_3
gen byte inwin_3 = (sample_3==1) ///
    & inrange(event_week_3,-$WPRE_3,$WPOST_3) ///
    & ($SAMPLE_OK_3) ///
    & !missing($Y_HET) ///
    & !missing(phat_topup)

cap drop phat_mu_w phat_sd_w phat_mu_wc phat_sd_wc

if "$PHAT_VAR_LEVEL"=="worker" {
    bysort worker_id: egen double phat_mu_w = mean(cond(inwin_3==1, phat_topup, .))
    bysort worker_id: egen double phat_sd_w = sd(  cond(inwin_3==1, phat_topup, .))
}
else {
    capture confirm variable wc_id
    if _rc {
        di as error "[3.8] wc_id not found -> falling back to worker moments"
        bysort worker_id: egen double phat_mu_w = mean(cond(inwin_3==1, phat_topup, .))
        bysort worker_id: egen double phat_sd_w = sd(  cond(inwin_3==1, phat_topup, .))
        global PHAT_VAR_LEVEL "worker"
    }
    else {
        bysort wc_id: egen double phat_mu_wc = mean(cond(inwin_3==1, phat_topup, .))
        bysort wc_id: egen double phat_sd_wc = sd(  cond(inwin_3==1, phat_topup, .))
    }
}

**********************************************************************
* 3.8B) VISUALIZE DISTRIBUTION OF phat (in window)
*   (simple, copy/paste friendly)
**********************************************************************
di as text "------------------------------------------------------------"
di as text "[3.8B] phat distribution (inwin_3==1):"
quietly count if inwin_3==1
di as text "   N(inwin)= " r(N)

if "$PHAT_VAR_LEVEL"=="worker" {
    summ phat_mu_w phat_sd_w if inwin_3==1
    histogram phat_mu_w if inwin_3==1, width(.02) frequency title("phat mean (worker), in window")
}
else {
    summ phat_mu_wc phat_sd_wc if inwin_3==1
    histogram phat_mu_wc if inwin_3==1, width(.02) frequency title("phat mean (wc), in window")
}
di as text "------------------------------------------------------------"

**********************************************************************
* 3.9) PARAMETRIC phat HETEROGENEITY DiD (obeys global toggles)
**********************************************************************
local mu = cond("$PHAT_VAR_LEVEL"=="worker","phat_mu_w","phat_mu_wc")
local sd = cond("$PHAT_VAR_LEVEL"=="worker","phat_sd_w","phat_sd_wc")

local rhs9 "1.treated_BC_3#1.post_3"
if `did_post_main'==1 local rhs9 "`rhs9' 1.post_3"

* Always include treated×post×phat slopes
local rhs9 "`rhs9' c.`mu'#1.treated_BC_3#1.post_3"
local rhs9 "`rhs9' c.`sd'#1.treated_BC_3#1.post_3"

* Optional: common post×phat slopes
if `did_postx'==1 {
    local rhs9 "`rhs9' c.`mu'#1.post_3"
    local rhs9 "`rhs9' c.`sd'#1.post_3"
}

reghdfe $Y_HET ///
    `rhs9' ///
    $X_3 $LIFECYCLE_3 ///
    if inwin_3==1 & !missing(`mu') & !missing(`sd'), ///
    absorb(`ABS_DID_3') ///
    $VCE_3
estimates store het_phat_param

/**********************************************************************
* 3.10) NONPARAMETRIC phat BINS (mu always; sd optional)
**********************************************************************/
cap drop phat_mu_bin phat_sd_bin

if "$PHAT_VAR_LEVEL"=="worker" {
    xtile phat_mu_bin = phat_mu_w if inwin_3==1, nq($PHAT_BINS_N)
    if $PHAT_SD_BINS==1 xtile phat_sd_bin = phat_sd_w if inwin_3==1, nq($PHAT_BINS_N)

    bysort worker_id: egen byte phat_mu_bin2 = max(phat_mu_bin)
    cap drop phat_mu_bin
    rename phat_mu_bin2 phat_mu_bin

    if $PHAT_SD_BINS==1 {
        bysort worker_id: egen byte phat_sd_bin2 = max(phat_sd_bin)
        cap drop phat_sd_bin
        rename phat_sd_bin2 phat_sd_bin
    }
}
else {
    xtile phat_mu_bin = phat_mu_wc if inwin_3==1, nq($PHAT_BINS_N)
    if $PHAT_SD_BINS==1 xtile phat_sd_bin = phat_sd_wc if inwin_3==1, nq($PHAT_BINS_N)

    bysort wc_id: egen byte phat_mu_bin2 = max(phat_mu_bin)
    cap drop phat_mu_bin
    rename phat_mu_bin2 phat_mu_bin

    if $PHAT_SD_BINS==1 {
        bysort wc_id: egen byte phat_sd_bin2 = max(phat_sd_bin)
        cap drop phat_sd_bin
        rename phat_sd_bin2 phat_sd_bin
    }
}

* RHS for bins model, obey global toggles:
local rhs10 "1.treated_BC_3#1.post_3"
if `did_post_main'==1 local rhs10 "`rhs10' 1.post_3"

* Always: treated×post×bin effects (incremental by bin)
local rhs10 "`rhs10' i.phat_mu_bin#1.treated_BC_3#1.post_3"

* Optional: common post×bin terms
if `did_postx'==1 local rhs10 "`rhs10' i.phat_mu_bin#1.post_3"

reghdfe $Y_HET ///
    `rhs10' ///
    $X_3 $LIFECYCLE_3 ///
    if inwin_3==1 & !missing(phat_mu_bin), ///
    absorb(`ABS_DID_3') ///
    $VCE_3
estimates store het_phat_bins

* Optional sd-bins model
if $PHAT_SD_BINS==1 {
    local rhs10s "1.treated_BC_3#1.post_3"
    if `did_post_main'==1 local rhs10s "`rhs10s' 1.post_3"
    local rhs10s "`rhs10s' i.phat_sd_bin#1.treated_BC_3#1.post_3"
    if `did_postx'==1 local rhs10s "`rhs10s' i.phat_sd_bin#1.post_3"

    reghdfe $Y_HET ///
        `rhs10s' ///
        $X_3 $LIFECYCLE_3 ///
        if inwin_3==1 & !missing(phat_sd_bin), ///
        absorb(`ABS_DID_3') ///
        $VCE_3
    estimates store het_phat_bins_sd
}

**********************************************************************
* 3.11) DIFFICULTY HETEROGENEITY: contract difficulty + weather shock
*   NOTE: if DID_CONTRACTFE_3==1, identification leans on spanning contracts.
**********************************************************************
cap drop est11
gen byte est11 = inwin_3==1 & !missing(contract_diff_z) & !missing(wx_shock)

local rhs11 "1.treated_BC_3#1.post_3"
if `did_post_main'==1 local rhs11 "`rhs11' 1.post_3"

* Always include treated×post×X heterogeneity terms
local rhs11 "`rhs11' c.contract_diff_z#1.treated_BC_3#1.post_3"
local rhs11 "`rhs11' c.wx_shock#1.treated_BC_3#1.post_3"

* Optional common post×X
if `did_postx'==1 {
    local rhs11 "`rhs11' c.contract_diff_z#1.post_3"
    local rhs11 "`rhs11' c.wx_shock#1.post_3"
}

reghdfe $Y_HET ///
    `rhs11' ///
    $X_3 $LIFECYCLE_3 ///
    if est11==1, ///
    absorb(`ABS_DID_3') ///
    $VCE_3
estimates store het_difficulty

**********************************************************************
* 3.12) ABILITY BINS HETEROGENEITY (lower priority; keep)
**********************************************************************
cap drop abil_bin_all
xtile abil_bin_all = ability_z if sample_3==1 & !missing(ability_z), nq(4)

local rhs12 "1.treated_BC_3#1.post_3"
if `did_post_main'==1 local rhs12 "`rhs12' 1.post_3"
local rhs12 "`rhs12' i.abil_bin_all#1.treated_BC_3#1.post_3"
if `did_postx'==1 local rhs12 "`rhs12' i.abil_bin_all#1.post_3"

reghdfe $Y_HET ///
    `rhs12' ///
    $X_3 $LIFECYCLE_3 ///
    if inwin_3==1 & !missing(ability_z) & !missing(abil_bin_all), ///
    absorb(`ABS_DID_3') ///
    $VCE_3
estimates store het_ability_bins

/**********************************************************************
  3.13) DESCRIPTIVE TABLE: ability bins × pre/post (BC only)
**********************************************************************/
table abil_bin_all post_3 if (year==$TY_3) & (province_s=="BC") ///
    & inrange(event_week_3,-$WPRE_3,$WPOST_3), ///
    statistic(mean piece) statistic(mean prod) statistic(mean tot_hrs_day) nformat(%9.3f)

/**********************************************************************
  3.14) QUICK PRINT: key estimates (copy/paste friendly)
**********************************************************************/
cap program drop _trylincom
program define _trylincom
    syntax , COEF(string) LABEL(string)
    capture noisily lincom `coef'
    if _rc==0 di as text "   -> `label'"
    else      di as text "   -> `label' (not in model / omitted / name mismatch)"
end

di as text "============================================================"
di as text "SECTION 3: KEY RESULTS (copy/paste)"
di as text "============================================================"

* (A) Parametric phat
capture quietly estimates restore het_phat_param
if _rc di as error "het_phat_param not found"
else {
    di as text "---- Parametric phat: het_phat_param ----"
    _trylincom, coef("1.treated_BC_3#1.post_3") label("Treat×Post")
    _trylincom, coef("c.`mu'#1.treated_BC_3#1.post_3") label("Treat×Post×phat_mu slope")
    _trylincom, coef("c.`sd'#1.treated_BC_3#1.post_3") label("Treat×Post×phat_sd slope")
}

di as text ""

* (B) Nonparametric phat bins
capture quietly estimates restore het_phat_bins
if _rc di as error "het_phat_bins not found"
else {
    di as text "---- phat mu-bins: het_phat_bins ----"
    _trylincom, coef("1.treated_BC_3#1.post_3") label("Treat×Post (bin 1 baseline)")
    foreach k of numlist 2/$PHAT_BINS_N {
        _trylincom, coef("`k'.phat_mu_bin#1.treated_BC_3#1.post_3") label("Bin `k' increment")
    }
}

di as text ""

* (C) Difficulty
capture quietly estimates restore het_difficulty
if _rc di as error "het_difficulty not found"
else {
    di as text "---- Difficulty: het_difficulty ----"
    _trylincom, coef("1.treated_BC_3#1.post_3") label("Treat×Post")
    _trylincom, coef("c.contract_diff_z#1.treated_BC_3#1.post_3") label("Treat×Post×contract_diff slope")
    _trylincom, coef("c.wx_shock#1.treated_BC_3#1.post_3")        label("Treat×Post×wx_shock slope")
}

di as text ""

* (D) Ability bins
capture quietly estimates restore het_ability_bins
if _rc di as error "het_ability_bins not found"
else {
    di as text "---- Ability bins: het_ability_bins ----"
    _trylincom, coef("1.treated_BC_3#1.post_3") label("Treat×Post (bin 1 baseline)")
    foreach k in 2 3 4 {
        _trylincom, coef("`k'.abil_bin_all#1.treated_BC_3#1.post_3") label("Ability bin `k' increment")
    }
}

di as text "============================================================"


/*
/**********************************************************************
  (3) ADDITIONAL EMPIRICAL ANALYSIS (SELF-CONTAINED)
**********************************************************************/

***********************************************************************
* 3.0) CLEAN RERUN STATE (SAFE)
***********************************************************************
cap estimates clear

foreach v in ///
    sample_3 treated_BC_3 post_3 event_week_3 event_week_3s ///
    ln_prod abil_sample resid_preprod ability_raw ability_z ///
    wx_shock ///
    cd_sample cd_sample_ab cd_sample_bc ///
    resid_cd resid_cd_ab resid_cd_bc ///
    contract_diff_raw contract_diff_z ///
    topup_ind topup_pre phat_topup ///
    phat_mu_w phat_sd_w phat_mu_wc phat_sd_wc ///
    phat_mu_bin phat_sd_bin ///
    abil_bin_all {
    cap drop `v'
}

foreach m in ///
    TY_3 CONTROL_3 Y_HET ///
    CLUSTER_3 VCE_3 ///
    USE_FLOORS_3 MIN_HRS_3 MIN_EARN_3 MIN_PROD_3 USE_MIN_PROD_3 ///
    USE_WEATHER_3 USE_EXPER_3 USE_HOURS_3 USE_DAYFLAGS_3 USE_LIFECYCLE_3 ///
    WPRE_3 WPOST_3 ///
    CD_POST_MODE CD_TIMEFE ///
    DID_WEEKFE_3 ///
    PHAT_VAR_LEVEL PHAT_BINS_N PHAT_SD_BINS ///
    X_3 LIFECYCLE_3 SAMPLE_OK_3 {
    cap macro drop `m'
}

***********************************************************************
* 3.1) SETTINGS (EDIT THESE)
***********************************************************************
global TY_3      2018
global CONTROL_3 "AB"
global Y_HET     "ln_incentive"

global CLUSTER_3 "worker_id"
global VCE_3     "vce(cluster $CLUSTER_3)"

global USE_WEATHER_3   1
global USE_EXPER_3     1
global USE_HOURS_3     0
global USE_DAYFLAGS_3  0
global USE_LIFECYCLE_3 1

global USE_FLOORS_3   0
global MIN_HRS_3      2
global MIN_EARN_3     30
global USE_MIN_PROD_3 1
global MIN_PROD_3     50

global WPRE_3  3
global WPOST_3 8

* Contract difficulty timing:
* 0 = PRE only (AB+BC)
* 1 = AB pre+post, BC pre only   (recommended)
* 2 = POST allowed for all (AB+BC)
global CD_POST_MODE 1

* Week FE in difficulty construction
global CD_TIMEFE 0

* Week FE in ALL DiD/heterogeneity regressions
global DID_WEEKFE_3 1

* phat moments level
global PHAT_VAR_LEVEL "wc"

* phat bins
global PHAT_BINS_N 4
global PHAT_SD_BINS 0

***********************************************************************
* 3.2) BUILD CONTROLS + FLOORS STRINGS
***********************************************************************
global X_3 ""
if $USE_HOURS_3==1    global X_3 "$X_3 tot_hrs_day"
if $USE_DAYFLAGS_3==1 global X_3 "$X_3 multi_contract_day"
if $USE_EXPER_3==1    global X_3 "$X_3 c.experience##c.experience"
if $USE_WEATHER_3==1  global X_3 "$X_3 c.temp##c.temp c.wind_spd##c.wind_spd precip0 c.precip_pos##c.precip_pos"

global LIFECYCLE_3 ""
if $USE_LIFECYCLE_3==1 {
    global LIFECYCLE_3 "$LIFECYCLE_3 c.contract_day##c.contract_day"
    global LIFECYCLE_3 "$LIFECYCLE_3 c.season_days_worked##c.season_days_worked"
    global LIFECYCLE_3 "$LIFECYCLE_3 c.season_days_worked##c.cum_seasons"
}

global SAMPLE_OK_3 "1==1"
if $USE_FLOORS_3==1 {
    global SAMPLE_OK_3 "$SAMPLE_OK_3 & tot_hrs_day >= $MIN_HRS_3"
    global SAMPLE_OK_3 "$SAMPLE_OK_3 & piece > 0 & prod > 0"
    global SAMPLE_OK_3 "$SAMPLE_OK_3 & (piece*prod) >= $MIN_EARN_3"
    if $USE_MIN_PROD_3==1 global SAMPLE_OK_3 "$SAMPLE_OK_3 & prod >= $MIN_PROD_3"
}

di as text "------------------------------------------------------------"
di as text "[3] TY_3:           $TY_3   CONTROL_3: $CONTROL_3   Y_HET: $Y_HET"
di as text "[3] Window:         [-$WPRE_3,$WPOST_3]"
di as text "[3] CD_POST_MODE:   $CD_POST_MODE   CD_TIMEFE: $CD_TIMEFE   DID_WEEKFE_3: $DID_WEEKFE_3"
di as text "[3] PHAT_VAR_LEVEL: $PHAT_VAR_LEVEL   PHAT_BINS_N: $PHAT_BINS_N   PHAT_SD_BINS: $PHAT_SD_BINS"
di as text "------------------------------------------------------------"

***********************************************************************
* 3.3) DEFINE SAMPLE + EVENT TIME
***********************************************************************
gen byte sample_3 = (year==$TY_3) & inlist(province_s,"BC","$CONTROL_3")
gen byte treated_BC_3 = .
replace treated_BC_3 = (province_s=="BC") if sample_3==1

local tdate_3 = mdy(6,1,$TY_3)
gen byte post_3       = (date >= `tdate_3') if sample_3==1
gen int  event_week_3 = floor((date - `tdate_3')/7) if sample_3==1
gen int  event_week_3s = event_week_3 + 100 if sample_3==1

tab province_s if sample_3==1
tab post_3     if sample_3==1

**********************************************************************
* 3.4) PRE-TREATMENT ABILITY (BC+AB, PRE ONLY)
**********************************************************************
gen double ln_prod = ln(prod) if prod>0

gen byte abil_sample = (sample_3==1) & (post_3==0) & !missing(ln_prod) & ($SAMPLE_OK_3)

reghdfe ln_prod $X_3 $LIFECYCLE_3 if abil_sample==1, absorb(contract_id) resid(resid_preprod)

bysort worker_id: egen double ability_raw = mean(resid_preprod)
quietly summarize ability_raw if !missing(ability_raw)
scalar abil_mu = r(mean)
scalar abil_sd = r(sd)

gen double ability_z = (ability_raw - abil_mu)/abil_sd if !missing(ability_raw) & abil_sd>0

**********************************************************************
* 3.5) DAILY WEATHER SHOCK (reduced-form)
**********************************************************************
quietly reg ln_prod c.temp##c.temp c.wind_spd##c.wind_spd precip0 c.precip_pos##c.precip_pos ///
    if sample_3==1 & !missing(ln_prod), vce(robust)
predict double wx_shock if sample_3==1, xb

/**********************************************************************
* 3.6) CONTRACT DIFFICULTY (AB + BC, consistent CD_POST_MODE)
**********************************************************************/
cap drop cd_sample_ab cd_sample_bc cd_sample resid_cd_ab resid_cd_bc resid_cd contract_diff_raw contract_diff_z

gen byte cd_sample_ab = (year==$TY_3) & (province_s=="AB") & !missing(ln_prod) & !missing(ability_z) & ($SAMPLE_OK_3)
gen byte cd_sample_bc = (year==$TY_3) & (province_s=="BC") & !missing(ln_prod) & !missing(ability_z) & ($SAMPLE_OK_3)

* Apply CD_POST_MODE (robust: copy global into local first)
local cdmode = $CD_POST_MODE

if `cdmode'==0 {
    replace cd_sample_ab = cd_sample_ab & (post_3==0)
    replace cd_sample_bc = cd_sample_bc & (post_3==0)
}
else if `cdmode'==1 {
    * AB pre+post allowed; BC pre only
    replace cd_sample_bc = cd_sample_bc & (post_3==0)
}
else if `cdmode'==2 {
    * allow post for all: do nothing
}

gen byte cd_sample = cd_sample_ab | cd_sample_bc

di as text "[3.6] cd_sample_ab N=" _N " (see counts below)"
count if cd_sample_ab==1
count if cd_sample_bc==1

* AB residualization
local abs_cd "worker_id"
if $CD_TIMEFE==1 local abs_cd "`abs_cd' event_week_3s"

reghdfe ln_prod c.ability_z $X_3 $LIFECYCLE_3 if cd_sample_ab==1, absorb(`abs_cd') resid(resid_cd_ab)
reghdfe ln_prod c.ability_z $X_3 $LIFECYCLE_3 if cd_sample_bc==1, absorb(`abs_cd') resid(resid_cd_bc)

gen double resid_cd = .
replace resid_cd = resid_cd_ab if !missing(resid_cd_ab)
replace resid_cd = resid_cd_bc if !missing(resid_cd_bc)

bysort contract_id: egen double contract_diff_raw = mean(resid_cd)

quietly summarize contract_diff_raw if !missing(contract_diff_raw)
scalar cd_mu = r(mean)
scalar cd_sd = r(sd)

gen double contract_diff_z = (contract_diff_raw - cd_mu)/cd_sd if !missing(contract_diff_raw) & cd_sd>0

di as text "[3.6] Contract difficulty coverage in sample_3:"
count if sample_3==1 & province_s=="AB" & !missing(contract_diff_z)
count if sample_3==1 & province_s=="BC" & !missing(contract_diff_z)

**********************************************************************
* 3.7) PREDICTED TOP-UP RISK phat_topup (PRE only; predict all)
**********************************************************************
gen byte topup_ind = (top_up > 0) if !missing(top_up)

gen byte topup_pre = (sample_3==1) & (post_3==0) ///
    & !missing(topup_ind) & !missing(ability_z) & !missing(mw) & !missing(wx_shock) ///
    & ($SAMPLE_OK_3)

reg topup_ind c.mw c.ability_z c.wx_shock $X_3 i.contract_id if topup_pre==1, vce(cluster $CLUSTER_3)
predict double phat_topup if sample_3==1, xb
replace phat_topup = max(0, min(1, phat_topup)) if sample_3==1 & !missing(phat_topup)

**********************************************************************
* 3.8) PHAT MOMENTS (within analysis window; worker or wc)
**********************************************************************
tempvar inwin
gen byte `inwin' = (sample_3==1) & inrange(event_week_3,-$WPRE_3,$WPOST_3) ///
    & ($SAMPLE_OK_3) & !missing($Y_HET) & !missing(phat_topup)

cap drop phat_mu_w phat_sd_w phat_mu_wc phat_sd_wc

if "$PHAT_VAR_LEVEL"=="worker" {
    bysort worker_id: egen double phat_mu_w = mean(cond(`inwin'==1, phat_topup, .))
    bysort worker_id: egen double phat_sd_w = sd(  cond(`inwin'==1, phat_topup, .))
}
else {
    capture confirm variable wc_id
    if _rc {
        di as error "[3.8] wc_id not found -> falling back to worker moments"
        bysort worker_id: egen double phat_mu_w = mean(cond(`inwin'==1, phat_topup, .))
        bysort worker_id: egen double phat_sd_w = sd(  cond(`inwin'==1, phat_topup, .))
        global PHAT_VAR_LEVEL "worker"
    }
    else {
        bysort wc_id: egen double phat_mu_wc = mean(cond(`inwin'==1, phat_topup, .))
        bysort wc_id: egen double phat_sd_wc = sd(  cond(`inwin'==1, phat_topup, .))
    }
}

**********************************************************************
* 3.9) PARAMETRIC phat HETEROGENEITY DiD (with optional week FE)
**********************************************************************
local abs_did "worker_id contract_id"
if $DID_WEEKFE_3==1 local abs_did "`abs_did' event_week_3s"

if "$PHAT_VAR_LEVEL"=="worker" {
    reghdfe $Y_HET ///
        i.treated_BC_3##i.post_3 ///
        c.phat_mu_w##i.treated_BC_3##i.post_3 ///
        c.phat_sd_w##i.treated_BC_3##i.post_3 ///
        $X_3 $LIFECYCLE_3 ///
        if `inwin'==1 & !missing(phat_mu_w) & !missing(phat_sd_w), ///
        absorb(`abs_did') ///
        $VCE_3
}
else {
    reghdfe $Y_HET ///
        i.treated_BC_3##i.post_3 ///
        c.phat_mu_wc##i.treated_BC_3##i.post_3 ///
        c.phat_sd_wc##i.treated_BC_3##i.post_3 ///
        $X_3 $LIFECYCLE_3 ///
        if `inwin'==1 & !missing(phat_mu_wc) & !missing(phat_sd_wc), ///
        absorb(`abs_did') ///
        $VCE_3
}
estimates store het_phat_param

/**********************************************************************
* 3.10) NONPARAMETRIC phat BINS (mu always; sd optional)
**********************************************************************/
cap drop phat_mu_bin phat_sd_bin

if "$PHAT_VAR_LEVEL"=="worker" {
    xtile phat_mu_bin = phat_mu_w if `inwin'==1, nq($PHAT_BINS_N)
    if $PHAT_SD_BINS==1 xtile phat_sd_bin = phat_sd_w if `inwin'==1, nq($PHAT_BINS_N)
}
else {
    xtile phat_mu_bin = phat_mu_wc if `inwin'==1, nq($PHAT_BINS_N)
    if $PHAT_SD_BINS==1 xtile phat_sd_bin = phat_sd_wc if `inwin'==1, nq($PHAT_BINS_N)
}

* Propagate bins to all rows within unit (so DiD estimation uses full inwin sample cleanly)
if "$PHAT_VAR_LEVEL"=="worker" {
    bysort worker_id: egen byte phat_mu_bin2 = max(phat_mu_bin)
    cap drop phat_mu_bin
    rename phat_mu_bin2 phat_mu_bin
    if $PHAT_SD_BINS==1 {
        bysort worker_id: egen byte phat_sd_bin2 = max(phat_sd_bin)
        cap drop phat_sd_bin
        rename phat_sd_bin2 phat_sd_bin
    }
}
else {
    bysort wc_id: egen byte phat_mu_bin2 = max(phat_mu_bin)
    cap drop phat_mu_bin
    rename phat_mu_bin2 phat_mu_bin
    if $PHAT_SD_BINS==1 {
        bysort wc_id: egen byte phat_sd_bin2 = max(phat_sd_bin)
        cap drop phat_sd_bin
        rename phat_sd_bin2 phat_sd_bin
    }
}

* mu bins regression
reghdfe $Y_HET ///
    i.treated_BC_3##i.post_3 ///
    i.phat_mu_bin##i.treated_BC_3##i.post_3 ///
    $X_3 $LIFECYCLE_3 ///
    if `inwin'==1 & !missing(phat_mu_bin), ///
    absorb(`abs_did') ///
    $VCE_3
estimates store het_phat_bins

* optional sd-bins regression
if $PHAT_SD_BINS==1 {
    reghdfe $Y_HET ///
        i.treated_BC_3##i.post_3 ///
        i.phat_sd_bin##i.treated_BC_3##i.post_3 ///
        $X_3 $LIFECYCLE_3 ///
        if `inwin'==1 & !missing(phat_sd_bin), ///
        absorb(`abs_did') ///
        $VCE_3
    estimates store het_phat_bins_sd
}


**********************************************************************
* 3.11) DIFFICULTY HETEROGENEITY: contract difficulty + weather shock
**********************************************************************

* replicate the exact estimation sample for 3.11
gen byte est11 = sample_3==1 ///
    & inrange(event_week_3,-$WPRE_3,$WPOST_3) ///
    & ($SAMPLE_OK_3) ///
    & !missing($Y_HET) ///
    & !missing(contract_diff_z) ///
    & !missing(wx_shock)

* does each contract have both pre and post within the window?
bys contract_id: egen byte minpost11 = min(post_3) if est11
bys contract_id: egen byte maxpost11 = max(post_3) if est11
gen byte spans11 = (minpost11==0 & maxpost11==1) if est11

* count spanning contracts by province
bys contract_id: egen byte any11 = max(est11)
bys contract_id: egen byte provBC = max(treated_BC_3) if any11

count if any11==1 & provBC==0 & spans11==1   // AB spanning contracts
count if any11==1 & provBC==1 & spans11==1   // BC spanning contracts

* optional: how many total contracts in each group?
count if any11==1 & provBC==0
count if any11==1 & provBC==1


reghdfe $Y_HET ///
    i.treated_BC_3##i.post_3 ///
    c.contract_diff_z#i.treated_BC_3#i.post_3 ///
    c.wx_shock#i.treated_BC_3#i.post_3 ///
    $X_3 $LIFECYCLE_3 ///
    if sample_3==1 ///
       & inrange(event_week_3,-$WPRE_3,$WPOST_3) ///
       & ($SAMPLE_OK_3) ///
       & !missing($Y_HET) ///
       & !missing(contract_diff_z) ///
       & !missing(wx_shock), ///
    absorb(worker_id contract_id `ABS_TIME_3') ///
    $VCE_3
estimates store het_difficulty

**********************************************************************
* 3.12) ABILITY BINS HETEROGENEITY (lower priority; keep)
**********************************************************************
cap drop abil_bin_all
xtile abil_bin_all = ability_z if sample_3==1 & !missing(ability_z), nq(4)

reghdfe $Y_HET ///
    i.treated_BC_3##i.post_3 ///
    i.abil_bin_all##i.treated_BC_3##i.post_3 ///
    $X_3 $LIFECYCLE_3 ///
    if sample_3==1 ///
       & inrange(event_week_3,-$WPRE_3,$WPOST_3) ///
       & ($SAMPLE_OK_3) ///
       & !missing($Y_HET) ///
       & !missing(ability_z), ///
    absorb(worker_id contract_id `ABS_TIME_3') ///
    $VCE_3
estimates store het_ability_bins

/**********************************************************************
  3.13) DESCRIPTIVE TABLE: ability bins × pre/post (BC only)
**********************************************************************/
table abil_bin_all post_3 if (year==$TY_3) & (province_s=="BC") ///
    & inrange(event_week_3,-$WPRE_3,$WPOST_3), ///
    statistic(mean piece) statistic(mean prod) statistic(mean tot_hrs_day) nformat(%9.3f)

**********************************************************************
  3.14) QUICK CHECK + QUICK PRINT (copy/paste friendly)
  - This is intentionally simple and robust:
    (i) restore the model
    (ii) attempt a small set of lincom calls
**********************************************************************/

cap program drop _trylincom
program define _trylincom
    syntax , COEF(string) LABEL(string)
    capture noisily lincom `coef'
    if _rc==0 di as text "   -> `label'"
    else      di as text "   -> `label' (not in model / omitted / name mismatch)"
end

di as text "============================================================"
di as text "SECTION 3: KEY RESULTS (copy/paste)"
di as text "============================================================"

* -----------------------------
* (A) Parametric phat model
* -----------------------------
capture quietly estimates restore het_phat_param
if _rc di as error "het_phat_param not found"
else {
    di as text "---- Parametric phat heterogeneity: het_phat_param ----"
    _trylincom, coef("1.treated_BC_3#1.post_3") label("Baseline DiD")

    * slopes: try wc first, then worker (only one will exist)
    _trylincom, coef("1.treated_BC_3#1.post_3#c.phat_mu_wc") label("DiD heterogeneity wrt mean phat (wc)")
    _trylincom, coef("1.treated_BC_3#1.post_3#c.phat_mu_w")  label("DiD heterogeneity wrt mean phat (worker)")

    _trylincom, coef("1.treated_BC_3#1.post_3#c.phat_sd_wc") label("DiD heterogeneity wrt sd phat (wc)")
    _trylincom, coef("1.treated_BC_3#1.post_3#c.phat_sd_w")  label("DiD heterogeneity wrt sd phat (worker)")

    di as text "   (sanity) moments:"
    capture noisily summ phat_mu_wc phat_sd_wc if sample_3==1
    capture noisily summ phat_mu_w  phat_sd_w  if sample_3==1
}

di as text ""

* -----------------------------
* (B) Nonparametric phat bins (mu) model
* -----------------------------
capture quietly estimates restore het_phat_bins
if _rc di as error "het_phat_bins not found"
else {
    di as text "---- Nonparametric phat mu-bins: het_phat_bins ----"
    _trylincom, coef("1.treated_BC_3#1.post_3") label("Baseline DiD (omitted mu-bin)")
    foreach k of numlist 2/6 {
        _trylincom, coef("`k'.phat_mu_bin#1.treated_BC_3#1.post_3") label("Mu-bin `k' increment")
    }
}

di as text ""

* Optional sd-bin model print
if $PHAT_SD_BINS==1 {
    capture quietly estimates restore het_phat_sdbins
    if _rc di as error "het_phat_sdbins not found"
    else {
        di as text "---- Nonparametric phat sd-bins: het_phat_sdbins ----"
        _trylincom, coef("1.treated_BC_3#1.post_3") label("Baseline DiD (omitted sd-bin)")
        foreach k of numlist 2/6 {
            _trylincom, coef("`k'.phat_sd_bin#1.treated_BC_3#1.post_3") label("SD-bin `k' increment")
        }
    }
    di as text ""
}

* -----------------------------
* (C) Difficulty model
* -----------------------------
capture quietly estimates restore het_difficulty
if _rc di as error "het_difficulty not found"
else {
    di as text "---- Difficulty heterogeneity: het_difficulty ----"
    _trylincom, coef("1.treated_BC_3#1.post_3") label("Baseline DiD")
    _trylincom, coef("1.treated_BC_3#1.post_3#c.contract_diff_z") label("DiD heterogeneity wrt contract difficulty")
    _trylincom, coef("1.treated_BC_3#1.post_3#c.wx_shock")        label("DiD heterogeneity wrt weather shock")
}

di as text ""

* -----------------------------
* (D) Ability bins model (lower priority)
* -----------------------------
capture quietly estimates restore het_ability_bins
if _rc di as error "het_ability_bins not found"
else {
    di as text "---- Ability-bin heterogeneity: het_ability_bins ----"
    _trylincom, coef("1.treated_BC_3#1.post_3") label("Baseline DiD (omitted ability bin)")
    foreach k in 2 3 4 {
        _trylincom, coef("`k'.abil_bin_all#1.treated_BC_3#1.post_3") label("Ability-bin `k' increment")
    }
}

di as text "============================================================"






















/*

/**********************************************************************
  (3) ADDITIONAL EMPIRICAL ANALYSIS (SELF-CONTAINED)
  Focus: BC vs AB in treated year (2A-style), heterogeneity by:
    (i) predicted top-up risk (MW-binding likelihood) + variance
        - parametric: phat_mu + phat_sd
        - nonparametric: bins (quartiles) of phat_mu (and optional phat_sd)
    (ii) job difficulty: contract-level + daily weather difficulty
    (iii) ability bins (kept, but lower emphasis)

  KEY DESIGN CHOICES:
    - Ability: estimated PRE only (post_3==0), propagated to all rows by worker.
    - Contract difficulty: estimated using AB-only data; toggle AB post inclusion.
    - phat_topup: time-varying by construction (mw + weather), model estimated PRE only.
    - phat moments: computed within analysis window at worker or wc_id level.
    - DiD time FE: toggle event-week FE in all DiD/heterogeneity regressions.
**********************************************************************/

***********************************************************************
* 3.0) CLEAN RERUN STATE (SAFE)
***********************************************************************
cap estimates clear

foreach v in ///
    sample_3 treated_BC_3 post_3 event_week_3 event_week_3s ///
    ln_prod abil_sample resid_preprod ability_raw ability_z ///
    topup_ind topup_pre phat_topup ///
    cd_sample resid_cd contract_diff_raw contract_diff_z ///
    wx_shock ///
    phat_mu_w phat_sd_w phat_mu_wc phat_sd_wc ///
    phat_mu_bin phat_sd_bin ///
    abil_bin_all {
    cap drop `v'
}

foreach m in ///
    TY_3 CONTROL_3 Y_HET CLUSTER_3 VCE_3 ///
    USE_FLOORS_3 MIN_HRS_3 MIN_EARN_3 MIN_PROD_3 USE_MIN_PROD_3 ///
    USE_WEATHER_3 USE_EXPER_3 USE_HOURS_3 USE_DAYFLAGS_3 USE_LIFECYCLE_3 ///
    WPRE_3 WPOST_3 ///
    CD_USE_AB_POST CD_TIMEFE ///
    DID_WEEKFE_3 ///
    PHAT_VAR_LEVEL ///
    PHAT_BINS_N ///
    PHAT_SD_BINS ///
    X_3 LIFECYCLE_3 SAMPLE_OK_3 {
    cap macro drop `m'
}

***********************************************************************
* 3.1) SETTINGS (EDIT THESE)
***********************************************************************
* Treated year and control province
global TY_3      2018
global CONTROL_3 "AB"

* Outcome for heterogeneity regressions
global Y_HET "lprod"

* Clustering
global CLUSTER_3 "worker_id"
global VCE_3 "vce(cluster $CLUSTER_3)"

* Controls toggles (mirror Section 2 logic)
global USE_WEATHER_3   1
global USE_EXPER_3     1
global USE_HOURS_3     0
global USE_DAYFLAGS_3  0
global USE_LIFECYCLE_3 1

* Floors toggle (0 = none)
global USE_FLOORS_3   0
global MIN_HRS_3      2
global MIN_EARN_3     30
global USE_MIN_PROD_3 1
global MIN_PROD_3     50

* Event window
global WPRE_3  3
global WPOST_3 8

* ---------------------------
* KEY TOGGLES
* ---------------------------

* (A) Contract difficulty estimation sample (AB-only):
*   0 = AB pre-only (conservative)
*   1 = AB pre+post (more support; assumes no spillovers into AB)
global CD_USE_AB_POST 1

* (A2) Add common time effects in AB-only difficulty estimation:
*   0 = none beyond controls
*   1 = include event-week FE (recommended)
global CD_TIMEFE 0

* (B) Add week fixed effects to ALL DiD/heterogeneity regressions:
*   0 = no week FE
*   1 = absorb(event_week_3s) in all DiD models
global DID_WEEKFE_3 1

* (C) Binding-risk moment level:
*   "worker" = moments by worker_id
*   "wc"     = moments by wc_id (worker×contract)
global PHAT_VAR_LEVEL "wc"

* (D) Nonparametric phat bins
global PHAT_BINS_N 4

* (E) Also run bins for phat_sd? (0/1)
global PHAT_SD_BINS 0


***********************************************************************
* 3.2) BUILD CONTROLS + FLOORS STRINGS
***********************************************************************
global X_3 ""
if $USE_HOURS_3==1    global X_3 "$X_3 tot_hrs_day"
if $USE_DAYFLAGS_3==1 global X_3 "$X_3 multi_contract_day"
if $USE_EXPER_3==1    global X_3 "$X_3 c.experience##c.experience"
if $USE_WEATHER_3==1  global X_3 "$X_3 c.temp##c.temp c.wind_spd##c.wind_spd precip0 c.precip_pos##c.precip_pos"

global LIFECYCLE_3 ""
if $USE_LIFECYCLE_3==1 {
    global LIFECYCLE_3 "$LIFECYCLE_3 c.contract_day##c.contract_day"
    global LIFECYCLE_3 "$LIFECYCLE_3 c.season_days_worked##c.season_days_worked"
    global LIFECYCLE_3 "$LIFECYCLE_3 c.season_days_worked##c.cum_seasons"
}

global SAMPLE_OK_3 "1==1"
if $USE_FLOORS_3==1 {
    global SAMPLE_OK_3 "$SAMPLE_OK_3 & tot_hrs_day >= $MIN_HRS_3"
    global SAMPLE_OK_3 "$SAMPLE_OK_3 & piece > 0 & prod > 0"
    global SAMPLE_OK_3 "$SAMPLE_OK_3 & (piece*prod) >= $MIN_EARN_3"
    if $USE_MIN_PROD_3==1 global SAMPLE_OK_3 "$SAMPLE_OK_3 & prod >= $MIN_PROD_3"
}

di as text "------------------------------------------------------------"
di as text "[3] TY_3:            $TY_3"
di as text "[3] CONTROL_3:       $CONTROL_3"
di as text "[3] Y_HET:           $Y_HET"
di as text "[3] CLUSTER:         $CLUSTER_3"
di as text "[3] Window:          [-$WPRE_3,$WPOST_3]"
di as text "[3] CD_USE_AB_POST:  $CD_USE_AB_POST"
di as text "[3] CD_TIMEFE:       $CD_TIMEFE"
di as text "[3] DID_WEEKFE_3:    $DID_WEEKFE_3"
di as text "[3] PHAT_VAR_LEVEL:  $PHAT_VAR_LEVEL"
di as text "[3] PHAT_BINS_N:     $PHAT_BINS_N"
di as text "[3] PHAT_SD_BINS:    $PHAT_SD_BINS"
di as text "------------------------------------------------------------"


***********************************************************************
* 3.3) DEFINE 2A-LIKE SAMPLE (BC vs AB in treated year) + EVENT TIME
***********************************************************************
gen byte sample_3 = (year==$TY_3) & inlist(province_s,"BC","$CONTROL_3")

gen byte treated_BC_3 = .
replace treated_BC_3 = (province_s=="BC") if sample_3==1

local tdate_3 = mdy(6,1,$TY_3)
gen byte post_3       = (date >= `tdate_3') if sample_3==1
gen int  event_week_3 = floor((date - `tdate_3')/7) if sample_3==1

* Stata factor vars cannot be negative: shift week index for FE use.
gen int event_week_3s = event_week_3 + 100 if sample_3==1

* Common time FE absorb block for DiD models
local ABS_TIME_3 ""
if $DID_WEEKFE_3==1 local ABS_TIME_3 "event_week_3s"

tab province_s if sample_3==1
tab post_3     if sample_3==1


/**********************************************************************
  3.4) PRE-TREATMENT ABILITY (BC+AB, PRE ONLY)
**********************************************************************/
cap drop ln_prod
gen double ln_prod = ln(prod) if prod>0

cap drop abil_sample
gen byte abil_sample = (sample_3==1) ///
    & (post_3==0) ///
    & !missing(ln_prod) ///
    & ($SAMPLE_OK_3)

di as text "[3.4] Ability PRE sample:"
count if abil_sample==1
tab province_s if abil_sample==1

cap drop resid_preprod
reghdfe ln_prod $X_3 $LIFECYCLE_3 if abil_sample==1, ///
    absorb(contract_id) ///
    resid(resid_preprod)

cap drop ability_raw ability_z
bysort worker_id: egen double ability_raw = mean(resid_preprod)

quietly summarize ability_raw if !missing(ability_raw)
scalar abil_mu = r(mean)
scalar abil_sd = r(sd)

gen double ability_z = (ability_raw - abil_mu)/abil_sd ///
    if !missing(ability_raw) & abil_sd>0

di as text "[3.4] Ability coverage:"
count if sample_3==1 & !missing(ability_z)
tab province_s if sample_3==1 & !missing(ability_z)


/**********************************************************************
  3.5) DAILY WEATHER SHOCK (not MW-contaminated)
**********************************************************************/
cap drop wx_shock
quietly reg ln_prod ///
    c.temp##c.temp c.wind_spd##c.wind_spd precip0 c.precip_pos##c.precip_pos ///
    if sample_3==1 & !missing(ln_prod), vce(robust)
predict double wx_shock if sample_3==1, xb
summ wx_shock if sample_3==1


/**********************************************************************
  3.6) CONTRACT DIFFICULTY (AB-only; toggle pre-only vs pre+post; optional week FE)
**********************************************************************/
cap drop cd_sample resid_cd contract_diff_raw contract_diff_z

gen byte cd_sample = (year==$TY_3) & (province_s=="AB") ///
    & !missing(ln_prod) ///
    & !missing(ability_z) ///
    & ($SAMPLE_OK_3)

if $CD_USE_AB_POST==0 replace cd_sample = cd_sample & (post_3==0)

di as text "[3.6] Contract difficulty sample (AB-only):"
count if cd_sample==1
tab post_3 if cd_sample==1

if $CD_TIMEFE==1 {
    reghdfe ln_prod ///
        c.ability_z ///
        $X_3 ///
        $LIFECYCLE_3 ///
        if cd_sample==1, ///
        absorb(worker_id event_week_3s) ///
        resid(resid_cd)
}
else {
    reghdfe ln_prod ///
        c.ability_z ///
        $X_3 ///
        $LIFECYCLE_3 ///
        if cd_sample==1, ///
        absorb(worker_id) ///
        resid(resid_cd)
}

bysort contract_id: egen double contract_diff_raw = mean(resid_cd)

quietly summarize contract_diff_raw if !missing(contract_diff_raw)
scalar cd_mu = r(mean)
scalar cd_sd = r(sd)

gen double contract_diff_z = (contract_diff_raw - cd_mu)/cd_sd ///
    if !missing(contract_diff_raw) & cd_sd>0

di as text "[3.6] Contract difficulty coverage (in sample_3):"
count if sample_3==1 & province_s=="AB" & !missing(contract_diff_z)
count if sample_3==1 & province_s=="BC" & !missing(contract_diff_z)


/**********************************************************************
  3.7) PREDICTED TOP-UP RISK phat_topup (time-varying; uses mw + wx_shock)
    - Estimate PRE only
    - Predict for ALL sample_3 days
**********************************************************************/
cap drop topup_ind topup_pre phat_topup
gen byte topup_ind = (top_up > 0) if !missing(top_up)

gen byte topup_pre = (sample_3==1) ///
    & (post_3==0) ///
    & !missing(topup_ind) ///
    & !missing(ability_z) ///
    & !missing(mw) ///
    & !missing(wx_shock) ///
    & ($SAMPLE_OK_3)

di as text "[3.7] Top-up PRE estimation sample:"
count if topup_pre==1
tab province_s if topup_pre==1
summ topup_ind if topup_pre==1

reg topup_ind ///
    c.mw ///
    c.ability_z ///
    c.wx_shock ///
    $X_3 ///
    i.contract_id ///
    if topup_pre==1, ///
    vce(cluster $CLUSTER_3)

predict double phat_topup if sample_3==1, xb
replace phat_topup = max(0, min(1, phat_topup)) if sample_3==1 & !missing(phat_topup)

di as text "[3.7] phat_topup coverage:"
summ phat_topup if sample_3==1
tab post_3 if sample_3==1, summarize(phat_topup)


/**********************************************************************
  3.8) BINDING-RISK MOMENTS (mean + sd) WITHIN ANALYSIS WINDOW
**********************************************************************/
cap drop phat_mu_w phat_sd_w phat_mu_wc phat_sd_wc

tempvar inwin
gen byte `inwin' = (sample_3==1) ///
    & inrange(event_week_3,-$WPRE_3,$WPOST_3) ///
    & ($SAMPLE_OK_3) ///
    & !missing($Y_HET) ///
    & !missing(phat_topup)

if "$PHAT_VAR_LEVEL"=="worker" {
    bysort worker_id: egen double phat_mu_w = mean(cond(`inwin'==1, phat_topup, .))
    bysort worker_id: egen double phat_sd_w = sd(  cond(`inwin'==1, phat_topup, .))
}
else {
    capture confirm variable wc_id
    if _rc {
        di as error "[3.8] wc_id not found; falling back to worker-level moments."
        global PHAT_VAR_LEVEL "worker"
        bysort worker_id: egen double phat_mu_w = mean(cond(`inwin'==1, phat_topup, .))
        bysort worker_id: egen double phat_sd_w = sd(  cond(`inwin'==1, phat_topup, .))
    }
    else {
        bysort wc_id: egen double phat_mu_wc = mean(cond(`inwin'==1, phat_topup, .))
        bysort wc_id: egen double phat_sd_wc = sd(  cond(`inwin'==1, phat_topup, .))
    }
}

di as text "[3.8] phat moment diagnostics:"
count if `inwin'==1
summ phat_topup if `inwin'==1
if "$PHAT_VAR_LEVEL"=="worker" summ phat_mu_w phat_sd_w if sample_3==1
else                          summ phat_mu_wc phat_sd_wc if sample_3==1


/**********************************************************************
  3.9) MAIN HETEROGENEITY DiD: MW effect vs binding-risk MEAN + VOLATILITY
  - Adds optional week FE (toggle via USE_WEEKFE_3) through `ABS_TIME_3'
  - Stores: het_phat_param
**********************************************************************/

* Guard: ensure SAMPLE_OK_3 is never empty (prevents & () parse errors)
capture confirm global SAMPLE_OK_3
if _rc | "$SAMPLE_OK_3"=="" global SAMPLE_OK_3 "1==1"

* (Optional) time FE absorbed in the DiD/heterogeneity regressions:
*   - reghdfe factor vars can't use negative event times, so we use event_week_3s.
*   - Set earlier in settings: global USE_WEEKFE_3 0/1
local ABS_TIME_3 ""
local ABS_TIME_3 ""
capture confirm global USE_WEEKFE_3
if !_rc {
    if $USE_WEEKFE_3==1 local ABS_TIME_3 "event_week_3s"
}

* Run parametric phat heterogeneity depending on moment level
if "$PHAT_VAR_LEVEL"=="worker" {

    reghdfe $Y_HET ///
        i.treated_BC_3##i.post_3 ///
        c.phat_mu_w##i.treated_BC_3##i.post_3 ///
        c.phat_sd_w##i.treated_BC_3##i.post_3 ///
        $X_3 $LIFECYCLE_3 ///
        if sample_3==1 ///
           & inrange(event_week_3,-$WPRE_3,$WPOST_3) ///
           & $SAMPLE_OK_3 ///
           & !missing($Y_HET) ///
           & !missing(phat_mu_w) ///
           & !missing(phat_sd_w), ///
        absorb(worker_id contract_id `ABS_TIME_3') ///
        $VCE_3

    estimates store het_phat_param
}
else {

    reghdfe $Y_HET ///
        i.treated_BC_3##i.post_3 ///
        c.phat_mu_wc##i.treated_BC_3##i.post_3 ///
        c.phat_sd_wc##i.treated_BC_3##i.post_3 ///
        $X_3 $LIFECYCLE_3 ///
        if sample_3==1 ///
           & inrange(event_week_3,-$WPRE_3,$WPOST_3) ///
           & $SAMPLE_OK_3 ///
           & !missing($Y_HET) ///
           & !missing(phat_mu_wc) ///
           & !missing(phat_sd_wc), ///
        absorb(worker_id contract_id `ABS_TIME_3') ///
        $VCE_3

    estimates store het_phat_param
}



/**********************************************************************
  3.10) PHAT HETEROGENEITY (NONPARAMETRIC BINS)
  - Bins for mean binding risk (phat_mu_*) and optionally volatility (phat_sd_*)
  - Uses *analysis window* only
**********************************************************************/

* --- settings (safe defaults)
capture confirm global PHAT_BINS_N
if _rc global PHAT_BINS_N 4   // default quartiles

* --- define window indicator locally (DO NOT rely on earlier tempvars)
tempvar inwin10
gen byte `inwin10' = (sample_3==1) ///
    & inrange(event_week_3,-$WPRE_3,$WPOST_3) ///
    & ($SAMPLE_OK_3) ///
    & !missing($Y_HET)

* --- choose which phat moments to bin (worker vs wc)
cap drop phat_mu_bin phat_sd_bin

if "$PHAT_VAR_LEVEL"=="worker" {

    xtile phat_mu_bin = phat_mu_w if `inwin10'==1 & !missing(phat_mu_w), nq($PHAT_BINS_N)
    xtile phat_sd_bin = phat_sd_w if `inwin10'==1 & !missing(phat_sd_w), nq($PHAT_BINS_N)

}
else {

    xtile phat_mu_bin = phat_mu_wc if `inwin10'==1 & !missing(phat_mu_wc), nq($PHAT_BINS_N)
    xtile phat_sd_bin = phat_sd_wc if `inwin10'==1 & !missing(phat_sd_wc), nq($PHAT_BINS_N)

}

* --- DiD with bin interactions (omit bin 1 as baseline)
reghdfe $Y_HET ///
    i.treated_BC_3##i.post_3 ///
    i.phat_mu_bin##i.treated_BC_3##i.post_3 ///
    i.phat_sd_bin##i.treated_BC_3##i.post_3 ///
    $X_3 $LIFECYCLE_3 ///
    if `inwin10'==1 ///
       & !missing(phat_mu_bin) ///
       & !missing(phat_sd_bin), ///
    absorb(worker_id contract_id `ABS_TIME_3') ///
    $VCE_3

estimates store het_phat_bins

/**********************************************************************
  3.11) DIFFICULTY HETEROGENEITY: contract difficulty + weather shock
**********************************************************************/
reghdfe $Y_HET ///
    i.treated_BC_3##i.post_3 ///
    c.contract_diff_z##i.treated_BC_3##i.post_3 ///
    c.wx_shock##i.treated_BC_3##i.post_3 ///
    $X_3 $LIFECYCLE_3 ///
    if sample_3==1 ///
       & inrange(event_week_3,-$WPRE_3,$WPOST_3) ///
       & ($SAMPLE_OK_3) ///
       & !missing($Y_HET) ///
       & !missing(contract_diff_z) ///
       & !missing(wx_shock), ///
    absorb(worker_id contract_id `ABS_TIME_3') ///
    $VCE_3
estimates store het_difficulty


/**********************************************************************
  3.12) ABILITY HETEROGENEITY (bins) — kept, lower emphasis
**********************************************************************/
cap drop abil_bin_all
xtile abil_bin_all = ability_z if sample_3==1 & !missing(ability_z), nq(4)

reghdfe $Y_HET ///
    i.treated_BC_3##i.post_3 ///
    i.abil_bin_all##i.treated_BC_3##i.post_3 ///
    $X_3 $LIFECYCLE_3 ///
    if sample_3==1 ///
       & inrange(event_week_3,-$WPRE_3,$WPOST_3) ///
       & ($SAMPLE_OK_3) ///
       & !missing($Y_HET) ///
       & !missing(ability_z), ///
    absorb(worker_id contract_id `ABS_TIME_3') ///
    $VCE_3
estimates store het_ability_bins


/**********************************************************************
  3.13) DESCRIPTIVE TABLE: ability bins × pre/post (BC only)
**********************************************************************/
table abil_bin_all post_3 if (year==$TY_3) & (province_s=="BC") ///
    & inrange(event_week_3,-$WPRE_3,$WPOST_3), ///
    statistic(mean piece) statistic(mean prod) statistic(mean tot_hrs_day) nformat(%9.3f)

/**********************************************************************
  3.14) QUICK PRINT: key heterogeneity estimates (copy/paste friendly)
  - Robust to: omitted coeffs, different base coding (1bn/0b), model variants
**********************************************************************/

capture program drop _trylincom
program define _trylincom
    syntax , COEF(string) [LABEL(string)]
    capture noisily lincom `coef'
    if _rc==0 {
        if "`label'"!="" di as text "   -> `label'"
    }
    else {
        if "`label'"!="" di as text "   -> `label' (not in model / omitted)"
        else di as text "   -> `coef' (not in model / omitted)"
    }
end

capture program drop _findcoef
program define _findcoef, rclass
    // Find first coefficient whose name contains all required substrings
    // and none of the excluded substrings.
    syntax , REQUIRE1(string) REQUIRE2(string) [REQUIRE3(string)] ///
             [EXCL1(string) EXCL2(string) EXCL3(string) EXCL4(string)]
    tempname b
    matrix `b' = e(b)
    local found ""
    local cnames : colnames `b'
    foreach cn of local cnames {
        if strpos("`cn'","`require1'")==0 continue
        if strpos("`cn'","`require2'")==0 continue
        if "`require3'"!="" & strpos("`cn'","`require3'")==0 continue

        if "`excl1'"!="" & strpos("`cn'","`excl1'")>0 continue
        if "`excl2'"!="" & strpos("`cn'","`excl2'")>0 continue
        if "`excl3'"!="" & strpos("`cn'","`excl3'")>0 continue
        if "`excl4'"!="" & strpos("`cn'","`excl4'")>0 continue

        local found "`cn'"
        continue, break
    }
    return local coef "`found'"
end

di as text "============================================================"
di as text "SECTION 3: KEY RESULTS (copy/paste)"
di as text "============================================================"

* -----------------------------
* (A) Parametric phat model
* -----------------------------
capture quietly estimates restore het_phat_param
if _rc {
    di as error "het_phat_param not found"
}
else {
    di as text "---- Parametric phat heterogeneity: het_phat_param ----"

    * baseline DiD term name (robust to 1bn/0b coding)
    quietly _findcoef, require1("treated_BC_3") require2("post_3") ///
        excl1("c.phat") excl2("phat_") excl3("abil_bin") excl4("contract_diff")
    local did = r(coef)
    _trylincom, coef("`did'") label("Baseline DiD")

    * phat mean slope: treated#post#c.phat_mu_*
    quietly _findcoef, require1("treated_BC_3") require2("post_3") require3("c.phat_mu_") ///
        excl1("") excl2("") excl3("") excl4("")
    local mu = r(coef)
    if "`mu'"!="" _trylincom, coef("`mu'") label("DiD heterogeneity wrt mean phat")

    * phat sd slope: treated#post#c.phat_sd_*
    quietly _findcoef, require1("treated_BC_3") require2("post_3") require3("c.phat_sd_") ///
        excl1("") excl2("") excl3("") excl4("")
    local sd = r(coef)
    if "`sd'"!="" _trylincom, coef("`sd'") label("DiD heterogeneity wrt sd phat")

    di as text "   (sanity) moments:"
    capture noisily summ phat_mu_wc phat_sd_wc if sample_3==1
    capture noisily summ phat_mu_w  phat_sd_w  if sample_3==1
}

di as text ""

* -----------------------------
* (B) phat bins model (nonparametric)
* -----------------------------
capture quietly estimates restore het_phat_bins
if _rc {
    di as error "het_phat_bins not found"
}
else {
    di as text "---- Nonparametric phat bins: het_phat_bins ----"

    quietly _findcoef, require1("treated_BC_3") require2("post_3") ///
        excl1("phat_mu_bin") excl2("phat_sd_bin") excl3("c.phat") excl4("abil_bin")
    local didb = r(coef)
    _trylincom, coef("`didb'") label("Baseline DiD (omitted phat bin)")

    * Print the bin increments that exist (2-4 are present in your run; 5 isn't)
    foreach k in 2 3 4 5 {
        _trylincom, coef("`k'.phat_mu_bin#1.treated_BC_3#1.post_3") ///
            label("Mu-bin `k' increment")
    }
    foreach k in 2 3 4 5 {
        _trylincom, coef("`k'.phat_sd_bin#1.treated_BC_3#1.post_3") ///
            label("SD-bin `k' increment")
    }
}

di as text ""

* -----------------------------
* (C) Difficulty model
* -----------------------------
capture quietly estimates restore het_difficulty
if _rc {
    di as error "het_difficulty not found"
}
else {
    di as text "---- Difficulty heterogeneity: het_difficulty ----"

    quietly _findcoef, require1("treated_BC_3") require2("post_3") ///
        excl1("contract_diff_z") excl2("wx_shock") excl3("c.") excl4("abil_bin")
    local didd = r(coef)
    _trylincom, coef("`didd'") label("Baseline DiD")

    quietly _findcoef, require1("treated_BC_3") require2("post_3") require3("contract_diff_z")
    local cd = r(coef)
    if "`cd'"!="" _trylincom, coef("`cd'") label("DiD heterogeneity wrt contract difficulty")

    quietly _findcoef, require1("treated_BC_3") require2("post_3") require3("wx_shock")
    local wx = r(coef)
    if "`wx'"!="" _trylincom, coef("`wx'") label("DiD heterogeneity wrt weather shock")
}

di as text ""

* -----------------------------
* (D) Ability bins model (lower priority)
* -----------------------------
capture quietly estimates restore het_ability_bins
if _rc {
    di as error "het_ability_bins not found"
}
else {
    di as text "---- Ability-bin heterogeneity: het_ability_bins ----"

    quietly _findcoef, require1("treated_BC_3") require2("post_3") ///
        excl1("abil_bin_all") excl2("c.") excl3("") excl4("")
    local dida = r(coef)
    _trylincom, coef("`dida'") label("Baseline DiD (omitted ability bin)")

    foreach k in 2 3 4 {
        _trylincom, coef("`k'.abil_bin_all#1.treated_BC_3#1.post_3") label("Ability-bin `k' increment")
    }
}

di as text "============================================================"





/*
/**********************************************************************
  3.14) QUICK INTERPRETATION BLOCK (ROBUST TO SPEC CHANGES)
**********************************************************************/

di as text "============================================================"
di as text "SECTION 3 – KEY RESULTS"
di as text "============================================================"

*******************************************************
* (1) PARAMETRIC PHAT MODEL
*******************************************************
capture noisily estimates restore het_phat_param
if _rc==0 {

    di as text "---- Parametric phat heterogeneity ----"

    * Find baseline DiD
    local tp ""
    foreach cn of local colnames : colnames e(b) {
        if strpos("`cn'","treated_BC_3#") & ///
           strpos("`cn'","#post_3") & ///
           !strpos("`cn'","#c.") {
            local tp "`cn'"
            continue, break
        }
    }

    if "`tp'"!="" {
        di as text "Baseline DiD:"
        lincom `tp'
    }

    * Find mu triple interaction automatically
    foreach cn of local colnames : colnames e(b) {
        if strpos("`cn'","#c.phat_mu") & ///
           strpos("`cn'","treated_BC_3") & ///
           strpos("`cn'","post_3") {
            di as text "Heterogeneity by mean binding risk:"
            lincom `cn'
        }
    }

    * Find sd triple interaction automatically
    foreach cn of local colnames : colnames e(b) {
        if strpos("`cn'","#c.phat_sd") & ///
           strpos("`cn'","treated_BC_3") & ///
           strpos("`cn'","post_3") {
            di as text "Heterogeneity by volatility of binding risk:"
            lincom `cn'
        }
    }

    di as text "N = " %9.0f e(N) "   Within R2 = " %6.3f e(r2_within)
}

*******************************************************
* (2) BINNED PHAT MODEL
*******************************************************
capture noisily estimates restore het_phat_bins
if _rc==0 {

    di as text "---- Nonparametric phat heterogeneity (bins) ----"

    foreach cn of local colnames : colnames e(b) {
        if strpos("`cn'","phat_") & ///
           strpos("`cn'","treated_BC_3") & ///
           strpos("`cn'","post_3") {
            di as text "`cn'"
            lincom `cn'
        }
    }

    di as text "N = " %9.0f e(N) "   Within R2 = " %6.3f e(r2_within)
}

*******************************************************
* (3) DIFFICULTY MODEL
*******************************************************
capture noisily estimates restore het_difficulty
if _rc==0 {

    di as text "---- Contract + Weather Difficulty ----"

    foreach cn of local colnames : colnames e(b) {
        if strpos("`cn'","contract_diff_z") & ///
           strpos("`cn'","treated_BC_3") & ///
           strpos("`cn'","post_3") {
            di as text "Contract difficulty heterogeneity:"
            lincom `cn'
        }
        if strpos("`cn'","wx_shock") & ///
           strpos("`cn'","treated_BC_3") & ///
           strpos("`cn'","post_3") {
            di as text "Weather difficulty heterogeneity:"
            lincom `cn'
        }
    }

    di as text "N = " %9.0f e(N) "   Within R2 = " %6.3f e(r2_within)
}

di as text "============================================================"




/*
/**********************************************************************
  (3) ADDITIONAL EMPIRICAL ANALYSIS (SELF-CONTAINED)
  Focus: BC vs AB in treated year (2A-style), heterogeneity by:
    (i) predicted top-up risk (MW-binding likelihood) + its variance
    (ii) job difficulty: contract-level + daily weather difficulty

  KEY DESIGN CHOICES (per your notes):
    - Ability: estimated PRE only (post_3==0), then propagated to all rows by worker.
    - Contract difficulty: estimated using AB-only data; toggle whether AB post is allowed.
      (AB is untreated, so using AB post can improve support if you believe no spillovers.)
    - phat_topup: time-varying by construction (mw + weather), estimated using PRE only.
    - Variance of phat_topup: compute within analysis window at worker or wc_id level.
**********************************************************************/

***********************************************************************
* 3.0) CLEAN RERUN STATE (SAFE)
***********************************************************************
cap estimates clear

foreach v in ///
    sample_3 treated_BC_3 post_3 event_week_3 event_week_3s ///
    ln_prod abil_sample resid_preprod ability_raw ability_z ///
    topup_ind topup_pre phat_topup ///
    cd_sample resid_cd contract_diff_raw contract_diff_z ///
    wx_shock ///
    phat_mu_w phat_sd_w phat_mu_wc phat_sd_wc ///
    abil_bin_all {
    cap drop `v'
}

foreach m in ///
    TY_3 CONTROL_3 Y_HET CLUSTER_3 VCE_3 ///
    USE_FLOORS_3 MIN_HRS_3 MIN_EARN_3 MIN_PROD_3 USE_MIN_PROD_3 ///
    USE_WEATHER_3 USE_EXPER_3 USE_HOURS_3 USE_DAYFLAGS_3 USE_LIFECYCLE_3 ///
    WPRE_3 WPOST_3 ///
    CD_USE_AB_POST CD_TIMEFE ///
    PHAT_VAR_LEVEL ///
    X_3 LIFECYCLE_3 SAMPLE_OK_3 {
    cap macro drop `m'
}

***********************************************************************
* 3.1) SETTINGS (EDIT THESE)
***********************************************************************
* Treated year and control province
global TY_3      2018
global CONTROL_3 "AB"

* Outcome for heterogeneity regressions
global Y_HET "ln_incentive"

* Clustering toggle: worker_id / contract_id / wc_id
global CLUSTER_3 "worker_id"
global VCE_3 "vce(cluster $CLUSTER_3)"

* Controls toggles (mirror Section 2 logic)
global USE_WEATHER_3   1
global USE_EXPER_3     1
global USE_HOURS_3     0
global USE_DAYFLAGS_3  0
global USE_LIFECYCLE_3 1

* Floors toggle (0 = none)
global USE_FLOORS_3   0
global MIN_HRS_3      2
global MIN_EARN_3     30
global USE_MIN_PROD_3 1
global MIN_PROD_3     50

* Event window (manual: keeps this section stable)
global WPRE_3  3
global WPOST_3 8

* ---------------------------
* NEW TOGGLES
* ---------------------------

* (A) Contract difficulty estimation sample (AB-only):
*   0 = AB pre-only (conservative)
*   1 = AB pre+post (uses more AB support; assumes no spillovers into AB)
global CD_USE_AB_POST 1

* (A2) Add common time effects in AB-only difficulty estimation:
*   0 = none beyond controls
*   1 = include event-week FE (recommended)
global CD_TIMEFE 1

* (B) Binding-risk "variance" level:
*   "worker" = mean/sd of phat_topup by worker_id
*   "wc"     = mean/sd by wc_id (worker×contract)
global PHAT_VAR_LEVEL "wc"

***********************************************************************
* 3.2) BUILD CONTROLS + FLOORS STRINGS
***********************************************************************
global X_3 ""
if $USE_HOURS_3==1    global X_3 "$X_3 tot_hrs_day"
if $USE_DAYFLAGS_3==1 global X_3 "$X_3 multi_contract_day"
if $USE_EXPER_3==1    global X_3 "$X_3 c.experience##c.experience"
if $USE_WEATHER_3==1  global X_3 "$X_3 c.temp##c.temp c.wind_spd##c.wind_spd precip0 c.precip_pos##c.precip_pos"

global LIFECYCLE_3 ""
if $USE_LIFECYCLE_3==1 {
    global LIFECYCLE_3 "$LIFECYCLE_3 c.contract_day##c.contract_day"
    global LIFECYCLE_3 "$LIFECYCLE_3 c.season_days_worked##c.season_days_worked"
    global LIFECYCLE_3 "$LIFECYCLE_3 c.season_days_worked##c.cum_seasons"
}

global SAMPLE_OK_3 "1==1"
if $USE_FLOORS_3==1 {
    global SAMPLE_OK_3 "$SAMPLE_OK_3 & tot_hrs_day >= $MIN_HRS_3"
    global SAMPLE_OK_3 "$SAMPLE_OK_3 & piece > 0 & prod > 0"
    global SAMPLE_OK_3 "$SAMPLE_OK_3 & (piece*prod) >= $MIN_EARN_3"
    if $USE_MIN_PROD_3==1 global SAMPLE_OK_3 "$SAMPLE_OK_3 & prod >= $MIN_PROD_3"
}

di as text "[3] TY_3:           $TY_3"
di as text "[3] CONTROL_3:      $CONTROL_3"
di as text "[3] Y_HET:          $Y_HET"
di as text "[3] CLUSTER:        $CLUSTER_3"
di as text "[3] Window:         [-$WPRE_3,$WPOST_3]"
di as text "[3] CD_USE_AB_POST: $CD_USE_AB_POST"
di as text "[3] CD_TIMEFE:      $CD_TIMEFE"
di as text "[3] PHAT_VAR_LEVEL: $PHAT_VAR_LEVEL"
di as text "[3] X_3:            $X_3"
di as text "[3] LIFECYCLE_3:    $LIFECYCLE_3"
di as text "[3] FLOORS:         $SAMPLE_OK_3"

***********************************************************************
* 3.3) DEFINE 2A-LIKE SAMPLE (BC vs AB in treated year) + EVENT TIME
***********************************************************************
gen byte sample_3 = (year==$TY_3) & inlist(province_s,"BC","$CONTROL_3")

gen byte treated_BC_3 = .
replace treated_BC_3 = (province_s=="BC") if sample_3==1

local tdate_3 = mdy(6,1,$TY_3)
gen byte post_3       = (date >= `tdate_3') if sample_3==1
gen int  event_week_3 = floor((date - `tdate_3')/7) if sample_3==1

* IMPORTANT: Stata factor vars cannot be negative; create shifted week for FE use.
gen int event_week_3s = event_week_3 + 100 if sample_3==1

tab province_s if sample_3==1
tab post_3     if sample_3==1

/**********************************************************************
  3.4) PRE-TREATMENT ABILITY (BC+AB, PRE ONLY)
  - ln(prod) residualized on controls + lifecycle + contract FE, PRE only
  - ability(worker) = mean residual, propagated to all rows of worker
**********************************************************************/
cap drop ln_prod
gen double ln_prod = ln(prod) if prod>0

cap drop abil_sample
gen byte abil_sample = (sample_3==1) ///
    & (post_3==0) ///
    & !missing(ln_prod) ///
    & ($SAMPLE_OK_3)

di as text "[3.4] Ability PRE sample:"
count if abil_sample==1
tab province_s if abil_sample==1

cap drop resid_preprod
reghdfe ln_prod $X_3 $LIFECYCLE_3 if abil_sample==1, ///
    absorb(contract_id) ///
    resid(resid_preprod)

cap drop ability_raw ability_z
bysort worker_id: egen double ability_raw = mean(resid_preprod)

quietly summarize ability_raw if !missing(ability_raw)
scalar abil_mu = r(mean)
scalar abil_sd = r(sd)

gen double ability_z = (ability_raw - abil_mu)/abil_sd ///
    if !missing(ability_raw) & abil_sd>0

di as text "[3.4] Ability coverage:"
count if sample_3==1 & !missing(ability_z)
tab province_s if sample_3==1 & !missing(ability_z)
tab post_3 if sample_3==1, summarize(ability_z)
bys worker_id: assert ability_z==ability_z[1] if !missing(ability_z)

**********************************************************************
* 3.5) DAILY WEATHER SHOCK (not MW-contaminated)
*   wx_shock = fitted weather component from a simple regression on weather only
*   (you can keep this intentionally "reduced form" / stable)
**********************************************************************
cap drop wx_shock
quietly reg ln_prod ///
    c.temp##c.temp c.wind_spd##c.wind_spd precip0 c.precip_pos##c.precip_pos ///
    if sample_3==1 & !missing(ln_prod), vce(robust)
predict double wx_shock if sample_3==1, xb
summ wx_shock if sample_3==1

/**********************************************************************
  3.6) CONTRACT DIFFICULTY (AB-only; toggle pre-only vs pre+post; optional week FE)
  Goal: get contract_diff_z with good AB support (so it's not missing),
        then carry it into the BC vs AB heterogeneity regression.

  Approach:
    - Use AB-only data (year==TY_3, province==AB), optionally incl post
    - Residualize ln_prod on ability + controls + lifecycle
    - Absorb worker_id (and optional event-week FE via event_week_3s)
    - Contract difficulty = mean residual by contract_id
**********************************************************************/
cap drop cd_sample resid_cd contract_diff_raw contract_diff_z

gen byte cd_sample = (year==$TY_3) & (province_s=="AB") ///
    & !missing(ln_prod) ///
    & !missing(ability_z) ///
    & ($SAMPLE_OK_3)

if $CD_USE_AB_POST==0 {
    replace cd_sample = cd_sample & (post_3==0)
}

di as text "[3.6] Contract-difficulty sample (AB-only):"
count if cd_sample==1
tab post_3 if cd_sample==1

* AB-only residualization (include common time FE via event_week_3s if toggled on)
if $CD_TIMEFE==1 {
    reghdfe ln_prod ///
        c.ability_z ///
        $X_3 ///
        $LIFECYCLE_3 ///
        if cd_sample==1, ///
        absorb(worker_id event_week_3s) ///
        resid(resid_cd)
}
else {
    reghdfe ln_prod ///
        c.ability_z ///
        $X_3 ///
        $LIFECYCLE_3 ///
        if cd_sample==1, ///
        absorb(worker_id) ///
        resid(resid_cd)
}

* Contract scalar from AB sample; egen propagates to all rows with that contract_id.
bysort contract_id: egen double contract_diff_raw = mean(resid_cd)

quietly summarize contract_diff_raw if !missing(contract_diff_raw)
scalar cd_mu = r(mean)
scalar cd_sd = r(sd)

gen double contract_diff_z = (contract_diff_raw - cd_mu)/cd_sd ///
    if !missing(contract_diff_raw) & cd_sd>0

di as text "[3.6] Contract difficulty coverage (by province in sample_3):"
count if sample_3==1 & province_s=="AB" & !missing(contract_diff_z)
count if sample_3==1 & province_s=="BC" & !missing(contract_diff_z)

**********************************************************************
* 3.7) PREDICTED TOP-UP RISK phat_topup (time-varying; uses mw + wx_shock)
*   - Estimate PRE only (avoid training on treatment)
*   - Predict for ALL sample_3 observations (so it updates with mw + weather)
**********************************************************************
cap drop topup_ind topup_pre phat_topup
gen byte topup_ind = (top_up > 0) if !missing(top_up)

gen byte topup_pre = (sample_3==1) ///
    & (post_3==0) ///
    & !missing(topup_ind) ///
    & !missing(ability_z) ///
    & !missing(mw) ///
    & !missing(wx_shock) ///
    & ($SAMPLE_OK_3)

di as text "[3.7] Top-up PRE estimation sample:"
count if topup_pre==1
tab province_s if topup_pre==1
summ topup_ind if topup_pre==1

* LPM is stable; include contract FE to capture baseline contract binding propensity
reg topup_ind ///
    c.mw ///
    c.ability_z ///
    c.wx_shock ///
    $X_3 ///
    i.contract_id ///
    if topup_pre==1, ///
    vce(cluster $CLUSTER_3)

predict double phat_topup if sample_3==1, xb
replace phat_topup = max(0, min(1, phat_topup)) if sample_3==1 & !missing(phat_topup)

di as text "[3.7] phat_topup coverage:"
summ phat_topup if sample_3==1
tab post_3 if sample_3==1, summarize(phat_topup)

**********************************************************************
* 3.8) MEAN + VOLATILITY OF BINDING RISK (within analysis window)
*   - compute on the same window you estimate heterogeneity on
*   - level toggle: worker vs wc_id
**********************************************************************
cap drop phat_mu_w phat_sd_w phat_mu_wc phat_sd_wc

tempvar inwin
gen byte `inwin' = (sample_3==1) ///
    & inrange(event_week_3,-$WPRE_3,$WPOST_3) ///
    & ($SAMPLE_OK_3) ///
    & !missing($Y_HET) ///
    & !missing(phat_topup)

if "$PHAT_VAR_LEVEL"=="worker" {
    bysort worker_id: egen double phat_mu_w = mean(cond(`inwin'==1, phat_topup, .))
    bysort worker_id: egen double phat_sd_w = sd(  cond(`inwin'==1, phat_topup, .))
}
else if "$PHAT_VAR_LEVEL"=="wc" {
    capture confirm variable wc_id
    if _rc {
        di as error "[3.8] wc_id not found; using worker-level moments instead."
        global PHAT_VAR_LEVEL "worker"
        bysort worker_id: egen double phat_mu_w = mean(cond(`inwin'==1, phat_topup, .))
        bysort worker_id: egen double phat_sd_w = sd(  cond(`inwin'==1, phat_topup, .))
    }
    else {
        bysort wc_id: egen double phat_mu_wc = mean(cond(`inwin'==1, phat_topup, .))
        bysort wc_id: egen double phat_sd_wc = sd(  cond(`inwin'==1, phat_topup, .))
    }
}

di as text "[3.8] phat moment diagnostics:"
count if `inwin'==1
summ phat_topup if `inwin'==1
if "$PHAT_VAR_LEVEL"=="worker" {
    summ phat_mu_w phat_sd_w if sample_3==1
}
else {
    summ phat_mu_wc phat_sd_wc if sample_3==1
}

**********************************************************************
* 3.9) MAIN HETEROGENEITY DiD: MW effect vs binding-risk mean + volatility
*   Key coefficients (depending on PHAT_VAR_LEVEL):
*     - treated#post#phat_mu_* : heterogeneity in DiD by average binding risk
*     - treated#post#phat_sd_* : heterogeneity in DiD by volatility of binding risk
**********************************************************************
if "$PHAT_VAR_LEVEL"=="worker" {

    reghdfe $Y_HET ///
        i.treated_BC_3##i.post_3 ///
        c.phat_mu_w##i.treated_BC_3##i.post_3 ///
        c.phat_sd_w##i.treated_BC_3##i.post_3 ///
        $X_3 ///
        $LIFECYCLE_3 ///
        if sample_3==1 ///
           & inrange(event_week_3,-$WPRE_3,$WPOST_3) ///
           & ($SAMPLE_OK_3) ///
           & !missing($Y_HET) ///
           & !missing(phat_mu_w) ///
           & !missing(phat_sd_w), ///
        absorb(worker_id contract_id) ///
        $VCE_3
    estimates store het_topup_mu_sd

}
else {

    reghdfe $Y_HET ///
        i.treated_BC_3##i.post_3 ///
        c.phat_mu_wc##i.treated_BC_3##i.post_3 ///
        c.phat_sd_wc##i.treated_BC_3##i.post_3 ///
        $X_3 ///
        $LIFECYCLE_3 ///
        if sample_3==1 ///
           & inrange(event_week_3,-$WPRE_3,$WPOST_3) ///
           & ($SAMPLE_OK_3) ///
           & !missing($Y_HET) ///
           & !missing(phat_mu_wc) ///
           & !missing(phat_sd_wc), ///
        absorb(worker_id contract_id) ///
        $VCE_3
    estimates store het_topup_mu_sd
}

**********************************************************************
* 3.10) DIFFICULTY HETEROGENEITY: contract difficulty + weather shock
*   NOTE:
*     - contract_diff_z level is contract-invariant (collinear with contract FE),
*       but treated#post#contract_diff_z is the heterogeneity object.
*     - If it is omitted: it's almost always because AB support collapses (missingness),
*       or because only 1 AB contract remains in e(sample).
**********************************************************************/
reghdfe $Y_HET ///
    i.treated_BC_3##i.post_3 ///
    c.contract_diff_z##i.treated_BC_3##i.post_3 ///
    c.wx_shock##i.treated_BC_3##i.post_3 ///
    $X_3 ///
    $LIFECYCLE_3 ///
    if sample_3==1 ///
       & inrange(event_week_3,-$WPRE_3,$WPOST_3) ///
       & ($SAMPLE_OK_3) ///
       & !missing($Y_HET) ///
       & !missing(contract_diff_z) ///
       & !missing(wx_shock), ///
    absorb(worker_id contract_id) ///
    $VCE_3
estimates store het_difficulty

**********************************************************************
* 3.11) NONLINEAR ABILITY HETEROGENEITY (BINS) — quick diagnostic
**********************************************************************
cap drop abil_bin_all
xtile abil_bin_all = ability_z if sample_3==1 & !missing(ability_z), nq(4)

reghdfe $Y_HET ///
    i.treated_BC_3##i.post_3 ///
    i.abil_bin_all##i.treated_BC_3##i.post_3 ///
    $X_3 ///
    $LIFECYCLE_3 ///
    if sample_3==1 ///
       & inrange(event_week_3,-$WPRE_3,$WPOST_3) ///
       & ($SAMPLE_OK_3) ///
       & !missing($Y_HET) ///
       & !missing(ability_z), ///
    absorb(worker_id contract_id) ///
    $VCE_3
estimates store het_bins

**********************************************************************
* 3.12) DESCRIPTIVE TABLE: ability bins × pre/post (BC only)
**********************************************************************
table abil_bin_all post_3 if (year==$TY_3) & (province_s=="BC") ///
    & inrange(event_week_3,-$WPRE_3,$WPOST_3), ///
    statistic(mean piece) statistic(mean prod) statistic(mean tot_hrs_day) nformat(%9.3f)
	
	

	
	
	
	/*
/**********************************************************************
  (3) ADDITIONAL EMPIRICAL ANALYSIS (SELF-CONTAINED)

  Focus: BC vs AB in treated year (2A-style), heterogeneity by:
    (i) MW-binding risk (predicted top-up probability) using PRE-only
        productivity components + contract-mean piece + weather variation
    (ii) job difficulty: contract-level + daily weather difficulty

  NOTE:
    - This section is self-contained (does NOT rely on Section 2 globals).
    - It is written to be safe to rerun (cap drop everywhere).
**********************************************************************/

***********************************************************************
* 3.0) CLEAN RERUN STATE (SAFE)
***********************************************************************
cap estimates clear

foreach v in ///
    sample_3 treated_BC_3 post_3 event_week_3 ///
    ln_prod abil_sample resid_preprod ability_raw ability_z ///
    resid_cd contract_diff_raw contract_diff_z ///
    piece_pre ///
    mw_obs mw_rate mw_rate_fallback ///
    earn_hat gap_hat ///
    topup_ind topup_pre phat_topup ///
    resid_wx wx_shock ///
    abil_bin_all {
    cap drop `v'
}

foreach m in ///
    TY_3 CONTROL_3 Y_HET CLUSTER_3 VCE_3 ///
    USE_FLOORS_3 MIN_HRS_3 MIN_EARN_3 MIN_PROD_3 USE_MIN_PROD_3 ///
    USE_WEATHER_3 USE_EXPER_3 USE_HOURS_3 USE_DAYFLAGS_3 USE_LIFECYCLE_3 ///
    WPRE_3 WPOST_3 ///
    X_3 LIFECYCLE_3 SAMPLE_OK_3 ///
    X_WX_3 {
    cap macro drop `m'
}

***********************************************************************
* 3.1) SETTINGS (EDIT THESE)
***********************************************************************
* Treated year and control province
global TY_3      2018
global CONTROL_3 "AB"

* Outcome for heterogeneity regressions
global Y_HET "ln_incentive"

* Clustering toggle: worker_id / contract_id / wc_id
global CLUSTER_3 "worker_id"

* Controls toggles (mirror Section 2 logic)
global USE_WEATHER_3   1
global USE_EXPER_3     1
global USE_HOURS_3     0
global USE_DAYFLAGS_3  0
global USE_LIFECYCLE_3 1

* Floors toggle (0 = none)
global USE_FLOORS_3   0
global MIN_HRS_3      2
global MIN_EARN_3     30
global USE_MIN_PROD_3 1
global MIN_PROD_3     50

* Event window (manual: keeps this section stable)
global WPRE_3  3
global WPOST_3 8

***********************************************************************
* 3.2) BUILD CONTROLS + FLOORS STRINGS + VCE STRING
***********************************************************************
global X_3 ""
if $USE_HOURS_3==1    global X_3 "$X_3 tot_hrs_day"
if $USE_DAYFLAGS_3==1 global X_3 "$X_3 multi_contract_day"
if $USE_EXPER_3==1    global X_3 "$X_3 c.experience##c.experience"
if $USE_WEATHER_3==1  global X_3 "$X_3 c.temp##c.temp c.wind_spd##c.wind_spd precip0 c.precip_pos##c.precip_pos"

* Weather-only block (handy for "wx shock" prediction later)
global X_WX_3 ""
if $USE_WEATHER_3==1  global X_WX_3 "$X_WX_3 c.temp##c.temp c.wind_spd##c.wind_spd precip0 c.precip_pos##c.precip_pos"

global LIFECYCLE_3 ""
if $USE_LIFECYCLE_3==1 {
    global LIFECYCLE_3 "$LIFECYCLE_3 c.contract_day##c.contract_day"
    global LIFECYCLE_3 "$LIFECYCLE_3 c.season_days_worked##c.season_days_worked"
    global LIFECYCLE_3 "$LIFECYCLE_3 c.season_days_worked##c.cum_seasons"
}

global SAMPLE_OK_3 "1==1"
if $USE_FLOORS_3==1 {
    global SAMPLE_OK_3 "$SAMPLE_OK_3 & tot_hrs_day >= $MIN_HRS_3"
    global SAMPLE_OK_3 "$SAMPLE_OK_3 & piece > 0 & prod > 0"
    global SAMPLE_OK_3 "$SAMPLE_OK_3 & (piece*prod) >= $MIN_EARN_3"
    if $USE_MIN_PROD_3==1 global SAMPLE_OK_3 "$SAMPLE_OK_3 & prod >= $MIN_PROD_3"
}

* Single place to toggle VCE
global VCE_3 "vce(cluster $CLUSTER_3)"

di as text "[3] TY_3:         $TY_3"
di as text "[3] CONTROL_3:    $CONTROL_3"
di as text "[3] Y_HET:        $Y_HET"
di as text "[3] CLUSTER:      $CLUSTER_3"
di as text "[3] X_3:          $X_3"
di as text "[3] LIFECYCLE_3:  $LIFECYCLE_3"
di as text "[3] FLOORS:       $SAMPLE_OK_3"
di as text "[3] Window:       [-$WPRE_3,$WPOST_3]"

***********************************************************************
* 3.3) DEFINE 2A-LIKE SAMPLE (BC vs AB in treated year) + EVENT TIME
***********************************************************************
cap drop sample_3 treated_BC_3 post_3 event_week_3
gen byte sample_3 = (year==$TY_3) & inlist(province_s,"BC","$CONTROL_3")

gen byte treated_BC_3 = .
replace treated_BC_3 = (province_s=="BC") if sample_3==1

local tdate_3 = mdy(6,1,$TY_3)
gen byte post_3       = (date >= `tdate_3') if sample_3==1
gen int  event_week_3 = floor((date - `tdate_3')/7) if sample_3==1

tab province_s if sample_3==1
tab post_3     if sample_3==1
summ event_week_3 if sample_3==1

/**********************************************************************
  3.4) PRE-TREATMENT ABILITY (BC+AB, PRE ONLY)
  - ln(prod) residualized on controls + lifecycle + contract FE
  - estimated in PRE only (post_3==0)
  - ability(worker) = mean residual (propagates to all rows of worker)
**********************************************************************/
cap drop ln_prod abil_sample resid_preprod ability_raw ability_z
gen double ln_prod = ln(prod) if prod>0

gen byte abil_sample = (sample_3==1) ///
    & (post_3==0) ///
    & !missing(ln_prod) ///
    & ($SAMPLE_OK_3)

di as text "[3.4] Ability PRE sample:"
count if abil_sample==1
tab province_s if abil_sample==1

reghdfe ln_prod $X_3 $LIFECYCLE_3 if abil_sample==1, ///
    absorb(contract_id) ///
    resid(resid_preprod)

bysort worker_id: egen double ability_raw = mean(resid_preprod)

quietly summarize ability_raw if !missing(ability_raw)
scalar abil_mu = r(mean)
scalar abil_sd = r(sd)

gen double ability_z = (ability_raw - abil_mu)/abil_sd ///
    if !missing(ability_raw) & abil_sd>0

di as text "[3.4] Ability coverage:"
count if sample_3==1 & !missing(ability_z)
tab province_s if sample_3==1 & !missing(ability_z)
tab post_3 if sample_3==1, summarize(ability_z)
bys worker_id: assert ability_z==ability_z[1] if !missing(ability_z)

***********************************************************************
* 3.5) CONTRACT DIFFICULTY INDEX (PRE ONLY)
* Contract difficulty should be PRE-determined, so estimate in PRE only.
*
* Construction:
*   - residualize ln_prod on ability_z + controls + lifecycle
*     absorbing worker FE (so residual is contract/day component)
*   - contract_diff_raw = mean residual by contract in PRE
*   - standardize to contract_diff_z (within sample_3 support)
***********************************************************************
cap drop resid_cd contract_diff_raw contract_diff_z

reghdfe ln_prod ///
    c.ability_z ///
    $X_3 ///
    $LIFECYCLE_3 ///
    if abil_sample==1 & !missing(ability_z), ///
    absorb(worker_id) ///
    resid(resid_cd)

bysort contract_id: egen double contract_diff_raw = mean(resid_cd)

quietly summarize contract_diff_raw if sample_3==1 & !missing(contract_diff_raw)
scalar cd_mu = r(mean)
scalar cd_sd = r(sd)

gen double contract_diff_z = (contract_diff_raw - cd_mu)/cd_sd ///
    if sample_3==1 & !missing(contract_diff_raw) & cd_sd>0

di as text "[3.5] Contract difficulty coverage:"
count if sample_3==1 & !missing(contract_diff_z)
summ contract_diff_z if sample_3==1

* What provinces/contracts are actually in the estimation sample?
tab province_s if e(sample)
tab treated_BC_3 if e(sample)

* Is contract_diff_z missing systematically by province?
count if sample_3==1 & province_s=="AB" & missing(contract_diff_z)
count if sample_3==1 & province_s=="AB" & !missing(contract_diff_z)
count if sample_3==1 & province_s=="BC" & !missing(contract_diff_z)

* How many contracts by province in e(sample)?
bys contract_id: gen byte one = (_n==1) if e(sample)
tab province_s if one==1
drop one


/**********************************************************************
  3.6) CONTRACT-MEAN PIECE RATE (PRE ONLY)  ***YOUR CHOICE***
  - piece_pre(contract) = mean(piece) among PRE observations in sample_3
  - propagate to all rows of that contract
**********************************************************************/
cap drop piece_pre
bysort contract_id: egen double piece_pre = mean(piece) ///
    if sample_3==1 & post_3==0 & !missing(piece)

* Propagate contract mean piece to all rows of contract:
bysort contract_id: replace piece_pre = piece_pre[_n-1] if missing(piece_pre) & _n>1
bysort contract_id: replace piece_pre = piece_pre[_n+1] if missing(piece_pre) & _n<_N

di as text "[3.6] piece_pre coverage (sample_3):"
count if sample_3==1 & !missing(piece_pre)
summ piece_pre if sample_3==1 & !missing(piece_pre)

/**********************************************************************
  3.7) PREDICTED PRODUCTIVITY (TIME-VARYING VIA WEATHER, PRE-ESTIMATED COEFS)
  Goal: build prod_hat_{it} using only PRE-estimated relationships:
    ln_prod = b0 + b1*ability_z + b2*contract_diff_z + b3*weather + b4*experience + lifecycle + e
  - Estimate on PRE only (abil_sample==1)
  - Predict for ALL observations in sample_3 (pre and post)
**********************************************************************/
cap drop earn_hat gap_hat

* PRE-only fit (no FE here; contract_diff_z and ability_z carry the PRE structure)
reg ln_prod ///
    c.ability_z ///
    c.contract_diff_z ///
    $X_3 ///
    $LIFECYCLE_3 ///
    if abil_sample==1 ///
       & !missing(ability_z) ///
       & !missing(contract_diff_z), ///
    vce(cluster contract_id)

predict double ln_prod_hat if sample_3==1, xb
gen double prod_hat = exp(ln_prod_hat) if sample_3==1 & !missing(ln_prod_hat)

* Predicted hourly incentive earnings from PRE piece rate and predicted productivity
* (piece is $/tree, prod is trees/hour => $/hour)
gen double earn_hat = piece_pre * prod_hat if sample_3==1 & !missing(piece_pre) & !missing(prod_hat)

di as text "[3.7] earn_hat coverage:"
count if sample_3==1 & !missing(earn_hat)
summ earn_hat if sample_3==1 & !missing(earn_hat)

/**********************************************************************
  3.8) MINIMUM WAGE (OBSERVED) + GAP CONSTRUCTION
  - Use mw directly (statutory minimum wage, $/hour)
  - earn_hat = piece_pre * prod_hat ($/hour) from PRE-only objects
  - gap_hat  = mw - earn_hat (positive => predicted binding risk)
**********************************************************************/
cap drop gap_hat

* sanity: mw should exist and be non-missing in sample_3
count if sample_3==1 & missing(mw)
summ mw if sample_3==1

gen double gap_hat = (mw - earn_hat) if sample_3==1 & !missing(mw) & !missing(earn_hat)

di as text "[3.8] gap_hat coverage:"
count if sample_3==1 & !missing(gap_hat)
summ gap_hat if sample_3==1 & !missing(gap_hat)
tab post_3 if sample_3==1, summarize(gap_hat)


/**********************************************************************
  3.9) TOP-UP RISK: phat_topup via PRE-only mapping from gap_hat
  - topup_ind = 1{top_up>0}
  - Fit mapping using PRE only (post_3==0) so we don't train on treatment
  - Predict for ALL observations in sample_3
**********************************************************************/
cap drop topup_ind topup_pre phat_topup

gen byte topup_ind = (top_up > 0) if sample_3==1 & !missing(top_up)

gen byte topup_pre = (sample_3==1) ///
    & (post_3==0) ///
    & !missing(topup_ind) ///
    & !missing(gap_hat) ///
    & ($SAMPLE_OK_3)

di as text "[3.9] Top-up PRE estimation sample:"
count if topup_pre==1
tab province_s if topup_pre==1
summ topup_ind gap_hat if topup_pre==1

* LPM on gap_hat (robust + stable; avoids separation headaches)
reg topup_ind c.gap_hat if topup_pre==1, vce(cluster $CLUSTER_3)

predict double phat_topup if sample_3==1 & !missing(gap_hat), xb
replace phat_topup = max(0, min(1, phat_topup)) if sample_3==1 & !missing(phat_topup)

di as text "[3.9] phat_topup coverage:"
summ phat_topup if sample_3==1
tab post_3 if sample_3==1, summarize(phat_topup)


/**********************************************************************
  3.10) DAILY WEATHER "DIFFICULTY SHOCK" (wx_shock) — PRE-ESTIMATED
  Goal: isolate time-varying weather component only, estimated in PRE.
  Approach:
    - PRE-only regression of ln_prod on weather terms (and nothing else),
      then wx_shock = fitted weather component (xb) for all sample_3 days.
  Note: This is intentionally "pure weather" for time-variation.
**********************************************************************/
cap drop wx_shock
if "$X_WX_3"=="" {
    gen double wx_shock = .   // weather not turned on
}
else {
    reg ln_prod $X_WX_3 if abil_sample==1, vce(cluster contract_id)
    predict double wx_shock if sample_3==1, xb
}

di as text "[3.10] wx_shock coverage:"
summ wx_shock if sample_3==1

/**********************************************************************
  3.11) MAIN HETEROGENEITY DiD: does MW effect vary with binding risk?
  Key coefficient: c.phat_topup#1.treated_BC_3#1.post_3
**********************************************************************/
reghdfe $Y_HET ///
    i.treated_BC_3##i.post_3 ///
    c.phat_topup##i.treated_BC_3##i.post_3 ///
    $X_3 ///
    $LIFECYCLE_3 ///
    if sample_3==1 ///
       & inrange(event_week_3,-$WPRE_3,$WPOST_3) ///
       & !missing($Y_HET) ///
       & !missing(phat_topup) ///
       & ($SAMPLE_OK_3), ///
    absorb(worker_id contract_id) ///
    $VCE_3
estimates store het_topup

/**********************************************************************
  3.12) DIFFICULTY HETEROGENEITY: contract difficulty + daily weather
  - contract_diff_z level is contract-invariant (collinear with contract FE),
    but interactions with post/treat are identified.
**********************************************************************/
reghdfe $Y_HET ///
    i.treated_BC_3##i.post_3 ///
    c.contract_diff_z##i.treated_BC_3##i.post_3 ///
    c.wx_shock##i.treated_BC_3##i.post_3 ///
    $X_3 ///
    $LIFECYCLE_3 ///
    if sample_3==1 ///
       & inrange(event_week_3,-$WPRE_3,$WPOST_3) ///
       & !missing($Y_HET) ///
       & !missing(contract_diff_z) ///
       & !missing(wx_shock) ///
       & ($SAMPLE_OK_3), ///
    absorb(worker_id contract_id) ///
    $VCE_3
estimates store het_difficulty

/**********************************************************************
  3.13) NONLINEAR ABILITY HETEROGENEITY (BINS) — quick diagnostic
**********************************************************************/
cap drop abil_bin_all
xtile abil_bin_all = ability_z if sample_3==1 & !missing(ability_z), nq(4)

reghdfe $Y_HET ///
    i.treated_BC_3##i.post_3 ///
    i.abil_bin_all##i.treated_BC_3##i.post_3 ///
    $X_3 ///
    $LIFECYCLE_3 ///
    if sample_3==1 ///
       & inrange(event_week_3,-$WPRE_3,$WPOST_3) ///
       & !missing($Y_HET) ///
       & !missing(ability_z) ///
       & ($SAMPLE_OK_3), ///
    absorb(worker_id contract_id) ///
    $VCE_3
estimates store het_bins

/**********************************************************************
  3.14) DESCRIPTIVE TABLE: ability bins × pre/post (BC only)
**********************************************************************/
table abil_bin_all post_3 if (year==$TY_3) & (province_s=="BC") ///
    & inrange(event_week_3,-$WPRE_3,$WPOST_3), ///
    statistic(mean piece) statistic(mean prod) statistic(mean tot_hrs_day) nformat(%9.3f)
	
	
	
	
	
	
	
	
	
	/*
	
	
	
	
	
	

/**********************************************************************
  (3) ADDITIONAL EMPIRICAL ANALYSIS (SELF-CONTAINED)
  Focus: BC vs AB in treated year (2A-style), heterogeneity by:
    (i) predicted top-up risk (MW-binding likelihood)
    (ii) job difficulty: contract-level + daily weather difficulty
**********************************************************************/

***********************************************************************
* 3.0) CLEAN RERUN STATE (SAFE)
***********************************************************************
cap estimates clear

foreach v in ///
    sample_3 treated_BC_3 post_3 event_week_3 ///
    ln_prod abil_sample resid_preprod ability_raw ability_z ///
    topup_ind topup_pre phat_topup ///
    resid_wx wx_shock ///
    resid_cd contract_diff_raw contract_diff_z ///
    abil_bin_all {
    cap drop `v'
}

foreach m in ///
    TY_3 CONTROL_3 Y_HET CLUSTER_3 VCE_3 ///
    USE_FLOORS_3 MIN_HRS_3 MIN_EARN_3 MIN_PROD_3 USE_MIN_PROD_3 ///
    USE_WEATHER_3 USE_EXPER_3 USE_HOURS_3 USE_DAYFLAGS_3 USE_LIFECYCLE_3 ///
    WPRE_3 WPOST_3 ///
    X_3 LIFECYCLE_3 SAMPLE_OK_3 {
    cap macro drop `m'
}


***********************************************************************
* 3.1) SETTINGS (EDIT THESE)
***********************************************************************
* Treated year and control province
global TY_3      2018
global CONTROL_3 "AB"

* Outcome for heterogeneity regressions
global Y_HET "ln_incentive"

* Clustering toggle: worker_id / contract_id / wc_id
global CLUSTER_3 "worker_id"

* Controls toggles (mirror your Section 2 logic)
global USE_WEATHER_3   1
global USE_EXPER_3     1
global USE_HOURS_3     0
global USE_DAYFLAGS_3  0
global USE_LIFECYCLE_3 1

* Floors toggle (0 = none)
global USE_FLOORS_3   0
global MIN_HRS_3      2
global MIN_EARN_3     30
global USE_MIN_PROD_3 1
global MIN_PROD_3     50

* Event window (manual: keeps this section stable)
global WPRE_3  3
global WPOST_3 8


***********************************************************************
* 3.2) BUILD CONTROLS + FLOORS STRINGS + VCE STRING
***********************************************************************
global X_3 ""
if $USE_HOURS_3==1    global X_3 "$X_3 tot_hrs_day"
if $USE_DAYFLAGS_3==1 global X_3 "$X_3 multi_contract_day"
if $USE_EXPER_3==1    global X_3 "$X_3 c.experience##c.experience"
if $USE_WEATHER_3==1  global X_3 "$X_3 c.temp##c.temp c.wind_spd##c.wind_spd precip0 c.precip_pos##c.precip_pos"

global LIFECYCLE_3 ""
if $USE_LIFECYCLE_3==1 {
    global LIFECYCLE_3 "$LIFECYCLE_3 c.contract_day##c.contract_day"
    global LIFECYCLE_3 "$LIFECYCLE_3 c.season_days_worked##c.season_days_worked"
    global LIFECYCLE_3 "$LIFECYCLE_3 c.season_days_worked##c.cum_seasons"
}

global SAMPLE_OK_3 "1==1"
if $USE_FLOORS_3==1 {
    global SAMPLE_OK_3 "$SAMPLE_OK_3 & tot_hrs_day >= $MIN_HRS_3"
    global SAMPLE_OK_3 "$SAMPLE_OK_3 & piece > 0 & prod > 0"
    global SAMPLE_OK_3 "$SAMPLE_OK_3 & (piece*prod) >= $MIN_EARN_3"
    if $USE_MIN_PROD_3==1 global SAMPLE_OK_3 "$SAMPLE_OK_3 & prod >= $MIN_PROD_3"
}

* VCE string used by reghdfe (single place to toggle)
global VCE_3 "vce(cluster $CLUSTER_3)"

di as text "[3] TY_3:         $TY_3"
di as text "[3] CONTROL_3:    $CONTROL_3"
di as text "[3] Y_HET:        $Y_HET"
di as text "[3] CLUSTER:      $CLUSTER_3"
di as text "[3] X_3:          $X_3"
di as text "[3] LIFECYCLE_3:  $LIFECYCLE_3"
di as text "[3] FLOORS:       $SAMPLE_OK_3"
di as text "[3] Window:       [-$WPRE_3,$WPOST_3]"


***********************************************************************
* 3.3) DEFINE 2A-LIKE SAMPLE (BC vs AB in treated year) + EVENT TIME
***********************************************************************
gen byte sample_3 = (year==$TY_3) & inlist(province_s,"BC","$CONTROL_3")

gen byte treated_BC_3 = .
replace treated_BC_3 = (province_s=="BC") if sample_3==1

local tdate_3 = mdy(6,1,$TY_3)
gen byte post_3       = (date >= `tdate_3') if sample_3==1
gen int  event_week_3 = floor((date - `tdate_3')/7) if sample_3==1

tab province_s if sample_3==1
tab post_3     if sample_3==1


/**********************************************************************
  3.4) PRE-TREATMENT ABILITY (BC+AB, PRE ONLY)
  - ln(prod) residualized on controls + lifecycle + contract FE
  - estimated in PRE only (post_3==0)
  - ability(worker) = mean residual, propagated to all rows of worker
**********************************************************************/
gen double ln_prod = ln(prod) if prod>0

gen byte abil_sample = (sample_3==1) ///
    & (post_3==0) ///
    & !missing(ln_prod) ///
    & ($SAMPLE_OK_3)

di as text "[3.4] Ability PRE sample size:"
count if abil_sample==1
tab province_s if abil_sample==1

cap drop resid_preprod
reghdfe ln_prod $X_3 $LIFECYCLE_3 if abil_sample==1, ///
    absorb(contract_id) ///
    resid(resid_preprod)

bysort worker_id: egen double ability_raw = mean(resid_preprod)

quietly summarize ability_raw if !missing(ability_raw)
scalar abil_mu = r(mean)
scalar abil_sd = r(sd)

cap drop ability_z
gen double ability_z = (ability_raw - abil_mu)/abil_sd ///
    if !missing(ability_raw) & abil_sd>0

di as text "[3.4] Ability coverage (should be >0 in both groups):"
count if sample_3==1 & !missing(ability_z)
tab province_s if sample_3==1 & !missing(ability_z)
tab post_3 if sample_3==1, summarize(ability_z)
bys worker_id: assert ability_z==ability_z[1] if !missing(ability_z)


/**********************************************************************
  3.5) TOP-UP RISK: P(topup=1 | X, ability, contract difficulty, weather)
  - Define topup indicator using top_up (already in your data)
  - Fit PRE only (post_3==0) so we don't "train on treatment"
  - Include i.contract_id to capture inherent difficulty
  - Predict phat_topup for ALL observations in sample_3
**********************************************************************/
cap drop topup_ind topup_pre phat_topup
gen byte topup_ind = (top_up > 0) if !missing(top_up)

gen byte topup_pre = (sample_3==1) ///
    & (post_3==0) ///
    & !missing(topup_ind) ///
    & !missing(ability_z) ///
    & ($SAMPLE_OK_3)

di as text "[3.5] Top-up PRE sample size:"
count if topup_pre==1
tab province_s if topup_pre==1
summ topup_ind if topup_pre==1

* LPM is robust/safe; avoids separation issues in logit with FE
reg topup_ind ///
    c.ability_z ///
    $X_3 ///
    i.contract_id ///
    if topup_pre==1, ///
    vce(cluster $CLUSTER_3)

predict double phat_topup if sample_3==1, xb
replace phat_topup = max(0, min(1, phat_topup)) if sample_3==1 & !missing(phat_topup)

di as text "[3.5] Predicted top-up risk coverage:"
summ phat_topup if sample_3==1
tab post_3 if sample_3==1, summarize(phat_topup)


/**********************************************************************
  3.6) JOB DIFFICULTY: contract-level difficulty index (PRE only)
  Idea:
    - residualize ln_prod on worker ability + controls + lifecycle
      absorbing worker_id, using PRE only
    - contract difficulty = mean residual by contract in PRE
      (higher residual => easier contract conditional on worker+weather)
**********************************************************************/
cap drop resid_cd contract_diff_raw contract_diff_z

reghdfe ln_prod ///
    c.ability_z ///
    $X_3 ///
    $LIFECYCLE_3 ///
    if abil_sample==1 & !missing(ability_z), ///
    absorb(worker_id) ///
    resid(resid_cd)

* contract difficulty computed from PRE residuals; propagate to all rows
bysort contract_id: egen double contract_diff_raw = mean(resid_cd)

quietly summarize contract_diff_raw if sample_3==1 & !missing(contract_diff_raw)
scalar cd_mu = r(mean)
scalar cd_sd = r(sd)

cap drop contract_diff_z
gen double contract_diff_z = (contract_diff_raw - cd_mu)/cd_sd ///
    if sample_3==1 & !missing(contract_diff_raw) & cd_sd>0

di as text "[3.6] Contract difficulty coverage:"
count if sample_3==1 & !missing(contract_diff_z)
summ contract_diff_z if sample_3==1


/**********************************************************************
  3.7) DAILY DIFFICULTY: weather-driven "difficulty shock" (PRE only)
  Idea:
    - isolate day-to-day weather component predicting ln_prod
    - wx_shock is predicted weather effect (higher => easier weather day)
**********************************************************************/
cap drop resid_wx wx_shock

reghdfe ln_prod ///
    c.ability_z ///
    $X_3 ///
    if abil_sample==1 & !missing(ability_z), ///
    absorb(worker_id contract_id) ///
    resid(resid_wx)

* residual includes everything not explained by FE+X; for a clean "weather only"
* we instead get wx_shock as the fitted value from X_3 alone:
predict double wx_shock if sample_3==1, xb

di as text "[3.7] Weather shock coverage:"
summ wx_shock if sample_3==1


/**********************************************************************
  3.8) MAIN HETEROGENEITY DiD: does MW effect vary with MW-binding risk?
  Key coefficient: c.phat_topup#1.treated_BC_3#1.post_3
**********************************************************************/
reghdfe $Y_HET ///
    i.treated_BC_3##i.post_3 ///
    c.phat_topup##i.treated_BC_3##i.post_3 ///
    $X_3 ///
    $LIFECYCLE_3 ///
    if sample_3==1 ///
       & inrange(event_week_3,-$WPRE_3,$WPOST_3) ///
       & !missing($Y_HET) ///
       & !missing(phat_topup) ///
       & ($SAMPLE_OK_3), ///
    absorb(worker_id contract_id) ///
    $VCE_3
estimates store het_topup


/**********************************************************************
  3.9) DIFFICULTY HETEROGENEITY: contract difficulty + daily weather
  - contract_diff_z is time-invariant within contract (captured by contract FE),
    so its LEVEL will be collinear, but interactions with post/treat still work.
  - wx_shock varies daily, so it can shift effects within contract.
**********************************************************************/
reghdfe $Y_HET ///
    i.treated_BC_3##i.post_3 ///
    c.contract_diff_z##i.treated_BC_3##i.post_3 ///
    c.wx_shock##i.treated_BC_3##i.post_3 ///
    $X_3 ///
    $LIFECYCLE_3 ///
    if sample_3==1 ///
       & inrange(event_week_3,-$WPRE_3,$WPOST_3) ///
       & !missing($Y_HET) ///
       & !missing(contract_diff_z) ///
       & !missing(wx_shock) ///
       & ($SAMPLE_OK_3), ///
    absorb(worker_id contract_id) ///
    $VCE_3
estimates store het_difficulty


/**********************************************************************
  3.10) NONLINEAR ABILITY HETEROGENEITY (BINS) — quick diagnostic
**********************************************************************/
cap drop abil_bin_all
xtile abil_bin_all = ability_z if sample_3==1 & !missing(ability_z), nq(4)

reghdfe $Y_HET ///
    i.treated_BC_3##i.post_3 ///
    i.abil_bin_all##i.treated_BC_3##i.post_3 ///
    $X_3 ///
    $LIFECYCLE_3 ///
    if sample_3==1 ///
       & inrange(event_week_3,-$WPRE_3,$WPOST_3) ///
       & !missing($Y_HET) ///
       & !missing(ability_z) ///
       & ($SAMPLE_OK_3), ///
    absorb(worker_id contract_id) ///
    $VCE_3
estimates store het_bins


/**********************************************************************
  3.11) DESCRIPTIVE TABLE: ability bins × pre/post (BC only)
**********************************************************************/
table abil_bin_all post_3 if (year==$TY_3) & (province_s=="BC") ///
    & inrange(event_week_3,-$WPRE_3,$WPOST_3), ///
    statistic(mean piece) statistic(mean prod) statistic(mean tot_hrs_day) nformat(%9.3f)
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
/*	
/***********************************************************************
  (3) ADDITIONAL EMPIRICAL ANALYSIS (SELF-CONTAINED)

  PURPOSE:
    - Pre-treatment ability construction
    - Heterogeneous MW effects (BC vs AB, treated year)
    - Nonlinear heterogeneity checks (quadratic + bins)
    - Selection diagnostic
    - Piece-rate sorting tests (linear + nonlinear)

  Requires:
    - Section 1 data prep already run
***********************************************************************/

***********************************************************************
* 3.0) CLEAN RERUN STATE
***********************************************************************
cap estimates clear
foreach v in sample_3 treated_BC_3 post_3 event_week_3 ///
            ln_prod resid_all ability_raw ability_z ability_z2 ///
            abil_sample abil_bin_all ln_piece {
    cap drop `v'
}
foreach m in TY_3 CONTROL_3 CLUSTER_3 VCE_3 ///
            USE_FLOORS_3 MIN_HRS_3 MIN_EARN_3 MIN_PROD_3 USE_MIN_PROD_3 ///
            USE_WEATHER_3 USE_EXPER_3 USE_HOURS_3 USE_DAYFLAGS_3 USE_LIFECYCLE_3 ///
            WPRE_3 WPOST_3 ///
            X_3 LIFECYCLE_3 SAMPLE_OK_3 Y_HET {
    cap macro drop `m'
}

***********************************************************************
* 3.1) SETTINGS
***********************************************************************
global TY_3 2018
global CONTROL_3 "AB"
global Y_HET "ln_incentive"
global CLUSTER_3 "worker_id"

* Controls
global USE_WEATHER_3   1
global USE_EXPER_3     1
global USE_HOURS_3     0
global USE_DAYFLAGS_3  0
global USE_LIFECYCLE_3 1

* Floors
global USE_FLOORS_3   0
global MIN_HRS_3      2
global MIN_EARN_3     30
global USE_MIN_PROD_3 1
global MIN_PROD_3     50

* Event window
global WPRE_3  3
global WPOST_3 8

global VCE_3 "vce(cluster $CLUSTER_3)"

***********************************************************************
* 3.2) BUILD CONTROL STRINGS
***********************************************************************
global X_3 ""
if $USE_HOURS_3==1    global X_3 "$X_3 tot_hrs_day"
if $USE_DAYFLAGS_3==1 global X_3 "$X_3 multi_contract_day"
if $USE_EXPER_3==1    global X_3 "$X_3 c.experience##c.experience"
if $USE_WEATHER_3==1  global X_3 "$X_3 c.temp##c.temp c.wind_spd##c.wind_spd precip0 c.precip_pos##c.precip_pos"

global LIFECYCLE_3 ""
if $USE_LIFECYCLE_3==1 {
    global LIFECYCLE_3 "$LIFECYCLE_3 c.contract_day##c.contract_day"
    global LIFECYCLE_3 "$LIFECYCLE_3 c.season_days_worked##c.season_days_worked"
    global LIFECYCLE_3 "$LIFECYCLE_3 c.season_days_worked##c.cum_seasons"
}

global SAMPLE_OK_3 "1==1"
if $USE_FLOORS_3==1 {
    global SAMPLE_OK_3 "$SAMPLE_OK_3 & tot_hrs_day >= $MIN_HRS_3"
    global SAMPLE_OK_3 "$SAMPLE_OK_3 & piece > 0 & prod > 0"
    global SAMPLE_OK_3 "$SAMPLE_OK_3 & (piece*prod) >= $MIN_EARN_3"
    if $USE_MIN_PROD_3==1 global SAMPLE_OK_3 "$SAMPLE_OK_3 & prod >= $MIN_PROD_3"
}

***********************************************************************
* 3.3) DEFINE SAMPLE (BC vs AB, TREATED YEAR)
***********************************************************************
gen byte sample_3 = (year==$TY_3) & inlist(province_s,"BC","$CONTROL_3")
gen byte treated_BC_3 = (province_s=="BC") if sample_3==1

local tdate_3 = mdy(6,1,$TY_3)
gen byte post_3       = (date >= `tdate_3') if sample_3==1
gen int  event_week_3 = floor((date - `tdate_3')/7) if sample_3==1

***********************************************************************
* 3.4) PRE-TREATMENT ABILITY (BC + AB, PRE-JUNE ONLY)
***********************************************************************
gen double ln_prod = ln(prod) if prod>0

gen byte abil_sample = sample_3==1 ///
    & post_3==0 ///
    & !missing(ln_prod) ///
    & ($SAMPLE_OK_3)

reghdfe ln_prod ///
    $X_3 ///
    $LIFECYCLE_3 ///
    if abil_sample==1, ///
    absorb(contract_id) ///
    resid(resid_all)

bys worker_id: egen double ability_raw = mean(resid_all)

quietly summarize ability_raw if !missing(ability_raw)
scalar abil_mu = r(mean)
scalar abil_sd = r(sd)

gen double ability_z = (ability_raw - abil_mu)/abil_sd ///
    if !missing(ability_raw) & abil_sd>0

gen double ability_z2 = ability_z^2 if !missing(ability_z)

***********************************************************************
* 3.5) LINEAR HETEROGENEITY DiD
***********************************************************************
reghdfe $Y_HET ///
    i.treated_BC_3##i.post_3 ///
    c.ability_z##i.treated_BC_3##i.post_3 ///
    $X_3 ///
    $LIFECYCLE_3 ///
    if sample_3==1 ///
       & inrange(event_week_3,-$WPRE_3,$WPOST_3) ///
       & !missing($Y_HET) ///
       & !missing(ability_z) ///
       & ($SAMPLE_OK_3), ///
    absorb(worker_id contract_id) ///
    $VCE_3
estimates store het_linear

***********************************************************************
* 3.5B) QUADRATIC (U-SHAPE) HETEROGENEITY
***********************************************************************
reghdfe $Y_HET ///
    i.treated_BC_3##i.post_3 ///
    c.ability_z##i.treated_BC_3##i.post_3 ///
    c.ability_z2##i.treated_BC_3##i.post_3 ///
    $X_3 ///
    $LIFECYCLE_3 ///
    if sample_3==1 ///
       & inrange(event_week_3,-$WPRE_3,$WPOST_3) ///
       & !missing($Y_HET) ///
       & !missing(ability_z) ///
       & ($SAMPLE_OK_3), ///
    absorb(worker_id contract_id) ///
    $VCE_3
estimates store het_quad

quietly test ///
    1.treated_BC_3#1.post_3#c.ability_z ///
    1.treated_BC_3#1.post_3#c.ability_z2
di as text "[3.5B] Nonlinear (quadratic) heterogeneity p = " %6.4f r(p)

***********************************************************************
* 3.5C) ABILITY BIN HETEROGENEITY
***********************************************************************
xtile abil_bin_all = ability_z if sample_3==1 & !missing(ability_z), nq(4)

reghdfe $Y_HET ///
    i.treated_BC_3##i.post_3 ///
    i.abil_bin_all##i.treated_BC_3##i.post_3 ///
    $X_3 ///
    $LIFECYCLE_3 ///
    if sample_3==1 ///
       & inrange(event_week_3,-$WPRE_3,$WPOST_3) ///
       & !missing($Y_HET) ///
       & !missing(ability_z) ///
       & ($SAMPLE_OK_3), ///
    absorb(worker_id contract_id) ///
    $VCE_3
estimates store het_bins

quietly test ///
    2.abil_bin_all#1.treated_BC_3#1.post_3 ///
    3.abil_bin_all#1.treated_BC_3#1.post_3 ///
    4.abil_bin_all#1.treated_BC_3#1.post_3
di as text "[3.5C] Bin heterogeneity joint p = " %6.4f r(p)

***********************************************************************
* 3.6) SELECTION DIAGNOSTIC (BC ONLY)
***********************************************************************
preserve
    keep if year==$TY_3 & province_s=="BC"
    bys worker_id: egen byte any_pre  = max(post_3==0)
    bys worker_id: egen byte any_post = max(post_3==1)
    bys worker_id: keep if _n==1
    keep if any_pre==1 & !missing(ability_z)
    reg any_post c.ability_z, vce(cluster worker_id)
restore

***********************************************************************
* 3.7) PIECE-RATE SORTING TEST (BC ONLY)
***********************************************************************
gen double ln_piece = ln(piece) if piece>0

reghdfe ln_piece ///
    i.post_3 ///
    c.ability_z##i.post_3 ///
    $X_3 ///
    if year==$TY_3 ///
       & province_s=="BC" ///
       & inrange(event_week_3,-$WPRE_3,$WPOST_3) ///
       & !missing(ln_piece) ///
       & !missing(ability_z), ///
    absorb(worker_id contract_id) ///
    $VCE_3
estimates store sort_linear

***********************************************************************
* 3.7B) NONLINEAR PIECE-RATE TEST
***********************************************************************
reghdfe ln_piece ///
    i.post_3 ///
    c.ability_z##i.post_3 ///
    c.ability_z2##i.post_3 ///
    $X_3 ///
    if year==$TY_3 ///
       & province_s=="BC" ///
       & inrange(event_week_3,-$WPRE_3,$WPOST_3) ///
       & !missing(ln_piece) ///
       & !missing(ability_z), ///
    absorb(worker_id contract_id) ///
    $VCE_3
estimates store sort_quad

quietly test 1.post_3#c.ability_z 1.post_3#c.ability_z2
di as text "[3.7B] Piece-rate nonlinear p = " %6.4f r(p)


/**********************************************************************
  3.8) DESCRIPTIVE ABILITY BINS (BC ONLY)
**********************************************************************/

cap drop abil_bin_bc
xtile abil_bin_bc = ability_z ///
    if year==$TY_3 & province_s=="BC" & !missing(ability_z), nq(4)

table abil_bin_bc post_3 ///
    if year==$TY_3 & province_s=="BC" ///
       & inrange(event_week_3,-$WPRE_3,$WPOST_3), ///
    statistic(mean piece) ///
    statistic(mean prod) ///
    statistic(mean tot_hrs_day)
	
	
	
	
	
	
	
	
	
	
	
	/*
/**********************************************************************
  (3) HETEROGENEOUS EFFECTS USING AB AS CONTROL (SELF-CONTAINED)

  What this does:
    3.1 Define 2A-style sample (BC vs AB in treated year) + event time
    3.2 Construct ability for BOTH BC and AB (NO pre-only restriction)
    3.3 Run heterogeneity DiD with worker+contract FE:
          Y = treated×post + ability×post + ability×treated×post + controls + FE

  Notes:
    - Ability is constructed from ln(prod) residualized on controls + contract FE.
    - Uses all observations in the 2A sample (as requested).
**********************************************************************/

***********************************************************************
* 3.0) CLEAN RERUN STATE (SAFE)
***********************************************************************
cap estimates clear
foreach v in sample_3 treated_BC_3 post_3 event_week_3 ///
            ln_prod resid_all ability_raw ability_z {
    cap drop `v'
}
foreach m in TY_3 CONTROL_3 CLUSTER_3 ///
            USE_FLOORS_3 MIN_HRS_3 MIN_EARN_3 MIN_PROD_3 USE_MIN_PROD_3 ///
            USE_WEATHER_3 USE_EXPER_3 USE_HOURS_3 USE_DAYFLAGS_3 USE_LIFECYCLE_3 ///
            TAU_3 MAXW_3 CORE0_3 CORE1_3 AUTO_WIN_3 WPRE_3 WPOST_3 ///
            X_3 LIFECYCLE_3 SAMPLE_OK_3 Y_HET {
    cap macro drop `m'
}

***********************************************************************
* 3.1) SETTINGS (EDIT THESE)
***********************************************************************
* Treated year and control province (AB control)
global TY_3 2019
global CONTROL_3 "AB"

* Outcome for heterogeneity (default to baseline outcome)
global Y_HET "ln_incentive"

* Clustering toggle (IMPORTANT with few contracts)
* Options: worker_id, contract_id, wc_id
global CLUSTER_3 "worker_id"

* Controls toggles (mirror Section 2 style)
global USE_WEATHER_3   1
global USE_EXPER_3     1
global USE_HOURS_3     0
global USE_DAYFLAGS_3  0
global USE_LIFECYCLE_3 1

* Floors toggle (0 = none)
global USE_FLOORS_3   0
global MIN_HRS_3      2
global MIN_EARN_3     30
global USE_MIN_PROD_3 1
global MIN_PROD_3     50

* Event window (manual here to stay robust / no extra moving parts)
global WPRE_3  3
global WPOST_3 8

***********************************************************************
* 3.2) BUILD CONTROLS + FLOORS STRINGS
***********************************************************************
global X_3 ""
if $USE_HOURS_3==1    global X_3 "$X_3 tot_hrs_day"
if $USE_DAYFLAGS_3==1 global X_3 "$X_3 multi_contract_day"
if $USE_EXPER_3==1    global X_3 "$X_3 c.experience##c.experience"
if $USE_WEATHER_3==1  global X_3 "$X_3 c.temp##c.temp c.wind_spd##c.wind_spd precip0 c.precip_pos##c.precip_pos"

global LIFECYCLE_3 ""
if $USE_LIFECYCLE_3==1 {
    global LIFECYCLE_3 "$LIFECYCLE_3 c.contract_day##c.contract_day"
    global LIFECYCLE_3 "$LIFECYCLE_3 c.season_days_worked##c.season_days_worked"
    global LIFECYCLE_3 "$LIFECYCLE_3 c.season_days_worked##c.cum_seasons"
}

global SAMPLE_OK_3 "1==1"
if $USE_FLOORS_3==1 {
    global SAMPLE_OK_3 "$SAMPLE_OK_3 & tot_hrs_day >= $MIN_HRS_3"
    global SAMPLE_OK_3 "$SAMPLE_OK_3 & piece > 0 & prod > 0"
    global SAMPLE_OK_3 "$SAMPLE_OK_3 & (piece*prod) >= $MIN_EARN_3"
    if $USE_MIN_PROD_3==1 global SAMPLE_OK_3 "$SAMPLE_OK_3 & prod >= $MIN_PROD_3"
}

di as text "[3] Y_HET:       $Y_HET"
di as text "[3] CLUSTER:     $CLUSTER_3"
di as text "[3] X_3:         $X_3"
di as text "[3] LIFECYCLE_3: $LIFECYCLE_3"
di as text "[3] FLOORS:      $SAMPLE_OK_3"
di as text "[3] Window:      [-$WPRE_3,$WPOST_3]"

***********************************************************************
* 3.3) DEFINE 2A-LIKE SAMPLE (BC vs AB IN TREATED YEAR) + EVENT TIME
***********************************************************************
gen byte sample_3 = (year==$TY_3) & inlist(province_s,"BC","$CONTROL_3")

gen byte treated_BC_3 = .
replace treated_BC_3 = (province_s=="BC") if sample_3==1

local tdate_3 = mdy(6,1,$TY_3)
gen byte post_3       = (date >= `tdate_3') if sample_3==1
gen int  event_week_3 = floor((date - `tdate_3')/7) if sample_3==1

tab province_s if sample_3==1
tab treated_BC_3 if sample_3==1
tab post_3 if sample_3==1

/**********************************************************************
  3.4) CONSTRUCT PRE-TREATMENT ABILITY (BC + AB, PRE-JUNE ONLY)

  Ability definition:
    - ln(prod) residualized on:
        * controls ($X_3)
        * lifecycle ($LIFECYCLE_3)
        * contract FE
    - estimated using PRE-TREATMENT observations only
    - ability(worker) = mean residual across pre-period days
    - standardized across workers
**********************************************************************/

*-------------------------------
* 3.4.0 Clean old objects
*-------------------------------
cap drop ln_prod resid_all ability_raw ability_z abil_sample

*-------------------------------
* 3.4.1 Define ln(prod)
*-------------------------------
gen double ln_prod = ln(prod) if prod>0

*-------------------------------
* 3.4.2 Define PRE-TREATMENT ability sample
*   - sample_3 ensures BC + AB in treated year
*   - post_3==0 ensures strictly pre-June 1
*-------------------------------
cap drop abil_sample
gen byte abil_sample = (sample_3==1) ///
    & (post_3==0) ///
    & !missing(ln_prod) ///
    & ($SAMPLE_OK_3)

di as text "[3.4] Pre-treatment ability sample:"
count if abil_sample==1
tab province_s if abil_sample==1
summ ln_prod if abil_sample==1

*-------------------------------
* 3.4.3 Residualize productivity in PRE period only
*-------------------------------
reghdfe ln_prod ///
    $X_3 ///
    $LIFECYCLE_3 ///
    if abil_sample==1, ///
    absorb(contract_id) ///
    resid(resid_all)

* Confirm residual coverage
count if abil_sample==1 & !missing(resid_all)

*-------------------------------
* 3.4.4 Worker-level ability
*   egen mean() propagates ability to all rows of that worker
*-------------------------------
bysort worker_id: egen double ability_raw = mean(resid_all)

*-------------------------------
* 3.4.5 Standardize across workers
*-------------------------------
quietly summarize ability_raw if !missing(ability_raw)
scalar abil_mu = r(mean)
scalar abil_sd = r(sd)

gen double ability_z = (ability_raw - abil_mu)/abil_sd ///
    if !missing(ability_raw) & abil_sd>0

*-------------------------------
* 3.4.6 Diagnostics
*-------------------------------
di as text "[3.4] Ability coverage diagnostics:"
count if sample_3==1 & !missing(ability_z)
tab province_s if sample_3==1 & !missing(ability_z)
tab post_3      if sample_3==1, summarize(ability_z)

* Ability must be time-invariant within worker
bys worker_id: assert ability_z==ability_z[1] if !missing(ability_z)

di as text "[3.4] DONE: Pre-treatment ability constructed."

***********************************************************************
* 3.5) HETEROGENEITY DiD (WORKER+CONTRACT FE)
*
* Key coefficient: c.ability_z#1.treated_BC_3#1.post_3
* Interpretation: how the BC-vs-AB post effect changes with ability.
***********************************************************************
reghdfe $Y_HET ///
    i.treated_BC_3##i.post_3 ///
    c.ability_z##i.treated_BC_3##i.post_3 ///
    $X_3 ///
    $LIFECYCLE_3 ///
    if sample_3==1 ///
       & inrange(event_week_3, -$WPRE_3, $WPOST_3) ///
       & !missing($Y_HET) ///
       & !missing(ability_z) ///
       & ($SAMPLE_OK_3), ///
    absorb(worker_id contract_id) ///
    vce(cluster $CLUSTER_3)

estimates store het_2A_AB

* Quick support diagnostics
tab post_3 if e(sample)
tab treated_BC_3 if e(sample)

* How much sample do we have in the heterogeneity window BEFORE ability?
count if sample_3==1 ///
    & inrange(event_week_3, -$WPRE_3, $WPOST_3) ///
    & !missing($Y_HET) ///
    & ($SAMPLE_OK_3)

tab province_s if sample_3==1 ///
    & inrange(event_week_3, -$WPRE_3, $WPOST_3) ///
    & !missing($Y_HET) ///
    & ($SAMPLE_OK_3)

* Now: how much survives once you require ability_z?
count if sample_3==1 ///
    & inrange(event_week_3, -$WPRE_3, $WPOST_3) ///
    & !missing($Y_HET) ///
    & !missing(ability_z) ///
    & ($SAMPLE_OK_3)

tab province_s if sample_3==1 ///
    & inrange(event_week_3, -$WPRE_3, $WPOST_3) ///
    & !missing($Y_HET) ///
    & !missing(ability_z) ///
    & ($SAMPLE_OK_3)
	
	count if sample_3==1 & !missing(ln_prod)

count if sample_3==1 ///
    & !missing(ln_prod) ///
    & inrange(event_week_3, -$WPRE_3, $WPOST_3)

summ ln_prod if sample_3==1

	
/**********************************************************************
  3.6) INTERPRETATION HELPERS (COMPUTE IMPLIED EFFECTS)
  - Baseline DiD effect for BC at ability=0 is coefficient on treat_post.
  - Heterogeneity: coefficient on ability_treat_post tells how the DiD
    effect changes per 1 SD of ability (within BC).
**********************************************************************/
di as text "------------------------------------------------------------"
di as text "[3.6] Key coefficients from het_2A:"
di as text "  - treat_post        : DiD effect at ability=0 (BC vs control)"
di as text "  - ability_treat_post: extra DiD effect per 1 SD higher ability (BC only)"
di as text "------------------------------------------------------------"


/**********************************************************************
  3.7) SELECTION DIAGNOSTIC (WITHIN BC, TREATED YEAR ONLY)
  Question:
    - Among BC workers observed pre, are higher-ability workers more likely
      to appear post within the treated year?
  Note:
    - This is NOT a causal "MW selection effect"; it's a descriptive check.
**********************************************************************/

preserve
    keep if (year==$TY_3) & (province_s=="BC")
    keep worker_id post_3 ability_z

    bysort worker_id: egen byte any_pre  = max(post_3==0)
    bysort worker_id: egen byte any_post = max(post_3==1)
    bysort worker_id: keep if _n==1

    keep if any_pre==1 & !missing(ability_z)

    reg any_post c.ability_z, vce(cluster worker_id)
restore


/**********************************************************************
  3.8) SORTING / MATCHING DIAGNOSTIC (BC ONLY, TREATED YEAR)
  Question:
    - Does the relationship between ability and piece rate differ post vs pre?
  Note:
    - Not causal. It's about assignment/matching changes.
**********************************************************************/

cap drop ln_piece
gen double ln_piece = ln(piece) if piece>0

reghdfe ln_piece ///
    post_3 c.ability_z##i.post_3 ///
    $X_3 ///
    if (year==$TY_3) ///
       & (province_s=="BC") ///
       & inrange(event_week_3, -$WPRE_3, $WPOST_3) ///
       & ($SAMPLE_OK_3) ///
       & !missing(ln_piece) ///
       & !missing(ability_z), ///
    absorb(worker_id contract_id) ///
    $VCE_3

estimates store sort_BC


/**********************************************************************
  3.9) SIMPLE DESCRIPTIVE TABLE: ABILITY BINS × PRE/POST (BC ONLY)
**********************************************************************/

cap drop abil_bin
xtile abil_bin = ability_z if (year==$TY_3) & (province_s=="BC") & !missing(ability_z), nq(4)

table abil_bin post_3 if (year==$TY_3) & (province_s=="BC") ///
    & inrange(event_week_3, -$WPRE_3, $WPOST_3), ///
    statistic(mean piece) statistic(mean prod) statistic(mean tot_hrs_day) nformat(%9.3f)
	
	
	/*
	

/**********************************************************************
  (3) ADDITIONAL EMPIRICAL ANALYSIS — CONSISTENT WITH SECTION 2A

  PURPOSE:
    - Construct predetermined worker ability using BC-only pre-June 1 data
      (treated year, same 2A window and sample filters)
    - Heterogeneous MW effects in 2A design (BC vs AB control)
    - Selection diagnostic within BC (who shows up post?)
    - Sorting/matching diagnostic: does ability–piece relationship shift post?

  REQUIREMENTS:
    - Section (1) DATA PREP already run
    - Section (2) BASELINE ANALYSIS already run THROUGH 2A.3 so these exist:
        sample_2A, treated_BC_2A, post_2A, event_week_2A
        globals: Y, X, LIFECYCLE, ABS, SAMPLE_OK, WPRE_2A, WPOST_2A
**********************************************************************/


/**********************************************************************
  3.0) FORCE CLEAR (SAFE RERUN OF SECTION 3 ONLY)
**********************************************************************/
cap estimates clear

foreach v in ///
    ln_prod abil_sample resid_pre ability_raw ability_z ///
    treat_post ability_post ability_treat ability_treat_post ///
    any_pre any_post ///
    ln_piece {
    cap drop `v'
}

foreach m in ///
    PRE_ABIL Y_HET CLUSTVAR {
    cap macro drop `m'
}


/**********************************************************************
  3.1) SETTINGS (EDIT THESE)
**********************************************************************/

* Ability construction window (pre only, within the 2A event-week window)
* Example: use weeks -12..-1 if available; will still be restricted by 2A window.
global PRE_ABIL -12

* Outcome for heterogeneity (default = same as baseline $Y)
capture confirm global Y_HET
if _rc global Y_HET "$Y"

* Cluster choice toggle:
*   contract_id is often too few clusters in your data (e.g., 29).
*   worker_id is usually safer.
global CLUSTVAR "worker_id"    // "worker_id" or "contract_id" (or "wc_id" if you prefer)


/**********************************************************************
  3.2) DEFENSIVE CHECKS (DON'T PROCEED IF SECTION 2A NOT RUN)
**********************************************************************/
capture confirm variable sample_2A
if _rc {
    di as error "[Section 3] sample_2A not found. Run Section 2A first."
    exit 198
}
capture confirm variable treated_BC_2A
if _rc {
    di as error "[Section 3] treated_BC_2A not found. Run Section 2A first."
    exit 198
}
capture confirm variable post_2A
if _rc {
    di as error "[Section 3] post_2A not found. Run Section 2A first."
    exit 198
}
capture confirm variable event_week_2A
if _rc {
    di as error "[Section 3] event_week_2A not found. Run Section 2A first."
    exit 198
}
capture confirm global WPRE_2A
if _rc {
    di as error "[Section 3] global WPRE_2A not found. Run Section 2A first."
    exit 198
}
capture confirm global WPOST_2A
if _rc {
    di as error "[Section 3] global WPOST_2A not found. Run Section 2A first."
    exit 198
}

di as text "[Section 3] Using:"
di as text "  Y_HET     = $Y_HET"
di as text "  Cluster   = $CLUSTVAR"
di as text "  Window    = [-$WPRE_2A,$WPOST_2A] (from Section 2A)"
di as text "  Floors    = $SAMPLE_OK"
di as text "  Controls  = $X"
di as text "  Lifecycle = $LIFECYCLE"
di as text "  FE        = $ABS"


/**********************************************************************
  3.3) CONSTRUCT PREDETERMINED ABILITY (BC-ONLY, PRE-JUNE 1, TREATED YEAR)

  Approach:
    - Use ln(prod) (lprod) as performance outcome
    - Residualize ln(prod) on the SAME control set ($X + $LIFECYCLE) + contract FE
      using ONLY BC pre observations in the 2A sample
    - Ability(worker) = mean residual across those pre observations
    - Standardize over BC workers with ability defined
    - Set ability_z = 0 for AB (so we do NOT "estimate ability for AB")
**********************************************************************/

* ln(prod) exists as lprod from Section 1
gen double ln_prod = lprod

* Ability sample: BC only, pre only, within ability window AND within 2A analysis window
gen byte abil_sample = (sample_2A==1) ///
    & (treated_BC_2A==1) ///
    & (post_2A==0) ///
    & inrange(event_week_2A, $PRE_ABIL, -1) ///
    & inrange(event_week_2A, -$WPRE_2A, $WPOST_2A) ///
    & ($SAMPLE_OK) ///
    & !missing(ln_prod)

di as text "[3.3] Ability sample counts"
count if abil_sample==1
tab event_week_2A if abil_sample==1

* Residualize: controls + lifecycle + contract FE (difficulty proxy)
* IMPORTANT: include resid option so predict works.
reghdfe ln_prod ///
    $X ///
    $LIFECYCLE ///
    if abil_sample==1, ///
    absorb(contract_id) ///
    vce(cluster $CLUSTVAR) ///
    resid

cap drop resid_pre
gen double resid_pre = _reghdfe_resid if abil_sample==1

* Worker-level ability = mean residual across pre observations
bys worker_id: egen double ability_raw = mean(resid_pre)

* Standardize over BC workers with ability defined
summ ability_raw if !missing(ability_raw), meanonly
cap drop ability_z
gen double ability_z = (ability_raw - r(mean)) / r(sd) if !missing(ability_raw)

* Set ability to zero for AB/control group (so heterogeneity is "within BC only")
replace ability_z = 0 if sample_2A==1 & treated_BC_2A==0

di as text "[3.3] Ability coverage checks (within 2A sample)"
count if sample_2A==1 & treated_BC_2A==1 & !missing(ability_raw)
count if sample_2A==1 & treated_BC_2A==1 & !missing(ability_z)
count if sample_2A==1 & treated_BC_2A==0 & ability_z==0


/**********************************************************************
  3.4) HETEROGENEOUS MW EFFECTS (2A DESIGN: BC vs AB, TREATED YEAR)

  Interpretation focus:
    - The core DiD effect is treat_post (BC × post)
    - The heterogeneity term of interest is ability_treat_post:
        does the BC×post effect vary with BC workers' pre ability?

  We build explicit interactions to keep it transparent and stable.
**********************************************************************/

cap drop treat_post ability_post ability_treat ability_treat_post
gen byte   treat_post        = (treated_BC_2A==1 & post_2A==1) if sample_2A==1
gen double ability_post      = ability_z * post_2A             if sample_2A==1
gen double ability_treat     = ability_z * treated_BC_2A       if sample_2A==1
gen double ability_treat_post= ability_z * treat_post          if sample_2A==1

* Main heterogeneity regression (same window/sample logic as Section 2A)
reghdfe $Y_HET ///
    treated_BC_2A post_2A treat_post ///
    ability_z ability_post ability_treat ability_treat_post ///
    $X ///
    $LIFECYCLE ///
    if sample_2A==1 ///
       & inrange(event_week_2A, -$WPRE_2A, $WPOST_2A) ///
       & ($SAMPLE_OK) ///
       & !missing($Y_HET) ///
       & !missing(ability_z), ///
    absorb($ABS) ///
    vce(cluster $CLUSTVAR)

estimates store het_2A

di as text "[3.4] Estimation-sample support checks"
tab post_2A if e(sample)
tab treated_BC_2A if e(sample)


/**********************************************************************
  3.5) SELECTION DIAGNOSTIC (WITHIN BC ONLY, TREATED YEAR)

  Question:
    - Among BC workers observed pre-June 1, do higher-ability workers
      have a different probability of being observed post-June 1?

  This is NOT a causal MW effect; it's checking compositional change.
**********************************************************************/

preserve
    keep if sample_2A==1 & treated_BC_2A==1
    keep worker_id post_2A event_week_2A ability_z $Y_HET

    * worker-level indicators
    bysort worker_id: egen byte any_pre  = max(post_2A==0)
    bysort worker_id: egen byte any_post = max(post_2A==1)

    bysort worker_id: keep if _n==1
    keep if any_pre==1

    reg any_post c.ability_z, vce(cluster worker_id)
restore


/**********************************************************************
  3.6) SORTING / MATCHING DIAGNOSTIC (PIECE RATE) — 2A SAMPLE

  Question:
    - Does the ability–piece relationship shift post in BC
      (relative to AB as a control environment)?

  Outcome:
    ln_piece = ln(piece)
  Key heterogeneity term:
    ability_treat_post_p  (BC×post shift in ability slope)
**********************************************************************/

cap drop ln_piece
gen double ln_piece = ln(piece) if piece>0

cap drop treat_post_p ability_post_p ability_treat_p ability_treat_post_p
gen byte   treat_post_p         = (treated_BC_2A==1 & post_2A==1) if sample_2A==1
gen double ability_post_p       = ability_z * post_2A             if sample_2A==1
gen double ability_treat_p      = ability_z * treated_BC_2A       if sample_2A==1
gen double ability_treat_post_p = ability_z * treat_post_p        if sample_2A==1

reghdfe ln_piece ///
    treated_BC_2A post_2A treat_post_p ///
    ability_z ability_post_p ability_treat_p ability_treat_post_p ///
    $X ///
    $LIFECYCLE ///
    if sample_2A==1 ///
       & inrange(event_week_2A, -$WPRE_2A, $WPOST_2A) ///
       & ($SAMPLE_OK) ///
       & !missing(ln_piece) ///
       & !missing(ability_z), ///
    absorb($ABS) ///
    vce(cluster $CLUSTVAR)

estimates store sort_2A

tab post_2A if e(sample)
tab treated_BC_2A if e(sample)


/**********************************************************************
  3.7) OPTIONAL: SIMPLE DESCRIPTIVE BINNING (BC ONLY)
**********************************************************************/

cap drop abil_bin
xtile abil_bin = ability_z if sample_2A==1 & treated_BC_2A==1 & !missing(ability_z), nq(4)

table abil_bin post_2A if sample_2A==1 & treated_BC_2A==1 ///
    & inrange(event_week_2A, -$WPRE_2A, $WPOST_2A), ///
    statistic(mean piece) statistic(mean prod) statistic(mean tot_hrs_day) nformat(%9.3f)
























/**********************************************************************
  MASTER DO-FILE: Planter MW Analysis

  STRUCTURE:
    (0) PACKAGES
    (1) DATA PREP (run once)
        1A) Import, types, IDs, outcomes
        1B) Safe string copies (avoid strL collapse/reshape failures)
        1C) Universal June 1 event-time (year-specific)
        1D) Contract lifecycle + within-season ramp + cross-season experience
    (2) BASELINE ANALYSIS (rerunnable control panel)
        2A) BC vs AB/ON (treated-year cross-province DiD + event study)
        2B) BC placebo-years (within-BC DiD + event study)

  FEATURES:
    - Floors to avoid tiny/partial-day artifacts
    - Auto event-window selection by TAU rule (2A and 2B separately)
      using support in the *actual regression sample*
    - Support-restricted event-study plots (coefplot uses only e(sample))

  REQUIREMENTS:
    - reghdfe (and ftools)
    - coefplot


clear all
set more off


/**********************************************************************
  (0) PACKAGES
**********************************************************************/
capture which reghdfe
if _rc {
    ssc install ftools, replace
    ssc install reghdfe, replace
}
capture which coefplot
if _rc {
    ssc install coefplot, replace
}


/**********************************************************************
  (1) DATA PREP (RUN ONCE)
**********************************************************************/

* --- Load data ---
import delimited using "/Users/danaandersen/Downloads/planter_mw.csv", clear varnames(1)

* --- Required vars check ---
confirm variable date year province contract name piece prod top_up tot_hrs
confirm variable temp total_precip wind_spd experience

* --- Date conversion to Stata daily ---
capture confirm numeric variable date
if _rc {
    gen double date_num = daily(date, "YMD")
    format date_num %td
    drop date
    rename date_num date
}
format date %td

* --- Destring core numeric vars if needed ---
foreach var in temp total_precip wind_spd experience tot_hrs piece prod top_up {
    capture confirm numeric variable `var'
    if _rc destring `var', replace ignore("NA" "." "" "null" "NULL")
}

* --- IDs (create once) ---
cap drop worker_id contract_id wc_id
encode name,     gen(worker_id)
encode contract, gen(contract_id)
egen wc_id = group(worker_id contract_id)

* --- Day-level aggregates to handle multi-contract days ---
cap drop tot_hrs_day n_contracts_day multi_contract_day long_day extreme_day
bysort worker_id date: egen double tot_hrs_day      = total(tot_hrs)
bysort worker_id date: egen int    n_contracts_day  = count(contract_id)
gen byte multi_contract_day = (n_contracts_day > 1)
gen byte long_day           = (tot_hrs_day > 14)
gen byte extreme_day        = (tot_hrs_day > 18)

* --- Outcomes (create once) ---
cap drop ln_incentive lprod ln_piece
gen double ln_incentive = ln(piece*prod) if piece>0 & prod>0
gen double lprod        = ln(prod)       if prod>0            // prod = trees/hour
gen double ln_piece     = ln(piece)      if piece>0

* --- Weather transforms (create once) ---
cap drop precip0 precip_pos temp2 wind2 precip2
gen byte   precip0    = (total_precip == 0) if !missing(total_precip)
gen double precip_pos = cond(total_precip > 0 & !missing(total_precip), total_precip, 0)
gen double temp2      = temp^2       if !missing(temp)
gen double wind2      = wind_spd^2   if !missing(wind_spd)
gen double precip2    = precip_pos^2 if !missing(precip_pos)

* --- Experience squared if not already in data ---
capture confirm variable experience2
if _rc gen double experience2 = experience^2 if !missing(experience)


/**********************************************************************
  (1B) SAFE STRING COPIES (AVOID strL ISSUES IN collapse/reshape)
**********************************************************************/
cap drop province_s contract_s name_s
gen str2   province_s = province
gen str40  contract_s = contract
gen str40  name_s     = name


/**********************************************************************
  (1C) UNIVERSAL YEAR-SPECIFIC JUNE 1 EVENT TIME (RUN ONCE)
**********************************************************************/
cap drop tdate_y post_y event_day_y event_week_y
gen double tdate_y       = mdy(6,1,year)
format tdate_y %td
gen byte   post_y        = (date >= tdate_y) if !missing(tdate_y)
gen int    event_day_y   = date - tdate_y if !missing(tdate_y)
gen int    event_week_y  = floor(event_day_y/7) if !missing(event_day_y)


/**********************************************************************
  (1D) CONTRACT LIFECYCLE + SEASON RAMP + CROSS-SEASON EXPERIENCE

  contract_day: days since contract start (within-year definition)
  season_days_worked: cumulative *distinct days* worked within season (worker×year)
  cum_seasons: cumulative count of seasons observed up to current (0 in first observed season)
**********************************************************************/

* --- Contract lifecycle (within year) ---
cap drop contract_start contract_day contract_day2
bysort contract_id year: egen double contract_start = min(date)
gen int contract_day  = date - contract_start
gen int contract_day2 = contract_day^2

* --- Within-season ramp: count distinct work days for worker×year ---
cap drop first_day season_days_worked season_days_worked2
bysort worker_id year (date): gen byte first_day = (date != date[_n-1])
bysort worker_id year (date): gen int  season_days_worked = sum(first_day)
gen int season_days_worked2 = season_days_worked^2

* --- Cross-season experience: number of seasons observed so far ---
cap drop yrflag cum_seasons
bysort worker_id year: gen byte yrflag = (_n==1)
bysort worker_id (year): gen int cum_seasons = sum(yrflag) - 1
drop yrflag


/**********************************************************************
  (2) BASELINE ANALYSIS (RERUNNABLE CONTROL PANEL)
**********************************************************************/

/**********************************************************************
  2.0) CLEAN RERUN STATE (SAFE)
**********************************************************************/
cap estimates clear
cap label drop ew_lbl ew_lbl2

foreach v in ///
    sample_2A treated_BC_2A post_2A event_week_2A ///
    in_esamp_2A bc_obs_2A ctrl_obs_2A bcN_2A ctrlN_2A okweek_2A event_week_s_2A ///
    sample_2B treated_year_2B placebo_year_2B post_y_2B event_week_y_2B ///
    in_esamp_2B t_obs_2B p_obs_2B tN_2B pN_2B okweek_2B event_week_sy_2B {
    cap drop `v'
}

foreach m in X ABS SAMPLE_OK LIFECYCLE WPRE_2A WPOST_2A PRE_MIN_2B POST_MAX_2B EW_BASE_2A EW_BASE_2B {
    cap macro drop `m'
}


/**********************************************************************
  2.1) CONTROL PANEL (EDIT THESE)
**********************************************************************/

* Outcome variable (created in Section 1)
global Y "ln_incentive"      // ln_incentive / lprod / ln_piece / tot_hrs_day

* Controls toggles
global USE_WEATHER     0
global USE_EXPER       0
global USE_HOURS       0
global USE_DAYFLAGS    0
global USE_LIFECYCLE   0     // contract ramp + season ramp + ramp×cum_seasons

* Fixed effects choice:
*   "worker_contract" -> absorb(worker_id contract_id)
*   "wc_only"         -> absorb(wc_id)
*   "worker_only"     -> absorb(worker_id)
*   "contract_only"   -> absorb(contract_id)
global FE "worker_contract"

* Floors
global USE_FLOORS   0
global MIN_HRS      2
global MIN_PROD     50
global MIN_EARN     30
global USE_MIN_PROD 1

* AUTO window selection (TAU rule)
global TAU   0.25
global CORE0 0
global CORE1 1
global MAXW  30

global AUTO_WIN_2A 1
global AUTO_WIN_2B 1

* Manual windows used if AUTO is off or fails
global WPRE_MAN_2A   3
global WPOST_MAN_2A  8

global PRE_MIN_MAN_2B  -6
global POST_MAX_MAN_2B  10


/**********************************************************************
  2.2) BUILD CONTROLS ($X), LIFECYCLE ($LIFECYCLE), FE ($ABS), FLOORS ($SAMPLE_OK)
**********************************************************************/

* ---- Main controls ----
global X ""
if $USE_HOURS==1    global X "$X tot_hrs_day"
if $USE_DAYFLAGS==1 global X "$X multi_contract_day"
if $USE_EXPER==1    global X "$X c.experience##c.experience"
if $USE_WEATHER==1  global X "$X c.temp##c.temp c.wind_spd##c.wind_spd precip0 c.precip_pos##c.precip_pos"

* ---- Lifecycle/ramp controls (built safely) ----
global LIFECYCLE ""
if $USE_LIFECYCLE==1 {
    global LIFECYCLE "$LIFECYCLE c.contract_day##c.contract_day"
    global LIFECYCLE "$LIFECYCLE c.season_days_worked##c.season_days_worked"
    global LIFECYCLE "$LIFECYCLE c.season_days_worked##c.cum_seasons"
}

* ---- Fixed effects ----
global ABS ""
if "$FE"=="worker_contract" global ABS "worker_id contract_id"
if "$FE"=="wc_only"         global ABS "wc_id"
if "$FE"=="worker_only"     global ABS "worker_id"
if "$FE"=="contract_only"   global ABS "contract_id"

* ---- Floors / inclusion ----
global SAMPLE_OK "1==1"
if $USE_FLOORS==1 {
    global SAMPLE_OK "$SAMPLE_OK & tot_hrs_day >= $MIN_HRS"
    global SAMPLE_OK "$SAMPLE_OK & piece > 0 & prod > 0"
    global SAMPLE_OK "$SAMPLE_OK & (piece*prod) >= $MIN_EARN"
    if $USE_MIN_PROD==1 global SAMPLE_OK "$SAMPLE_OK & prod >= $MIN_PROD"
}

di as text "[DEBUG] Y:          $Y"
di as text "[DEBUG] X:          $X"
di as text "[DEBUG] LIFECYCLE:  $LIFECYCLE"
di as text "[DEBUG] ABS:        $ABS"
di as text "[DEBUG] FLOORS:     $SAMPLE_OK"


/**********************************************************************
  2A) BC vs AB/ON — treated year only
**********************************************************************/

* --------------------------
* 2A.1 SETTINGS
* --------------------------
global TY_2A 2018
global CONTROL_2A "AB"      // "AB" or "ON" (ABON allowed but auto-window off)

if "$CONTROL_2A"=="ABON" global AUTO_WIN_2A 0

* --------------------------
* 2A.2 DEFINE SAMPLE + EVENT TIME (ROBUST; NO $IF MACROS)
* --------------------------
cap drop sample_2A treated_BC_2A post_2A event_week_2A

gen byte sample_2A = (year==$TY_2A) & inlist(province_s,"BC","$CONTROL_2A")
if "$CONTROL_2A"=="ABON" replace sample_2A = (year==$TY_2A) & inlist(province_s,"BC","AB","ON")

gen byte treated_BC_2A = .
replace treated_BC_2A = (province_s=="BC") if sample_2A==1

local tdate_2A = mdy(6,1,$TY_2A)
gen byte post_2A       = (date >= `tdate_2A') if sample_2A==1
gen int  event_week_2A = floor((date - `tdate_2A')/7) if sample_2A==1

* sanity checks
tab province_s if sample_2A==1
tab treated_BC_2A if sample_2A==1
tab post_2A if sample_2A==1


/**********************************************************************
  2A.3 AUTO WINDOW SELECTION (TAU) OR MANUAL FALLBACK
  Uses actual regression sample: sample_2A & floors & !missing(Y)
**********************************************************************/
cap macro drop WPRE_2A WPOST_2A
global WPRE_2A  .
global WPOST_2A .

if $AUTO_WIN_2A==1 {

    local oksetup = 1
    preserve
        keep if sample_2A==1 ///
            & ($SAMPLE_OK) ///
            & !missing($Y) ///
            & inrange(event_week_2A, -$MAXW, $MAXW)

        keep if inlist(province_s,"BC","$CONTROL_2A")

        collapse (count) N=date, by(event_week_2A province_s)
        reshape wide N, i(event_week_2A) j(province_s) string

        capture confirm variable NBC
        if _rc local oksetup = 0
        local vctrl = "N$CONTROL_2A"
        capture confirm variable `vctrl'
        if _rc local oksetup = 0

        if `oksetup'==1 {
            quietly summarize NBC if inrange(event_week_2A,$CORE0,$CORE1)
            scalar coreBC = r(mean)
            quietly summarize `vctrl' if inrange(event_week_2A,$CORE0,$CORE1)
            scalar coreCT = r(mean)

            scalar thrBC = $TAU*coreBC
            scalar thrCT = $TAU*coreCT

            gen byte pass = (NBC>=thrBC & `vctrl'>=thrCT)
            sort event_week_2A

            gen byte cont_pre = .
            replace cont_pre = pass if event_week_2A==-1
            forvalues k=2/$MAXW {
                local wk = -`k'
                quietly replace cont_pre = pass & cont_pre[_n+1] if event_week_2A==`wk'
            }

            gen byte cont_post = .
            replace cont_post = pass if event_week_2A==0
            forvalues k=1/$MAXW {
                quietly replace cont_post = pass & cont_post[_n-1] if event_week_2A==`k'
            }

            quietly summarize event_week_2A if cont_pre==1, meanonly
            scalar PRE_sug = r(min)
            quietly summarize event_week_2A if cont_post==1, meanonly
            scalar POST_sug = r(max)

            quietly count if cont_pre==1
            scalar npre = r(N)
            quietly count if cont_post==1
            scalar npost = r(N)
        }
    restore

    if (`oksetup'==0) {
        di as error "[2A AUTO] FAIL: could not build counts (missing NBC or N$CONTROL_2A). Using manual."
        global WPRE_2A  $WPRE_MAN_2A
        global WPOST_2A $WPOST_MAN_2A
    }
    else if missing(PRE_sug) | npre==0 {
        di as error "[2A AUTO] FAIL: no PRE window meets TAU=" %6.4f $TAU " ending at -1. Using manual."
        global WPRE_2A  $WPRE_MAN_2A
        global WPOST_2A $WPOST_MAN_2A
    }
    else if missing(POST_sug) | npost==0 {
        di as error "[2A AUTO] FAIL: no POST window meets TAU=" %6.4f $TAU " starting at 0. Using manual."
        global WPRE_2A  $WPRE_MAN_2A
        global WPOST_2A $WPOST_MAN_2A
    }
    else {
        global WPRE_2A  = -PRE_sug
        global WPOST_2A =  POST_sug
    }

}
else {
    global WPRE_2A  $WPRE_MAN_2A
    global WPOST_2A $WPOST_MAN_2A
}

di as text "[2A] FINAL window: [-$WPRE_2A,$WPOST_2A]"


/**********************************************************************
  2A.4 DiD regression
**********************************************************************/
reghdfe $Y ///
    i.treated_BC_2A##i.post_2A ///
    $X ///
    $LIFECYCLE ///
    if sample_2A==1 ///
       & inrange(event_week_2A, -$WPRE_2A, $WPOST_2A) ///
       & ($SAMPLE_OK) ///
       & !missing($Y), ///
    absorb($ABS) ///
    vce(cluster contract_id)
estimates store did_2A


/**********************************************************************
  2A.5 Event study (support-restricted weeks) + plot
  - Restrict the factor variable's scope by requiring !missing(event_week_s_2A)
  - Plot only coefficients present in e(sample)
**********************************************************************/
cap drop in_esamp_2A bc_obs_2A ctrl_obs_2A bcN_2A ctrlN_2A okweek_2A event_week_s_2A
gen byte in_esamp_2A = (sample_2A==1) ///
    & inrange(event_week_2A, -$WPRE_2A, $WPOST_2A) ///
    & ($SAMPLE_OK) ///
    & !missing($Y) ///
    & inlist(province_s,"BC","$CONTROL_2A")

gen byte bc_obs_2A   = in_esamp_2A & (province_s=="BC")
gen byte ctrl_obs_2A = in_esamp_2A & (province_s=="$CONTROL_2A")

bysort event_week_2A: egen int bcN_2A   = total(bc_obs_2A)   if !missing(event_week_2A)
bysort event_week_2A: egen int ctrlN_2A = total(ctrl_obs_2A) if !missing(event_week_2A)

gen byte okweek_2A = in_esamp_2A & (bcN_2A>0) & (ctrlN_2A>0)
gen int  event_week_s_2A = event_week_2A + 100 if okweek_2A==1

quietly count if in_esamp_2A==1 & event_week_s_2A==99
if r(N)==0 {
    di as error "[2A] FAIL: week -1 not available after support restriction. Adjust windows."
    exit 198
}
global EW_BASE_2A 99

reghdfe $Y ///
    ib($EW_BASE_2A).event_week_s_2A##i.treated_BC_2A ///
    $X ///
    $LIFECYCLE ///
    if in_esamp_2A==1 & !missing(event_week_s_2A), ///
    absorb($ABS) ///
    vce(cluster worker_id)
estimates store es_2A

* coefplot: only coefficients in e(sample)
levelsof event_week_s_2A if e(sample), local(actual_levs_2A)

local keepcoef_2A ""
local clabel_2A   ""
foreach k of local actual_levs_2A {
    local keepcoef_2A `keepcoef_2A' `k'.event_week_s_2A#1.treated_BC_2A
    local val = `k' - 100
    local clabel_2A `clabel_2A' `k'.event_week_s_2A#1.treated_BC_2A = "`val'"
}

* vertical line between -1 and 0 based on the plotted sequence
preserve
    keep if e(sample)
    keep event_week_s_2A
    duplicates drop
    sort event_week_s_2A
    gen idx = _n
    quietly summarize idx if event_week_s_2A==100, meanonly
    local line_pos_2A = r(mean) - 0.5
restore

coefplot es_2A, ///
    keep(`keepcoef_2A') ///
    vertical ///
    omitted baselevels ///
    coeflabel(`clabel_2A') ///
    yline(0, lpattern(dash)) ///
    xline(`line_pos_2A', lwidth(vthin)) ///
    title("Event study: $Y ($TY_2A)", size(medium)) ///
    subtitle("BC vs $CONTROL_2A | Window [-$WPRE_2A,$WPOST_2A] | TAU=$TAU", size(small)) ///
    ytitle("Effect relative to week -1") ///
    xtitle("Weeks relative to June 1 cutoff") ///
    graphregion(color(white))

graph export "event_2A_${Y}_${TY_2A}_${CONTROL_2A}.png", replace



/**********************************************************************
  2B) BC placebo-years — within BC
**********************************************************************/

* --------------------------
* 2B.1 SETTINGS
* --------------------------
global TY_2B 2018
global BCPLACEBO_2B "2016 2017"

* --------------------------
* 2B.2 DEFINE SAMPLE + GROUPS (NO $IF MACROS)
* --------------------------
cap drop treated_year_2B placebo_year_2B sample_2B post_y_2B event_week_y_2B

gen byte treated_year_2B = (year==$TY_2B)

gen byte placebo_year_2B = 0
foreach yy of numlist $BCPLACEBO_2B {
    replace placebo_year_2B = 1 if year==`yy'
}

gen byte sample_2B = (province_s=="BC") & (treated_year_2B==1 | placebo_year_2B==1)

* year-specific event time already exists as post_y / event_week_y
gen byte post_y_2B = post_y if sample_2B==1
gen int  event_week_y_2B = event_week_y if sample_2B==1


/**********************************************************************
  2B.3 AUTO WINDOW SELECTION (TAU) OR MANUAL FALLBACK
  Uses actual regression sample: sample_2B & floors & !missing(Y)
**********************************************************************/
cap macro drop PRE_MIN_2B POST_MAX_2B
global PRE_MIN_2B  .
global POST_MAX_2B .

if $AUTO_WIN_2B==1 {

    local oksetup = 1
    preserve
        keep if sample_2B==1 ///
            & ($SAMPLE_OK) ///
            & !missing($Y) ///
            & inrange(event_week_y_2B, -$MAXW, $MAXW)

        collapse (count) N=date, by(event_week_y_2B treated_year_2B)
        reshape wide N, i(event_week_y_2B) j(treated_year_2B)

        capture confirm variable N0
        if _rc local oksetup = 0
        capture confirm variable N1
        if _rc local oksetup = 0

        if `oksetup'==1 {
            quietly summarize N1 if inrange(event_week_y_2B,$CORE0,$CORE1)
            scalar coreT = r(mean)
            quietly summarize N0 if inrange(event_week_y_2B,$CORE0,$CORE1)
            scalar coreP = r(mean)

            scalar thrT = $TAU*coreT
            scalar thrP = $TAU*coreP

            gen byte pass = (N1>=thrT & N0>=thrP)
            sort event_week_y_2B

            gen byte cont_pre = .
            replace cont_pre = pass if event_week_y_2B==-1
            forvalues k=2/$MAXW {
                local wk = -`k'
                quietly replace cont_pre = pass & cont_pre[_n+1] if event_week_y_2B==`wk'
            }

            gen byte cont_post = .
            replace cont_post = pass if event_week_y_2B==0
            forvalues k=1/$MAXW {
                quietly replace cont_post = pass & cont_post[_n-1] if event_week_y_2B==`k'
            }

            quietly summarize event_week_y_2B if cont_pre==1, meanonly
            scalar PRE_sug = r(min)
            quietly summarize event_week_y_2B if cont_post==1, meanonly
            scalar POST_sug = r(max)

            quietly count if cont_pre==1
            scalar npre = r(N)
            quietly count if cont_post==1
            scalar npost = r(N)
        }
    restore

    if (`oksetup'==0) {
        di as error "[2B AUTO] FAIL: could not build treated/placebo week counts. Using manual."
        global PRE_MIN_2B  $PRE_MIN_MAN_2B
        global POST_MAX_2B $POST_MAX_MAN_2B
    }
    else if missing(PRE_sug) | npre==0 {
        di as error "[2B AUTO] FAIL: no PRE window meets TAU=" %6.4f $TAU " ending at -1. Using manual."
        global PRE_MIN_2B  $PRE_MIN_MAN_2B
        global POST_MAX_2B $POST_MAX_MAN_2B
    }
    else if missing(POST_sug) | npost==0 {
        di as error "[2B AUTO] FAIL: no POST window meets TAU=" %6.4f $TAU " starting at 0. Using manual."
        global PRE_MIN_2B  $PRE_MIN_MAN_2B
        global POST_MAX_2B $POST_MAX_MAN_2B
    }
    else {
        global PRE_MIN_2B  = PRE_sug
        global POST_MAX_2B = POST_sug
    }

}
else {
    global PRE_MIN_2B  $PRE_MIN_MAN_2B
    global POST_MAX_2B $POST_MAX_MAN_2B
}

di as text "[2B] FINAL window: [$PRE_MIN_2B,$POST_MAX_2B]"


/**********************************************************************
  2B.4 DiD regression (within BC)
**********************************************************************/
reghdfe $Y ///
    i.treated_year_2B##i.post_y_2B ///
    $X ///
    $LIFECYCLE ///
    if sample_2B==1 ///
       & inrange(event_week_y_2B, $PRE_MIN_2B, $POST_MAX_2B) ///
       & ($SAMPLE_OK) ///
       & !missing($Y), ///
    absorb($ABS) ///
    vce(cluster contract_id)
estimates store did_2B


/**********************************************************************
  2B.5 Event study (support-restricted weeks) + plot
**********************************************************************/
cap drop in_esamp_2B t_obs_2B p_obs_2B tN_2B pN_2B okweek_2B event_week_sy_2B
gen byte in_esamp_2B = (sample_2B==1) ///
    & inrange(event_week_y_2B, $PRE_MIN_2B, $POST_MAX_2B) ///
    & ($SAMPLE_OK) ///
    & !missing($Y)

gen byte t_obs_2B = in_esamp_2B & (treated_year_2B==1)
gen byte p_obs_2B = in_esamp_2B & (treated_year_2B==0)

bysort event_week_y_2B: egen int tN_2B = total(t_obs_2B) if !missing(event_week_y_2B)
bysort event_week_y_2B: egen int pN_2B = total(p_obs_2B) if !missing(event_week_y_2B)

gen byte okweek_2B = in_esamp_2B & (tN_2B>0) & (pN_2B>0)
gen int  event_week_sy_2B = event_week_y_2B + 100 if okweek_2B==1

quietly count if in_esamp_2B==1 & event_week_sy_2B==99
if r(N)==0 {
    di as error "[2B] FAIL: week -1 not available after support restriction. Adjust windows."
    exit 198
}
global EW_BASE_2B 99

reghdfe $Y ///
    ib($EW_BASE_2B).event_week_sy_2B##i.treated_year_2B ///
    $X ///
    $LIFECYCLE ///
    if in_esamp_2B==1 & !missing(event_week_sy_2B), ///
    absorb($ABS) ///
    vce(cluster contract_id)
estimates store es_2B

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
    quietly summarize idx if event_week_sy_2B==100, meanonly
    local line_pos_2B = r(mean) - 0.5
restore

coefplot es_2B, ///
    keep(`keepcoef_2B') ///
    vertical ///
    omitted baselevels ///
    coeflabel(`clabel_2B') ///
    yline(0, lpattern(dash)) ///
    xline(`line_pos_2B', lwidth(vthin)) ///
    title("Event study: $Y (BC only)", size(medium)) ///
    subtitle("Treated $TY_2B vs placebo ($BCPLACEBO_2B) | Window [$PRE_MIN_2B,$POST_MAX_2B] | TAU=$TAU", size(small)) ///
    ytitle("Effect relative to week -1") ///
    xtitle("Weeks relative to June 1 cutoff") ///
    graphregion(color(white))

graph export "event_2B_${Y}_${TY_2B}.png", replace

/**********************************************************************
  (3) ADDITIONAL EMPIRICAL ANALYSIS  (SELF-CONTAINED)

  PURPOSE:
    - Ability construction (BC-only, pre-June 1, predetermined)
    - Heterogeneous MW effects (primarily BC placebo design)
    - Optional: cross-province heterogeneity (BC vs AB/ON) using BC ability
    - Selection tests (who shows up post?)
    - Sorting tests (do high-ability workers shift into higher piece-rate contracts?)

  REQUIREMENTS:
    - Section (1) DATA PREP has already been run in this session so that:
        date, year, province, worker_id, contract_id, wc_id exist
        ln_incentive, lprod exist
        tot_hrs_day, multi_contract_day exist
        temp, wind_spd, precip0, precip_pos exist
        experience, experience2 exist (or experience only; we use quadratic)
**********************************************************************/

***********************************************************************
* 3.0) FORCE CLEAR (SAFE RERUN OF SECTION 3 ONLY)
***********************************************************************
cap estimates clear

foreach v in treated_year placebo_year bc_placebo_sample ///
            tdate_y post_y event_week_y ///
            ln_prod abil_sample resid_pre ability_raw ability_z ///
            treated_BC_3X post_3X event_week_3X ///
            any_pre any_post spans ///
            ln_piece {
    cap drop `v'
}
foreach m in TY_BC BCPLACEBO PRE_ABIL PRE_MAIN POST_MAIN ///
            CONTROL_3X TY_3X WPRE_3X WPOST_3X ///
            USE_FLOORS_3 MIN_HRS_3 MIN_EARN_3 MIN_PROD_3 USE_MIN_PROD_3 ///
            X_3 X_ABIL SAMPLE_OK_3 {
    cap macro drop `m'
}


***********************************************************************
* 3.1) USER SETTINGS (EDIT THESE)
***********************************************************************

* --- Treated MW year in BC for placebo design ---
global TY_BC 2019

* --- Placebo years in BC (no MW change) ---
global BCPLACEBO "2016 2017"

* --- Ability window (how far back pre-June 1 to use) ---
* (This is ONLY for estimating ability; it can be longer than your main window.)
global PRE_ABIL -12

* --- Main heterogeneity window for BC placebo design ---
* (This is the window used in the heterogeneity regression.)
global PRE_MAIN  -8
global POST_MAIN  10

* --- Floors (keep consistent with your baseline; edit as needed) ---
global USE_FLOORS_3 1
global MIN_HRS_3    2
global MIN_EARN_3   30
global USE_MIN_PROD_3 1
global MIN_PROD_3   50     // prod assumed = trees/hour

* --- Optional cross-province heterogeneity design (BC vs AB/ON in treated year) ---
* If you don't want this, set RUN_3X = 0.
global RUN_3X 0
global TY_3X 2019
global CONTROL_3X "AB"
global WPRE_3X 3
global WPOST_3X 8


***********************************************************************
* 3.2) BUILD CONTROLS + SAMPLE FLOORS (USED THROUGHOUT SECTION 3)
***********************************************************************

* Controls used in heterogeneity/selection/sorting regressions
global X_3 "tot_hrs_day multi_contract_day"
global X_3 "$X_3 c.experience##c.experience"
global X_3 "$X_3 c.temp##c.temp c.wind_spd##c.wind_spd"
global X_3 "$X_3 precip0 c.precip_pos##c.precip_pos"

* Controls used in ABILITY regression (keep aligned)
global X_ABIL "$X_3"

* Floors (as a single string condition)
global SAMPLE_OK_3 "1==1"
if $USE_FLOORS_3==1 {
    global SAMPLE_OK_3 "$SAMPLE_OK_3 & tot_hrs_day >= $MIN_HRS_3"
    global SAMPLE_OK_3 "$SAMPLE_OK_3 & piece > 0 & prod > 0"
    global SAMPLE_OK_3 "$SAMPLE_OK_3 & (piece*prod) >= $MIN_EARN_3"
    if $USE_MIN_PROD_3==1 global SAMPLE_OK_3 "$SAMPLE_OK_3 & prod >= $MIN_PROD_3"
}

di as text "[Section 3] X_3 = $X_3"
di as text "[Section 3] SAMPLE_OK_3 = $SAMPLE_OK_3"


/**********************************************************************
  3.3) DEFINE BC PLACEBO DESIGN OBJECTS (SELF-CONTAINED)
**********************************************************************/

* Treated vs placebo indicators
gen byte treated_year = (year == $TY_BC)

gen byte placebo_year = 0
foreach yy of numlist $BCPLACEBO {
    replace placebo_year = 1 if year == `yy'
}

* BC placebo sample: BC only, treated year + placebo years
gen byte bc_placebo_sample = (province=="BC") & (treated_year==1 | placebo_year==1)

* Year-specific June 1 cutoff and event time within BC placebo sample
gen double tdate_y = mdy(6,1,year)
format tdate_y %td

gen byte post_y = (date >= tdate_y) if bc_placebo_sample
gen int  event_week_y = floor((date - tdate_y)/7) if bc_placebo_sample

* Quick checks
di as text "[3.3] BC placebo sample counts"
tab year if bc_placebo_sample
tab post_y if bc_placebo_sample
summ event_week_y if bc_placebo_sample


/**********************************************************************
  3.4) CONSTRUCT PREDETERMINED ABILITY (BC-ONLY, PRE-JUNE 1)
  Approach:
    - Use ln(prod) as outcome (prod assumed trees/hour)
    - Residualize ln(prod) on controls + contract FE (difficulty)
    - Ability(worker) = mean residual over pre-period days
**********************************************************************/

* Outcome for ability
gen double ln_prod = lprod   // lprod created in Section 1 as ln(prod) if prod>0

* Ability estimation sample: BC placebo sample, pre only, within PRE_ABIL..-1
gen byte abil_sample = bc_placebo_sample ///
    & inrange(event_week_y, $PRE_ABIL, -1) ///
    & (post_y==0) ///
    & !missing(ln_prod) ///
    & ($SAMPLE_OK_3)

di as text "[3.4] Ability sample counts"
count if abil_sample==1
tab year if abil_sample==1
summ event_week_y if abil_sample==1

* Residualize on controls + contract FE using ONLY ability sample
reghdfe ln_prod $X_ABIL ///
    if abil_sample==1, ///
    absorb(contract_id) ///
    resid

* Now residuals are available directly as _reghdfe_resid
cap drop resid_pre
gen double resid_pre = _reghdfe_resid if abil_sample==1


* Worker-level ability = mean residual across their pre-period days
bys worker_id: egen double ability_raw = mean(resid_pre)

* Standardize (over workers with ability defined)
summ ability_raw if !missing(ability_raw)
gen double ability_z = (ability_raw - r(mean)) / r(sd) if !missing(ability_raw)

di as text "[3.4] Ability coverage checks"
count if bc_placebo_sample & year==$TY_BC & !missing(ability_z)
count if bc_placebo_sample & !missing(ability_z)
tab treated_year if bc_placebo_sample & !missing(ability_z)


/**********************************************************************
  3.5) HETEROGENEOUS MW EFFECTS (BC PLACEBO DESIGN)
  Model (within BC):
    Y = treated_year + post_y + treated_year×post_y
        + ability + ability×post + ability×treated_year + ability×treated_year×post
        + controls + FE
  Notes:
    - Use $Y_HET for outcome; if you want to reuse baseline $Y, set it here.
**********************************************************************/

* Choose outcome for heterogeneity (edit)
capture confirm global Y_HET
if _rc global Y_HET "ln_incentive"

* Explicit interactions (transparent)
cap drop ty_post ability_post ability_ty ability_ty_post
gen byte   ty_post        = (treated_year==1 & post_y==1) if bc_placebo_sample
gen double ability_post   = ability_z * post_y           if bc_placebo_sample & !missing(ability_z)
gen double ability_ty     = ability_z * treated_year     if bc_placebo_sample & !missing(ability_z)
gen double ability_ty_post= ability_z * ty_post          if bc_placebo_sample & !missing(ability_z)

* Heterogeneity regression
reghdfe $Y_HET ///
    treated_year post_y ty_post ///
    ability_z ability_post ability_ty ability_ty_post ///
    $X_3 ///
    if bc_placebo_sample ///
       & inrange(event_week_y, $PRE_MAIN, $POST_MAIN) ///
       & !missing($Y_HET) ///
       & ($SAMPLE_OK_3) ///
       & !missing(ability_z), ///
    absorb(worker_id contract_id) ///
    vce(cluster contract_id)

estimates store het_2B

* Diagnostics: ensure both pre/post are in the estimation sample
tab post_y if e(sample)
tab treated_year if e(sample)


/**********************************************************************
  3.6) OPTIONAL: CROSS-PROVINCE HETEROGENEITY (BC vs AB/ON, TREATED YEAR)
  Important:
    - ability_z is only defined for BC workers (constructed from BC-only pre data).
    - This regression uses ability heterogeneity ONLY within BC via interactions.
    - If you set RUN_3X=0, this block is skipped.
**********************************************************************/

if $RUN_3X==1 {

    * Define cross-province sample in treated year
    gen byte treated_BC_3X = (province=="BC") if year==$TY_3X
    gen byte post_3X       = (date >= td(01jun$TY_3X)) if year==$TY_3X
    gen int  event_week_3X = floor((date - td(01jun$TY_3X))/7) if year==$TY_3X

    cap drop sample_3X
    gen byte sample_3X = (year==$TY_3X) ///
        & (province=="BC" | province=="$CONTROL_3X") ///
        & inrange(event_week_3X, -$WPRE_3X, $WPOST_3X) ///
        & ($SAMPLE_OK_3) ///
        & !missing($Y_HET)

    * Ability interactions only apply to BC (treated_BC_3X==1)
    cap drop treat_post_3X ability_post_3X ability_treat_3X ability_treat_post_3X
    gen byte   treat_post_3X         = (treated_BC_3X==1 & post_3X==1) if sample_3X
    gen double ability_post_3X       = ability_z * post_3X            if sample_3X & treated_BC_3X==1 & !missing(ability_z)
    gen double ability_treat_3X      = ability_z * treated_BC_3X      if sample_3X & treated_BC_3X==1 & !missing(ability_z)
    gen double ability_treat_post_3X = ability_z * treat_post_3X      if sample_3X & treated_BC_3X==1 & !missing(ability_z)

    reghdfe $Y_HET ///
        treated_BC_3X post_3X treat_post_3X ///
        ability_z ability_post_3X ability_treat_3X ability_treat_post_3X ///
        $X_3 ///
        if sample_3X, ///
        absorb(worker_id contract_id) ///
        vce(cluster contract_id)

    estimates store het_2A

    tab post_3X if e(sample)
    tab treated_BC_3X if e(sample)
}


/**********************************************************************
  3.7) SELECTION TESTS
  Question:
    - Are low-ability workers less likely to be observed post-June 1
      in treated year relative to placebo years (within BC)?

  Construct worker-year panel among BC placebo sample:
    any_pre: worker worked pre in that year
    any_post: worker worked post in that year
  Then regress any_post on ability × treated_year among those with any_pre==1.
**********************************************************************/

preserve
    keep if bc_placebo_sample==1
    keep worker_id year post_y ability_z treated_year placebo_year

    bysort worker_id year: egen byte any_pre  = max(post_y==0)
    bysort worker_id year: egen byte any_post = max(post_y==1)

    bysort worker_id year: keep if _n==1
    keep if any_pre==1

    * LPM selection test (cluster worker)
    reg any_post c.ability_z##i.treated_year, vce(cluster worker_id)

restore


/**********************************************************************
  3.8) SORTING / MATCHING TESTS (DIAGNOSTIC)
  Question:
    - Does the relationship between ability and piece rate change post
      in treated year relative to placebo years?

  Outcome:
    ln_piece = ln(piece)
  Model:
    ln_piece = treated_year + post_y + treated_year×post_y
               + ability + ability×post + ability×treated_year + ability×treated_year×post
               + controls + FE

  Note:
    - This is NOT a causal wage equation; it's a matching/sorting diagnostic.
**********************************************************************/

gen double ln_piece = ln(piece) if piece>0

cap drop ty_post_p ability_post_p ability_ty_p ability_ty_post_p
gen byte   ty_post_p         = (treated_year==1 & post_y==1) if bc_placebo_sample
gen double ability_post_p    = ability_z * post_y           if bc_placebo_sample & !missing(ability_z)
gen double ability_ty_p      = ability_z * treated_year     if bc_placebo_sample & !missing(ability_z)
gen double ability_ty_post_p = ability_z * ty_post_p        if bc_placebo_sample & !missing(ability_z)

reghdfe ln_piece ///
    treated_year post_y ty_post_p ///
    ability_z ability_post_p ability_ty_p ability_ty_post_p ///
    $X_3 ///
    if bc_placebo_sample ///
       & inrange(event_week_y, $PRE_MAIN, $POST_MAIN) ///
       & !missing(ln_piece) ///
       & ($SAMPLE_OK_3) ///
       & !missing(ability_z), ///
    absorb(worker_id contract_id) ///
    vce(cluster contract_id)

estimates store sort_2B

tab post_y if e(sample)
tab treated_year if e(sample)


/**********************************************************************
  3.9) OPTIONAL: SIMPLE CONTRACT SORTING BY ABILITY (DESCRIPTIVE)
  Useful quick diagnostic:
    - Compare mean piece rate / mean productivity by ability bins pre vs post.
**********************************************************************/

cap drop abil_bin
xtile abil_bin = ability_z if bc_placebo_sample & !missing(ability_z), nq(4)

table abil_bin post_y if bc_placebo_sample & inrange(event_week_y,$PRE_MAIN,$POST_MAIN), ///
    statistic(mean piece) statistic(mean prod) statistic(mean tot_hrs_day) nformat(%9.3f)
	
	
	
	
	
	

	
	/*
	
/**********************************************************************
  (3) ADDITIONAL ANALYSIS: BC-ONLY ABILITY + HETEROGENEITY USING BC PLACEBOS

  GOAL:
    - Construct worker ability using BC-only PRE data (pre-determined)
    - Estimate heterogeneous MW effects using BC treated-year vs BC placebo years
      (diff-in-discontinuities around June 1)

  THIS SECTION IS SELF-CONTAINED (does NOT require Section 2).

  REQUIREMENTS (from Section 1 data prep):
    - date (Stata daily), year, province, contract_id, worker_id
    - outcomes: lprod, ln_incentive, etc.
    - controls: experience, temp, wind_spd, precip0, precip_pos, tot_hrs_day, multi_contract_day
**********************************************************************/

***********************************************************************
* 3.0) USER CHOICES
***********************************************************************

* Treated MW year in BC
global TY_BC 2019

* BC placebo years with no MW change (space separated)
global BCPLACEBO "2016 2017"

* Window for MAIN heterogeneity estimation (BC has lots of mass)
global PRE_MAIN   -8
global POST_MAIN   10

* Window used ONLY for ability estimation (can be longer than PRE_MAIN)
* (This is BC-only so you can usually go longer without support problems.)
global PRE_ABILITY -12

* Outcome for heterogeneity regressions (must exist)
* Options you use: ln_incentive, lprod, lpiece, tot_hrs, etc.
global Y_HET "ln_incentive"

* Floors (reuse your Section 2 logic if you want; here is a simple version)
global USE_FLOORS 1
global MIN_HRS   2
global MIN_EARN  30
global MIN_PROD  50   // if prod is trees/hour; adjust if needed

/**********************************************************************
  3.1) PRE-DETERMINED WORKER ABILITY (BC-ONLY, PRE-JUNE 1)
  Goal:
    - Construct ability using ONLY BC observations in:
        treated year (TY_BC) + placebo years
      and ONLY pre-cutoff weeks (event_week_y < 0).
    - Ability is a worker-level scalar, then attached to ALL rows for that worker.

  Method (transparent, robust):
    1) In pre-period BC-only sample, residualize ln_prod on:
         - controls (same style as main)
         - contract fixed effects (absorbs difficulty)
    2) Ability_raw(worker) = mean residual across that worker's pre-period days
    3) Standardize to ability_z

  REQUIREMENTS (from Section 1 + Section 2B):
    - worker_id contract_id
    - lprod (or ln_prod) outcome
    - bc_placebo_sample, post_y, event_week_y
**********************************************************************/

* --- choose the productivity outcome used for ability ---
* If you want ability based on trees/hour productivity:
cap drop ln_prod
gen double ln_prod = lprod   // lprod already = ln(prod) if prod>0

* --- define ability window (you can change this) ---
* use as wide a pre window as you want for ability; does NOT need to match main window
cap macro drop PRE_ABIL
global PRE_ABIL -12   // e.g., use up to 12 weeks pre; adjust as desired

* --- pre-period ability sample (BC placebo design, pre only) ---
cap drop abil_sample
gen byte abil_sample = (bc_placebo_sample==1) ///
    & (event_week_y < 0) ///
    & inrange(event_week_y, $PRE_ABIL, -1) ///
    & !missing(ln_prod)

di as text "[ABILITY] Pre sample counts:"
count if abil_sample==1
tab year if abil_sample==1

* --- controls for ability regression (build it explicitly; no multi-line globals) ---
cap macro drop X_ABIL
global X_ABIL "tot_hrs_day multi_contract_day"
global X_ABIL "$X_ABIL c.experience##c.experience"
global X_ABIL "$X_ABIL c.temp##c.temp c.wind_spd##c.wind_spd"
global X_ABIL "$X_ABIL precip0 c.precip_pos##c.precip_pos"

* --- residualize ln_prod on controls + contract FE using ONLY pre data ---
cap drop resid_pre
reghdfe ln_prod $X_ABIL if abil_sample==1, ///
    absorb(contract_id) ///
    vce(cluster contract_id)

predict double resid_pre if e(sample), resid

* --- worker-level ability = mean residual over pre-period days ---
cap drop ability_raw ability_z
bys worker_id: egen double ability_raw = mean(resid_pre)

* Propagate to all observations for that worker (ability_raw already constant within worker
* once computed, but will be missing for workers who never appear in abil_sample)
* No further action needed—egen computed mean over missing/nonmissing appropriately.

* --- standardize over workers with ability defined ---
quietly summ ability_raw if !missing(ability_raw)
gen double ability_z = (ability_raw - r(mean)) / r(sd) if !missing(ability_raw)

* --- diagnostics: should now be non-missing for many BC treated-year workers ---
di as text "[ABILITY] Coverage checks:"
count if province=="BC" & year==$TY_BC & !missing(ability_z)
tab year if province=="BC" & !missing(ability_z)

* (optional) sanity: ability is constant within worker
bys worker_id: egen byte _nabil = nvals(ability_z)
summ _nabil if !missing(ability_z)
drop _nabil

***********************************************************************
* 3.2) DEFINE A SINGLE "OK SAMPLE" FILTER (NO DROPPING OBS)
***********************************************************************

cap macro drop SAMPLE_OK_3
global SAMPLE_OK_3 "1==1"

if $USE_FLOORS==1 {
    * require hours/day
    global SAMPLE_OK_3 "$SAMPLE_OK_3 & tot_hrs_day >= $MIN_HRS"

    * require positive piece & prod for ln_incentive
    * (safe for other outcomes too)
    global SAMPLE_OK_3 "$SAMPLE_OK_3 & piece>0 & prod>0"

    * incentive earnings floor (piece*prod)
    global SAMPLE_OK_3 "$SAMPLE_OK_3 & (piece*prod) >= $MIN_EARN"

    * productivity floor (if prod = trees/hour)
    global SAMPLE_OK_3 "$SAMPLE_OK_3 & prod >= $MIN_PROD"
}

display as text "[Section 3] SAMPLE_OK_3 = $SAMPLE_OK_3"

***********************************************************************
* 3.3) PRE-DETERMINED ABILITY (BC ONLY, PRE PERIOD ONLY)
*
* We avoid fragile reghdfe savefe naming by using:
*   ability_raw = worker mean of residual productivity in pre sample,
*   where residuals come from regressing lprod on controls + contract FE.
*
* Interpretation:
*   ability_raw = persistent worker productivity component net of
*                 contract difficulty and observables, measured pre-policy.
***********************************************************************

cap drop abil_sample resid_pre ability_raw ability_z ability_hi

* Ability estimation sample: BC only, treated+placebos, PRE weeks only
gen byte abil_sample = bc_placebo_sample ///
    & inrange(event_week_y, $PRE_ABILITY, -1) ///
    & !missing(lprod) ///
    & ($SAMPLE_OK_3)

* Controls used in ability regression (keep consistent with main work)
* (You can add/remove stuff, but keep contract FE in ability estimation.)
global X_ABIL "tot_hrs_day multi_contract_day"
global X_ABIL "$X_ABIL c.experience##c.experience"
global X_ABIL "$X_ABIL c.temp##c.temp c.wind_spd##c.wind_spd"
global X_ABIL "$X_ABIL precip0 c.precip_pos##c.precip_pos"

* Ability regression: PRE only, BC only
reg lprod $X_ABIL i.contract_id if abil_sample

* Residual productivity net of controls and contract
predict double resid_pre if abil_sample, resid

* Worker ability = worker mean residual in pre sample
bysort worker_id: egen double ability_raw = mean(resid_pre) if abil_sample

* Spread to all observations for that worker (within BC sample):
bysort worker_id: egen double ability_raw_all = max(ability_raw)
drop ability_raw
rename ability_raw_all ability_raw

* Standardize ability using the ABILITY estimation sample (pre only)
summ ability_raw if abil_sample & !missing(ability_raw), meanonly
gen double ability_z = (ability_raw - r(mean)) / r(sd) if !missing(ability_raw)

* Median split (optional; useful for plots / group event studies)
xtile ability_hi = ability_raw if abil_sample & !missing(ability_raw), nq(2)
replace ability_hi = ability_hi - 1   // 0=low, 1=high

display as text "========== Ability coverage =========="
count if bc_placebo_sample & !missing(ability_z)
tab treated_year if bc_placebo_sample & !missing(ability_z)
summ ability_raw ability_z if bc_placebo_sample & !missing(ability_z)
display as text "======================================"

***********************************************************************
* 3.4) HETEROGENEITY: BC PLACEBO DIFF-IN-DISC (CONTINUOUS ABILITY)
*
* Model (BC only):
*   Y = treated_year + post_y + treated_year×post_y
*       + ability + ability×post_y + ability×treated_year
*       + ability×treated_year×post_y
*       + controls + FE
*
* Key coefficient:
*   ability×treated_year×post_y  (heterogeneity of the policy discontinuity)
***********************************************************************
/**********************************************************************
  3.4 DIAGNOSTICS — robust sample flag (no macro quoting issues)
**********************************************************************/

cap drop sample34
gen byte sample34 = 0

* Put your intended 3.4 logic here directly (edit provinces/years/vars as needed)
replace sample34 = 1 if ///
    (year == $TY_2A) ///
    & (province == "BC" | province == "$CONTROL_2A") ///
    & inrange(event_week_2A, -$WPRE_2A, $WPOST_2A)

* Now diagnostics never fail:
count if sample34==1
count if sample34==1 & !missing($Y)
count if sample34==1 & !missing($Y) & ($SAMPLE_OK)

* If 3.4 needs ability:
count if sample34==1 & !missing($Y) & ($SAMPLE_OK) & !missing(ability_z)

* Pre/post support:
tab post_2A if sample34==1 & !missing($Y) & ($SAMPLE_OK) & !missing(ability_z)
tab treated_BC_2A if sample34==1 & !missing($Y) & ($SAMPLE_OK) & !missing(ability_z)
* Main regression window for the heterogeneity estimand
cap noi reghdfe $Y_HET ///
    i.treated_year##i.post_y##c.ability_z ///
    tot_hrs_day multi_contract_day ///
    c.experience##c.experience ///
    c.temp##c.temp c.wind_spd##c.wind_spd ///
    precip0 c.precip_pos##c.precip_pos ///
    if bc_placebo_sample ///
       & inrange(event_week_y, $PRE_MAIN, $POST_MAIN) ///
       & ($SAMPLE_OK_3) ///
       & !missing($Y_HET) ///
       & !missing(ability_z), ///
    absorb(worker_id contract_id) ///
    vce(cluster contract_id)

estimates store het_bc_placebo

* Quick support check in estimating sample
tab post_y if e(sample)
tab treated_year if e(sample)

***********************************************************************
* 3.5) EVENT STUDY (BC PLACEBO DESIGN) — BY ABILITY GROUP (LOW vs HIGH)
*
* This is often the cleanest way to "show no pretrends" separately by ability:
* Run the same event study within BC, separately for ability_hi=0 and ability_hi=1,
* then plot the treated-year coefficients.
***********************************************************************

capture which coefplot
if _rc ssc install coefplot, replace

* Shifted event-week restricted to the MAIN window for plotting
cap drop event_week_sy_main
gen int event_week_sy_main = event_week_y + 100 if bc_placebo_sample ///
    & inrange(event_week_y, $PRE_MAIN, $POST_MAIN)

* Base week = -1 => shifted 99
global EW_BASE_BC 99

* Label only the weeks we actually use
cap label drop ew_lbl
levelsof event_week_sy_main if bc_placebo_sample & !missing(event_week_sy_main), local(levs_bc)
foreach k of local levs_bc {
    local v = `k' - 100
    label define ew_lbl `k' "`v'", add
}
label values event_week_sy_main ew_lbl

* LOW ability group
reghdfe $Y_HET ///
    ib($EW_BASE_BC).event_week_sy_main##i.treated_year ///
    tot_hrs_day multi_contract_day ///
    c.experience##c.experience ///
    c.temp##c.temp c.wind_spd##c.wind_spd ///
    precip0 c.precip_pos##c.precip_pos ///
    if bc_placebo_sample ///
       & ability_hi==0 ///
       & !missing(event_week_sy_main) ///
       & ($SAMPLE_OK_3) ///
       & !missing($Y_HET), ///
    absorb(worker_id contract_id) ///
    vce(cluster contract_id)
estimates store es_low

* HIGH ability group
reghdfe $Y_HET ///
    ib($EW_BASE_BC).event_week_sy_main##i.treated_year ///
    tot_hrs_day multi_contract_day ///
    c.experience##c.experience ///
    c.temp##c.temp c.wind_spd##c.wind_spd ///
    precip0 c.precip_pos##c.precip_pos ///
    if bc_placebo_sample ///
       & ability_hi==1 ///
       & !missing(event_week_sy_main) ///
       & ($SAMPLE_OK_3) ///
       & !missing($Y_HET), ///
    absorb(worker_id contract_id) ///
    vce(cluster contract_id)
estimates store es_high

* Build coef labels for treated-year coefficients only
local clabel ""
foreach k of local levs_bc {
    local v = `k' - 100
    local clabel `clabel' `k'.event_week_sy_main#1.treated_year = "`v'"
}

* Compute xline location between -1 and 0 based on actual plotted weeks
preserve
    keep if bc_placebo_sample & !missing(event_week_sy_main)
    keep event_week_sy_main
    duplicates drop
    sort event_week_sy_main
    gen idx = _n
    quietly summarize idx if event_week_sy_main==100, meanonly
    local line_pos = r(mean) - 0.5
restore

coefplot ///
    (es_low,  keep(*.event_week_sy_main#1.treated_year) label("Low ability")) ///
    (es_high, keep(*.event_week_sy_main#1.treated_year) label("High ability")) ///
    , vertical ///
      coeflabel(`clabel') ///
      yline(0, lpattern(dash)) ///
      xline(`line_pos', lwidth(vthin)) ///
      title("BC placebo event study by ability: $Y_HET", size(medium)) ///
      subtitle("Treated year $TY_BC vs placebo ($BCPLACEBO) | Window [$PRE_MAIN,$POST_MAIN]", size(small)) ///
      ytitle("Effect relative to week -1") ///
      xtitle("Weeks relative to June 1 cutoff") ///
      graphregion(color(white))

graph export "es_bc_placebo_by_ability_${Y_HET}_${TY_BC}.png", replace













/*
/**********************************************************************
  (3) ADDITIONAL EMPIRICAL ANALYSIS
  Purpose:
    3.1) Estimate PRE-treatment worker ability (AKM-style) using only pre–June 1 data
         and worker + contract FE. Then attach ability to all observations.
    3.2) Heterogeneous MW effects:
         (A) Cross-province (BC vs AB/ON) in treated year (2A design)
         (B) Within-BC placebo design (treated year vs placebo years) (2B design)
         (C) Optional triple-diff (BC vs AB) × (treated vs placebo) × (post)
    3.3) Selection diagnostics (do low-ability workers disappear post?)
    3.4) Sorting diagnostics (does ability predict piece-rate/contract allocation differently post?)

  Notes:
    - This section is written to be runnable AFTER Section 1 and Section 2.
    - It does NOT drop observations permanently; any temporary reshaping is preserve/restore.
**********************************************************************/


/**********************************************************************
  3.0) FORCE CLEAR (safe rerun)
**********************************************************************/
cap estimates clear
foreach v in ability_raw ability_z ability_hi ///
            ln_prod ln_piece ///
            tdate_all post_all event_week_all ///
            abil_sample ///
            fe_worker_pre fe_contract_pre ///
            treated_year_3B placebo_year_3B triplediff_sample_3C ///
            has_pre has_post spans any_pre any_post {
    cap drop `v'
}



/**********************************************************************
  3.1) PRE-TREATMENT WORKER ABILITY (AKM-style)
  Key idea:
    - Use ONLY pre–June 1 observations (post_all==0)
    - Pool across years/provinces to maximize ability coverage
    - Regress ln_prod on controls, absorbing worker_id and contract_id
    - Save worker FE as ability_raw, standardize to ability_z
**********************************************************************/

* ---- 3.1A Construct a generic June 1 cutoff for EACH observation's year ----
* (works for any year, any province)
cap drop tdate_all post_all event_week_all
gen double tdate_all = mdy(6,1,year)
format tdate_all %td

gen byte post_all = (date >= tdate_all) if !missing(date) & !missing(tdate_all)
gen int  event_week_all = floor((date - tdate_all)/7) if !missing(date) & !missing(tdate_all)

* ---- 3.1B Outcome for ability model: ln(prod) ----
* (you said prod is trees/hour; if instead prod is trees/day adjust here)
cap drop ln_prod
gen double ln_prod = ln(prod) if prod>0

* ---- 3.1C Define ability estimation sample ----
* Conservative default: use PRE only, apply SAME floors as baseline if you want.
* You can tighten/loosen ABILITY_PREW if desired.
cap macro drop ABILITY_PREW
global ABILITY_PREW 12    // use up to 12 pre weeks (relative to June 1) for ability

cap drop abil_sample
gen byte abil_sample = (post_all==0) ///
    & inrange(event_week_all, -$ABILITY_PREW, -1) ///
    & !missing(ln_prod)

* Optional: apply the same floors as baseline if they exist
* (If $SAMPLE_OK exists from Section 2, use it; otherwise do nothing.)
capture confirm global SAMPLE_OK
if _rc==0 {
    replace abil_sample = abil_sample & ($SAMPLE_OK)
}

count if abil_sample
summ ln_prod if abil_sample

* ---- 3.1D Estimate ability (worker FE) using pre sample only ----
* IMPORTANT: absorb order matters for savefe mapping: worker first, contract second.
*
* We include the same flexible controls you've standardized in Section 1/2.
* If you want zero-rain dummy + quadratic rain + quadratic temp/wind: already built.
*
* NOTE: savefe will create __hdfe1__ (worker FE) and __hdfe2__ (contract FE)
* because worker_id is first in absorb().
cap noi reghdfe ln_prod ///
    c.experience##c.experience ///
    c.temp##c.temp c.wind_spd##c.wind_spd ///
    precip0 c.precip_pos##c.precip_pos ///
    if abil_sample, ///
    absorb(worker_id contract_id, savefe) ///
    vce(cluster contract_id)

* ---- 3.1E Pull saved FE safely ----
* Worker FE should be __hdfe1__ given absorb(worker_id contract_id, savefe)
cap drop fe_worker_pre fe_contract_pre ability_raw
cap confirm variable __hdfe1__
if _rc {
    di as error "[Section 3.1] ERROR: reghdfe did not produce __hdfe1__. Check reghdfe version/savefe output."
    exit 198
}
gen double fe_worker_pre = __hdfe1__

cap confirm variable __hdfe2__
if !_rc gen double fe_contract_pre = __hdfe2__

gen double ability_raw = fe_worker_pre

* Standardize
summ ability_raw if !missing(ability_raw)
gen double ability_z = (ability_raw - r(mean)) / r(sd) if !missing(ability_raw)

* Optional binary split (median)
cap drop ability_hi
xtile ability_hi = ability_raw if !missing(ability_raw), nq(2)
replace ability_hi = ability_hi - 1

summ ability_raw ability_z
tab ability_hi


/**********************************************************************
  3.2) HETEROGENEOUS MW EFFECTS
  We implement heterogeneity with explicitly-generated interaction terms.
  This avoids factor-variable corner cases and makes omissions diagnosable.

  3.2A) Cross-province (BC vs AB/ON) in treated year (2A design)
        Outcome: $Y (set in Section 2 control panel)
        Model: Y = treat + post + treat×post + ability + ability×post + ability×treat + ability×treat×post + controls + FE

  REQUIREMENTS:
    - from Section 2A: treated_BC_2A, post_2A, event_week_2A, $IF_2A, $SAMPLE_OK, $X, $ABS
**********************************************************************/

* ---- 3.2A.0 Recreate 2A objects if missing (defensive) ----
capture confirm variable treated_BC_2A
if _rc {
    di as error "[Section 3.2A] treated_BC_2A not found. Run Section 2A first (or define treated_BC_2A/post_2A/event_week_2A)."
    exit 198
}
capture confirm variable post_2A
if _rc {
    di as error "[Section 3.2A] post_2A not found. Run Section 2A first."
    exit 198
}
capture confirm global IF_2A
if _rc {
    di as error "[Section 3.2A] global IF_2A not found. Run Section 2A first."
    exit 198
}

* ---- 3.2A.1 Build explicit interactions ----
cap drop treat_post ability_post ability_treat ability_treat_post
gen byte   treat_post        = (treated_BC_2A==1 & post_2A==1) if $IF_2A
gen double ability_post      = ability_z * post_2A            if $IF_2A & !missing(ability_z)
gen double ability_treat     = ability_z * treated_BC_2A      if $IF_2A & !missing(ability_z)
gen double ability_treat_post= ability_z * treat_post         if $IF_2A & !missing(ability_z)

* ---- 3.2A.2 Run heterogeneity regression (AB/ON control) ----
* IMPORTANT: Do not accidentally kill post support.
* We include ability_z but do NOT require ability_z for *everyone* unless you want that restriction.
* If you want to require it, keep !missing(ability_z); otherwise you can include missing as zero (not recommended).
cap noi reghdfe $Y ///
    treated_BC_2A post_2A treat_post ///
    ability_z ability_post ability_treat ability_treat_post ///
    $X ///
    if $IF_2A ///
       & inrange(event_week_2A, -$WPRE_2A, $WPOST_2A) ///
       & !missing($Y) ///
       & ($SAMPLE_OK) ///
       & !missing(ability_z), ///
    absorb($ABS) ///
    vce(cluster contract_id)
estimates store het_2A

* Quick diagnostic: confirm BOTH pre and post exist in estimating sample
tab post_2A if e(sample)
tab treated_BC_2A if e(sample)


/**********************************************************************
  3.2B) Within-BC placebo design heterogeneity (2B design)
  Model: Y = treated_year + post_y + treated_year×post_y
            + ability + ability×post + ability×treated_year + ability×treated_year×post
            + controls + FE
**********************************************************************/

* Defensive checks
capture confirm variable bc_placebo_sample
if _rc {
    di as error "[Section 3.2B] bc_placebo_sample not found. Run Section 2B first."
    exit 198
}
capture confirm variable treated_year
if _rc {
    di as error "[Section 3.2B] treated_year not found. Run Section 2B first."
    exit 198
}
capture confirm variable post_y
if _rc {
    di as error "[Section 3.2B] post_y not found. Run Section 2B first."
    exit 198
}
capture confirm variable event_week_y
if _rc {
    di as error "[Section 3.2B] event_week_y not found. Run Section 2B first."
    exit 198
}

cap drop ty_post ability_post_BC ability_ty ability_ty_post
gen byte   ty_post        = (treated_year==1 & post_y==1) if bc_placebo_sample
gen double ability_post_BC= ability_z * post_y           if bc_placebo_sample & !missing(ability_z)
gen double ability_ty     = ability_z * treated_year     if bc_placebo_sample & !missing(ability_z)
gen double ability_ty_post= ability_z * ty_post          if bc_placebo_sample & !missing(ability_z)

cap noi reghdfe $Y ///
    treated_year post_y ty_post ///
    ability_z ability_post_BC ability_ty ability_ty_post ///
    $X ///
    if bc_placebo_sample ///
       & inrange(event_week_y, $PRE_MIN_2B, $POST_MAX_2B) ///
       & !missing($Y) ///
       & ($SAMPLE_OK) ///
       & !missing(ability_z), ///
    absorb($ABS) ///
    vce(cluster contract_id)
estimates store het_2B

tab post_y if e(sample)
tab treated_year if e(sample)


/**********************************************************************
  3.2C) OPTIONAL TRIPLE DIFF heterogeneity
  Sample: BC + control province, treated year + placebo years
  Model includes:
    treated_BC × treated_year × post_y and ability interactions
**********************************************************************/

* You only want AB control in triple-diff typically; adjust if needed.
* Build the combined sample indicator using what Section 2 already defined.
cap drop triplediff_sample_3C treated_BC_3C
gen byte triplediff_sample_3C = (province=="BC" | province=="$CONTROL_2A") ///
    & (treated_year==1 | placebo_year==1)

gen byte treated_BC_3C = (province=="BC") if triplediff_sample_3C

* Need a post variable defined by year-specific cutoff (post_y already is for BC years only in 2B);
* For triple-diff, define post_all (already built in 3.1) and use that.
cap drop post_3C
gen byte post_3C = post_all if triplediff_sample_3C

* Core triple interaction
cap drop DDD
gen byte DDD = treated_BC_3C*treated_year*post_3C if triplediff_sample_3C

* Ability × DDD
cap drop ability_DDD
gen double ability_DDD = ability_z*DDD if triplediff_sample_3C & !missing(ability_z)

cap noi reghdfe $Y ///
    i.treated_BC_3C##i.treated_year##i.post_3C ///
    c.ability_z##i.treated_BC_3C##i.treated_year##i.post_3C ///
    $X ///
    if triplediff_sample_3C ///
       & inrange(event_week_all, -$ABILITY_PREW, $POST_MAX_2B) ///
       & !missing($Y) ///
       & ($SAMPLE_OK) ///
       & !missing(ability_z), ///
    absorb($ABS) ///
    vce(cluster contract_id)
estimates store het_DDD


/**********************************************************************
  3.3) SELECTION TESTS
  A) Within-BC placebo: does ability predict working post June 1 more
     in treated year vs placebo years?

  B) Cross-province treated year: does ability predict working post more
     in BC vs AB (or ON)?
**********************************************************************/

* ---- 3.3A Within-BC placebo selection (worker-year) ----
preserve
    keep if bc_placebo_sample
    keep worker_id year post_y ability_z

    bysort worker_id year: egen any_post = max(post_y==1)
    bysort worker_id year: egen any_pre  = max(post_y==0)
    bysort worker_id year: keep if _n==1
    keep if any_pre==1
    * treated_year already exists
    reg any_post c.ability_z##i.treated_year, vce(cluster worker_id)
restore

* ---- 3.3B Cross-province treated-year selection (worker-province) ----
preserve
    keep if $IF_2A
    keep worker_id province post_2A ability_z

    bysort worker_id province: egen any_post = max(post_2A==1)
    bysort worker_id province: egen any_pre  = max(post_2A==0)
    bysort worker_id province: keep if _n==1
    keep if any_pre==1

    gen byte treatedBC = (province=="BC")
    reg any_post c.ability_z##i.treatedBC, vce(cluster worker_id)
	
	
	
	
	










/*

/**********************************************************************
  (2) BASELINE ANALYSIS (RERUNNABLE)
  Two baseline designs:

    2A) Cross-province: BC vs AB/ON in treated year
        - single cutoff (June 1 of treated year)
        - DiD + event study + plot

    2B) Within-BC placebo: treated year vs placebo years (BC only)
        - year-specific cutoff (June 1 each year)
        - DiD + event study + plot

  Key features:
    - One consistent control panel
    - Floors (hours, prod, earnings proxy) to avoid tiny/partial-day artifacts
    - Auto event-window choice using TAU rule (default 2/3 of "core" mass)
**********************************************************************/


/**********************************************************************
  2.0) RERUN SAFETY: clear estimates + dynamic vars created in Section 2
**********************************************************************/
cap estimates clear
cap label drop ew_lbl

foreach v in sample_2A treated_BC post event_week event_week_s ///
            placebo_year treated_year bc_placebo_sample group_str group_s event_week_sy {
    cap drop `v'
}
foreach m in X ABS SAMPLE_OK IFBASE TDATE EW_BASE EW_BASE_BC ///
            WPRE WPOST PRE_MIN_BC POST_MAX_BC {
    cap macro drop `m'
}


/**********************************************************************
  2.1) CONTROL PANEL — EDIT THESE ONLY
**********************************************************************/

* --------- Outcome ---------
* Options created in Section 1: ln_incentive  ln_piece  lprod  tot_hrs_day
global Y "ln_incentive"

* --------- Controls toggles ---------
global USE_WEATHER  1
global USE_EXPER    1
global USE_HOURS    1
global USE_DAYFLAGS 1

* --------- Fixed effects ---------
* "worker_contract" -> worker_id + contract_id (2WFE)
* "wc_only"         -> wc_id (worker×contract)
* "worker_only"     -> worker only
* "contract_only"   -> contract only
global FE "worker_contract"

* --------- Floors (avoid tiny-day artifacts) ---------
global USE_FLOORS 1
global MIN_HRS   2
global MIN_PROD  50
global MIN_EARN  30
global USE_MIN_PROD 1

* --------- Auto window selection (TAU rule) ---------
global TAU   0.6667       // default = 2/3
global CORE0 0            // "core" window for mass benchmark (weeks)
global CORE1 1
global MAXW  30           // max weeks to search outward

* Manual fallbacks if you turn auto off
global AUTO_WIN_2A 1
global AUTO_WIN_2B 1
global WPRE_MAN  3
global WPOST_MAN 8
global PRE_MIN_BC_MAN  -10
global POST_MAX_BC_MAN  12

* --------- 2A design choices ---------
global TY_2A 2019
global CONTROL_2A "AB"     // AB, ON, or ABON (ABON uses manual window)

* --------- 2B design choices ---------
global TY_2B 2019
global BCPLACEBO "2016 2017"


/**********************************************************************
  2.2) BUILD CONTROL VECTOR ($X), FE VECTOR ($ABS), SAMPLE FLOORS ($SAMPLE_OK)
**********************************************************************/

* --- Controls ---
global X ""
if $USE_HOURS==1    global X "$X tot_hrs_day"
if $USE_DAYFLAGS==1 global X "$X multi_contract_day"
if $USE_EXPER==1    global X "$X c.experience##c.experience"
if $USE_WEATHER==1  global X "$X c.temp##c.temp c.wind_spd##c.wind_spd precip0 c.precip_pos##c.precip_pos"

* --- Fixed effects ---
global ABS ""
if "$FE"=="worker_contract" global ABS "worker_id contract_id"
if "$FE"=="wc_only"         global ABS "wc_id"
if "$FE"=="worker_only"     global ABS "worker_id"
if "$FE"=="contract_only"   global ABS "contract_id"

* --- Floors (as an if-condition; no dropping) ---
global SAMPLE_OK "1==1"
if $USE_FLOORS==1 {
    global SAMPLE_OK "$SAMPLE_OK & tot_hrs_day >= $MIN_HRS"
    global SAMPLE_OK "$SAMPLE_OK & piece > 0 & prod > 0"
    global SAMPLE_OK "$SAMPLE_OK & (piece*prod) >= $MIN_EARN"
    if $USE_MIN_PROD==1 global SAMPLE_OK "$SAMPLE_OK & prod >= $MIN_PROD"
}

display as text "[DEBUG] Outcome Y:     $Y"
display as text "[DEBUG] Controls X:    $X"
display as text "[DEBUG] Absorb ABS:    $ABS"
display as text "[DEBUG] Floors OK:     $SAMPLE_OK"


/**********************************************************************
  2.3) DEPENDENCIES FOR WINDOW SELECTION
  We require a short province string for collapse/reshape robustness.
  (province_s created in Section 1; recreate safely if not present)
**********************************************************************/
capture confirm variable province_s
if _rc {
    cap drop province_s
    gen str2 province_s = province
}
/**********************************************************************
  2A) CROSS-PROVINCE BASELINE: BC vs AB/ON in treated year
      - Self-contained: defines its own sample (no reliance on sample_2A)
      - Auto window selection using TAU rule
      - DiD + Event study + coefplot export
**********************************************************************/

* --------------------------
* 2A.0 SETTINGS
* --------------------------
global TY_2A 2019
global CONTROL_2A "AB"      // AB, ON, or ABON

* Auto-window parameters (2-group only; ABON uses manual)
global AUTO_WIN_2A 1
global TAU   0.6667
global CORE0 0
global CORE1 1
global MAXW  30

* Manual fallback
global WPRE_MAN  3
global WPOST_MAN 8


* --------------------------
* 2A.1 DEFINE SAMPLE FILTER (STRING) + EVENT TIME
* --------------------------
cap macro drop IF_2A
global IF_2A "(year==$TY_2A)"

if "$CONTROL_2A"=="AB"   global IF_2A `"$IF_2A & (province=="BC" | province=="AB")"'
if "$CONTROL_2A"=="ON"   global IF_2A `"$IF_2A & (province=="BC" | province=="ON")"'
if "$CONTROL_2A"=="ABON" global IF_2A `"$IF_2A & (province=="BC" | province=="AB" | province=="ON")"'

global TDATE_2A = td(01jun$TY_2A)

cap drop treated_BC_2A post_2A event_week_2A
gen byte treated_BC_2A = (province=="BC") if $IF_2A
gen byte post_2A       = (date >= $TDATE_2A) if $IF_2A
gen int  event_week_2A = floor((date - $TDATE_2A)/7) if $IF_2A

* quick support check
tab province if $IF_2A
summ event_week_2A if $IF_2A


/**********************************************************************
  2A.2 AUTO WINDOW SELECTION (TAU rule) OR MANUAL WINDOW
  - Goal: choose the largest symmetric-ish window around June 1 such that
    weekly counts in BOTH groups stay above TAU × (core mean count)
    for a contiguous pre block ending at -1 and post block starting at 0.

  USER CONTROLS:
    global AUTO_WIN_2A   1  (set 0 to force manual)
    global TAU           0.3333   (default 1/3; change as desired)
    global CORE0         0        (core window start)
    global CORE1         1        (core window end)
    global MAXW          30       (max weeks to search outward)
    global WPRE_MAN      3        (manual pre weeks if auto off/fails)
    global WPOST_MAN     8        (manual post weeks if auto off/fails)

  NOTE:
    - This block is self-contained and relies only on:
      event_week_2A, province, $IF_2A, $CONTROL_2A
**********************************************************************/

* ABON is 3-group; this auto-window is written for 2 groups. Force manual.
if "$CONTROL_2A"=="ABON" global AUTO_WIN_2A 0

* If user did not define TAU elsewhere, default to 1/3
capture confirm global TAU
if _rc global TAU 0.3333

* If user did not define manual windows, define defaults
capture confirm global WPRE_MAN
if _rc global WPRE_MAN 3
capture confirm global WPOST_MAN
if _rc global WPOST_MAN 8

* Main branch
if $AUTO_WIN_2A==1 {

    preserve
        keep if $IF_2A
        keep if inrange(event_week_2A, -$MAXW, $MAXW)
        keep if province=="BC" | province=="$CONTROL_2A"

        * Counts per week×province (count date = always non-missing)
        collapse (count) N=date, by(event_week_2A province)
        reshape wide N, i(event_week_2A) j(province) string

        * Identify the two count series
        capture confirm variable NBC
        if _rc {
            di as error "[2A AUTO WIN] FAIL: After reshape, NBC not found. Check province values."
            list in 1/20
            restore
            global WPRE_2A  $WPRE_MAN
            global WPOST_2A $WPOST_MAN
            exit
        }

        local vctrl = "N$CONTROL_2A"
        capture confirm variable `vctrl'
        if _rc {
            di as error "[2A AUTO WIN] FAIL: After reshape, `vctrl' not found. CONTROL_2A=$CONTROL_2A"
            list in 1/20
            restore
            global WPRE_2A  $WPRE_MAN
            global WPOST_2A $WPOST_MAN
            exit
        }

        * Core means
        quietly summarize NBC if inrange(event_week_2A,$CORE0,$CORE1)
        scalar coreBC = r(mean)

        quietly summarize `vctrl' if inrange(event_week_2A,$CORE0,$CORE1)
        scalar coreCT = r(mean)

        * Thresholds
        scalar thrBC = $TAU * coreBC
        scalar thrCT = $TAU * coreCT

        di as text "[2A AUTO WIN] Params: TAU=" %6.4f $TAU " CORE=[" $CORE0 "," $CORE1 "] MAXW=" $MAXW
        di as text "[2A AUTO WIN] Core means: BC=" %9.2f coreBC "  CTRL=" %9.2f coreCT
        di as text "[2A AUTO WIN] Thresholds: BC>=" %9.2f thrBC "  CTRL>=" %9.2f thrCT

        * Passing weeks: both groups above threshold
        gen byte pass = (NBC >= thrBC & `vctrl' >= thrCT)
        sort event_week_2A

        * PRE contiguous block ending at -1
        gen byte okpre = pass & (event_week_2A<=-1)
        gen byte cont_pre = .
        replace cont_pre = okpre if event_week_2A==-1
        forvalues k=2/$MAXW {
            local wk = -`k'
            quietly replace cont_pre = okpre & cont_pre[_n+1] if event_week_2A==`wk'
        }

        * POST contiguous block starting at 0
        gen byte okpost = pass & (event_week_2A>=0)
        gen byte cont_post = .
        replace cont_post = okpost if event_week_2A==0
        forvalues k=1/$MAXW {
            quietly replace cont_post = okpost & cont_post[_n-1] if event_week_2A==`k'
        }

        quietly summarize event_week_2A if cont_pre==1, meanonly
        scalar PRE_sug  = r(min)

        quietly summarize event_week_2A if cont_post==1, meanonly
        scalar POST_sug = r(max)

        * Diagnostics for failure modes
        quietly count if cont_pre==1
        scalar n_pre_ok = r(N)
        quietly count if cont_post==1
        scalar n_post_ok = r(N)

        di as text "[2A AUTO WIN] Contiguous support: pre_weeks_ok=" n_pre_ok " post_weeks_ok=" n_post_ok
        di as text "[2A AUTO WIN] Suggested bounds: PRE_sug=" PRE_sug " POST_sug=" POST_sug

    restore

    * ---- SAFE ASSIGNMENT WITH INFORMATIVE MESSAGES ----
    if missing(PRE_sug) {
        di as error "[2A AUTO WIN] FAIL: No PRE window meets TAU=" %6.4f $TAU ///
                    " for a contiguous block ending at week -1 (BC vs $CONTROL_2A)."
        di as error "             -> Using MANUAL window instead: [-$WPRE_MAN, $WPOST_MAN]."
        global WPRE_2A  $WPRE_MAN
        global WPOST_2A $WPOST_MAN
    }
    else if missing(POST_sug) {
        di as error "[2A AUTO WIN] FAIL: No POST window meets TAU=" %6.4f $TAU ///
                    " for a contiguous block starting at week 0 (BC vs $CONTROL_2A)."
        di as error "             -> Using MANUAL window instead: [-$WPRE_MAN, $WPOST_MAN]."
        global WPRE_2A  $WPRE_MAN
        global WPOST_2A $WPOST_MAN
    }
    else {
        global WPRE_2A  = -PRE_sug
        global WPOST_2A =  POST_sug
        di as text  "[2A AUTO WIN] SUCCESS: Using window [-$WPRE_2A, $WPOST_2A]."
    }
}
else {
    * Manual window requested
    global WPRE_2A  $WPRE_MAN
    global WPOST_2A $WPOST_MAN
    di as text "[2A AUTO WIN] OFF: Using manual window [-$WPRE_2A, $WPOST_2A]."
}

display as text "[2A] FINAL window: [-$WPRE_2A, $WPOST_2A]  (AUTO=$AUTO_WIN_2A, TAU=$TAU)"
/**********************************************************************
  2A EVENT STUDY — HARD SUPPORT FILTER (GUARANTEED)
  Drops any event_week where the control province has zero support
  in the *actual regression sample*.
**********************************************************************/

* 2A: shifted week
cap drop event_week_s_2A okweek_2A
gen int event_week_s_2A = event_week_2A + 100 if $IF_2A ///
    & inrange(event_week_2A, -$WPRE_2A, $WPOST_2A)

* Build week-support table inside the EXACT regression sample
preserve
    keep if $IF_2A ///
        & inrange(event_week_2A, -$WPRE_2A, $WPOST_2A) ///
        & ($SAMPLE_OK) ///
        & !missing($Y)

    keep if province=="BC" | province=="$CONTROL_2A"

    collapse (count) N=date, by(event_week_2A province)
    reshape wide N, i(event_week_2A) j(province) string

    gen byte ok = (NBC>0 & N$CONTROL_2A>0)

    keep event_week_2A ok
    tempfile wkOK
    save `wkOK', replace
restore

* Merge and enforce: any week not ok => event_week_s_2A missing
merge m:1 event_week_2A using `wkOK'
gen byte okweek_2A = (_merge==3 & ok==1)
drop _merge ok
replace event_week_s_2A = . if okweek_2A==0

* Confirm: there should be NO weeks with AB=0 remaining
tab event_week_2A province if $IF_2A & !missing(event_week_s_2A) & ($SAMPLE_OK) & !missing($Y)

* Base week: use -1 if present
cap macro drop EW_BASE_2A
quietly count if $IF_2A & !missing(event_week_s_2A) & event_week_2A==-1
if r(N)>0 global EW_BASE_2A 99
else {
    di as error "[2A] ERROR: week -1 not in identified support. Pick a different window."
    exit 198
}

* Labels for included weeks only
cap label drop ew_lbl
levelsof event_week_s_2A if !missing(event_week_s_2A), local(levs)
foreach k of local levs {
    local val = `k' - 100
    label define ew_lbl `k' "`val'", add
}
label values event_week_s_2A ew_lbl
/**********************************************************************
  2A EVENT STUDY — SIMPLE SUPPORT FILTER (NO PRESERVE/COLLAPSE/MERGE)
  Goal: Drop any event-week where the CONTROL group has zero observations
        in the *actual estimating sample* (after floors + !missing(Y)).
  This fixes the -6..-4 problem when AB has 0 obs there.
**********************************************************************/

* 0) Clean up
cap drop in_esamp_2A bc_obs_2A ctrl_obs_2A bcN_2A ctrlN_2A okweek_2A
cap drop event_week_s_2A

* 1) Define the *actual event-study estimation sample* row-by-row
gen byte in_esamp_2A = ///
    ($IF_2A) ///
    & inrange(event_week_2A, -$WPRE_2A, $WPOST_2A) ///
    & ($SAMPLE_OK) ///
    & !missing($Y) ///
    & (province=="BC" | province=="$CONTROL_2A")

* 2) Count BC and CONTROL observations by event week within that sample
gen byte bc_obs_2A   = in_esamp_2A & (province=="BC")
gen byte ctrl_obs_2A = in_esamp_2A & (province=="$CONTROL_2A")

bysort event_week_2A: egen int bcN_2A   = total(bc_obs_2A)
bysort event_week_2A: egen int ctrlN_2A = total(ctrl_obs_2A)

* 3) Identify usable weeks: both groups have >0 obs
cap drop okweek_2A
gen byte okweek_2A = in_esamp_2A & (bcN_2A>0) & (ctrlN_2A>0)

* 4) Create shifted event time ONLY for usable weeks
gen int event_week_s_2A = event_week_2A + 100 if okweek_2A

* 5) Debug check: this table MUST show no AB=0 weeks remaining
tab event_week_2A province if !missing(event_week_s_2A)

* 6) Base week (week -1) must exist among usable weeks
cap macro drop EW_BASE_2A
quietly count if !missing(event_week_s_2A) & event_week_2A==-1
if r(N)==0 {
    di as error "[2A] FAIL: week -1 not available after support filter."
    di as error "     -> Reduce WPRE_2A or loosen floors, or switch to BC-placebo design."
    exit 198
}
global EW_BASE_2A 99

* 7) Labels only for included weeks
cap label drop ew_lbl
levelsof event_week_s_2A if !missing(event_week_s_2A), local(levs)
foreach k of local levs {
    local val = `k' - 100
    label define ew_lbl `k' "`val'", add
}
label values event_week_s_2A ew_lbl
/**********************************************************************
  2A.5 Event study regression (FIXED)
**********************************************************************/

* Enforce restricted-support weeks
replace event_week_s_2A = . if okweek_2A==0

* Guard: ensure base week (-1 => 99) exists
quietly count if $IF_2A & !missing(event_week_s_2A) & event_week_s_2A==99 ///
    & (province=="BC" | province=="$CONTROL_2A") & ($SAMPLE_OK) & !missing($Y)
if r(N)==0 {
    di as error "[2A] FAIL: Base week -1 not present after support restriction."
    exit 198
}

* Event study regression: restrict factor-variable scope via if-condition
reghdfe $Y ///
    ib(99).event_week_s_2A##i.treated_BC_2A ///
    $X ///
    if $IF_2A ///
       & (province=="BC" | province=="$CONTROL_2A") ///
       & !missing(event_week_s_2A) ///
       & ($SAMPLE_OK) ///
       & !missing($Y), ///
    absorb($ABS) ///
    vce(cluster contract_id)
estimates store es_2A


/**********************************************************************
  2A.6 Plot (FIXED)
**********************************************************************/

* Keep only weeks actually used in estimation
levelsof event_week_s_2A if e(sample), local(actual_levs)

local keepcoef ""
local clabel   ""
foreach k of local actual_levs {
    local keepcoef `keepcoef' `k'.event_week_s_2A#1.treated_BC_2A
    local val = `k' - 100
    local clabel `clabel' `k'.event_week_s_2A#1.treated_BC_2A = "`val'"
}

* Vertical line between -1 and 0 based on plotted sequence
preserve
    keep if e(sample)
    keep event_week_s_2A
    duplicates drop
    sort event_week_s_2A
    gen idx = _n
    quietly summarize idx if event_week_s_2A==100, meanonly
    local line_pos = r(mean) - 0.5
restore

coefplot es_2A, ///
    keep(`keepcoef') ///
    vertical ///
    omitted baselevels ///
    coeflabel(`clabel') ///
    yline(0, lpattern(dash)) ///
    xline(`line_pos', lwidth(vthin)) ///
    title("Event study: $Y ($TY_2A)", size(medium)) ///
    subtitle("BC vs $CONTROL_2A | Restricted Support Only", size(small)) ///
    graphregion(color(white))
	
	
/**********************************************************************
  2B) WITHIN-BC PLACEBO: treated year vs placebo years (BC only)
**********************************************************************/

* --------------------------
* 2B.1 Construct sample + group labels
* --------------------------
cap drop placebo_year treated_year bc_placebo_sample group_str group_s event_week_sy
gen byte placebo_year = 0
foreach y of numlist $BCPLACEBO {
    replace placebo_year = 1 if year==`y'
}
gen byte treated_year = (year==$TY_2B)

gen byte bc_placebo_sample = (province_s=="BC") & (treated_year==1 | placebo_year==1)

gen str8 group_str = cond(treated_year==1, "treated", "placebo") if bc_placebo_sample==1
gen str8 group_s   = group_str if bc_placebo_sample==1   // safe copy for reshape

* --------------------------
* 2B.2 Choose window (TAU rule) or manual
* Uses universal event_week_y defined in Section 1
* --------------------------
if $AUTO_WIN_2B==1 {

    preserve
        keep if bc_placebo_sample==1
        keep if inrange(event_week_y, -$MAXW, $MAXW)

        collapse (count) N=event_week_y, by(event_week_y group_s)
        reshape wide N, i(event_week_y) j(group_s) string

        quietly summarize Ntreated if inrange(event_week_y,$CORE0,$CORE1)
        scalar coreT = r(mean)
        quietly summarize Nplacebo if inrange(event_week_y,$CORE0,$CORE1)
        scalar coreP = r(mean)

        scalar thrT = $TAU * coreT
        scalar thrP = $TAU * coreP

        gen byte pass = (Ntreated >= thrT & Nplacebo >= thrP)

        sort event_week_y

        * PRE contiguous block ending at -1
        gen byte okpre = pass & (event_week_y<=-1)
        gen byte cont_pre = .
        replace cont_pre = okpre if event_week_y==-1
        forvalues k=2/$MAXW {
            local wk = -`k'
            quietly replace cont_pre = okpre & cont_pre[_n+1] if event_week_y==`wk'
        }

        * POST contiguous block starting at 0
        gen byte okpost = pass & (event_week_y>=0)
        gen byte cont_post = .
        replace cont_post = okpost if event_week_y==0
        forvalues k=1/$MAXW {
            quietly replace cont_post = okpost & cont_post[_n-1] if event_week_y==`k'
        }

        quietly summarize event_week_y if cont_pre==1, meanonly
        scalar PRE_sug_BC = r(min)
        quietly summarize event_week_y if cont_post==1, meanonly
        scalar POST_sug_BC = r(max)

    restore

    global PRE_MIN_BC  = PRE_sug_BC
    global POST_MAX_BC = POST_sug_BC
}
else {
    global PRE_MIN_BC  $PRE_MIN_BC_MAN
    global POST_MAX_BC $POST_MAX_BC_MAN
}

display as text "[2B] Using window [$PRE_MIN_BC, $POST_MAX_BC] (TY=$TY_2B vs placebo=$BCPLACEBO, TAU=$TAU)"

* Shifted week var (avoid negatives in factor vars)
gen int event_week_sy = event_week_y + 100 if bc_placebo_sample==1 & inrange(event_week_y,$PRE_MIN_BC,$POST_MAX_BC)

* Labels
cap label drop ew_lbl
local lo2 = 100 + $PRE_MIN_BC
local hi2 = 100 + $POST_MAX_BC
forval i = `lo2'/`hi2' {
    local val = `i' - 100
    label define ew_lbl `i' "`val'", add
}
label values event_week_sy ew_lbl
global EW_BASE_BC 99

* --------------------------
* 2B.3 DiD regression (within BC)
* --------------------------
reghdfe $Y ///
    i.treated_year##i.post_y ///
    $X ///
    if bc_placebo_sample==1 & inrange(event_week_y,$PRE_MIN_BC,$POST_MAX_BC) & !missing($Y) & ($SAMPLE_OK), ///
    absorb($ABS) ///
    vce(cluster contract_id)
estimates store did_2B

* --------------------------
* 2B.4 Event-study regression (within BC)
* --------------------------
reghdfe $Y ///
    ib($EW_BASE_BC).event_week_sy##i.treated_year ///
    $X ///
    if bc_placebo_sample==1 & !missing(event_week_sy) & !missing($Y) & ($SAMPLE_OK), ///
    absorb($ABS) ///
    vce(cluster contract_id)
estimates store es_2B

* --------------------------
* 2B.5 Event-study plot export
* --------------------------
local clabel2 ""
forval i = `lo2'/`hi2' {
    local val = `i' - 100
    local clabel2 `clabel2' `i'.event_week_sy#1.treated_year = "`val'" 1.treated_year#`i'.event_week_sy = "`val'"
}
local wpre_bc = -$PRE_MIN_BC
local line_pos2 = `wpre_bc' + 0.5

coefplot, ///
    keep(*.event_week_sy#1.treated_year 1.treated_year#*.event_week_sy) ///
    omitted baselevels vertical ///
    coeflabel(`clabel2') ///
    yline(0, lpattern(dash)) ///
    xline(`line_pos2', lwidth(medium)) ///
    title("Event study: $Y (BC only)", size(medium)) ///
    subtitle("Treated $TY_2B vs placebo ($BCPLACEBO) | Window [$PRE_MIN_BC,$POST_MAX_BC] | TAU=$TAU", size(small)) ///
    ytitle("Effect relative to week -1") ///
    xtitle("Weeks relative to June 1 cutoff") ///
    graphregion(color(white)) xsize(12) ysize(7) ///
    xlabel(, labsize(small)) grid(none)

graph export "event_2B_${Y}_${TY_2B}.png", replace






*/

/**********************************************************************
  (2) BASELINE ANALYSIS (RERUNNABLE CONTROL PANEL)
**********************************************************************/

***********************************************************************
* 2.0) FORCE CLEAR DYNAMIC VARS / MACROS (SAFE RERUN)
***********************************************************************/
cap estimates clear
cap label drop ew_lbl

foreach v in treated_BC post event_week event_week_s ///
            placebo_year treated_year bc_placebo_sample group_str event_week_sy sample_2A {
    cap drop `v'
}
foreach m in IFBASE TDATE EW_BASE EW_BASE_BC X ABS SAMPLE_OK {
    cap macro drop `m'
}


***********************************************************************
* 2.1) USER SWITCHES (EDIT THESE)
***********************************************************************

* Dependent variable options: ln_incentive  ln_piece  lprod  tot_hrs_day
global Y "ln_incentive"

* Controls toggles
global USE_WEATHER  1
global USE_EXPER    1
global USE_HOURS    1
global USE_DAYFLAGS 1

* Fixed effects choice:
*   "worker_contract" -> absorb(worker_id contract_id)
*   "wc_only"         -> absorb(wc_id)
*   "worker_only"     -> absorb(worker_id)
*   "contract_only"   -> absorb(contract_id)
global FE "worker_contract"

* Floors to avoid tiny/partial-day artifacts
global USE_FLOORS 1
global MIN_HRS   2      // minimum total hours worked that day (tot_hrs_day)
global MIN_PROD  50     // minimum productivity (prod = trees/hour)
global MIN_EARN  30     // minimum incentive earnings proxy: piece*prod
global USE_MIN_PROD 1   // keep prod floor on/off

* Automatic window selection by TAU rule
global AUTO_WIN_2A 1
global AUTO_WIN_2B 1

* TAU rule settings
global TAU   0.6667     // change freely (e.g., 0.50, 0.75)
global CORE0 0          // core weeks for baseline mass (e.g., 0..1)
global CORE1 1
global MAXW  30         // search range for suggested window

* Manual windows if AUTO_WIN_* = 0
global WPRE_MAN  3
global WPOST_MAN 8

global PRE_MIN_BC_MAN  -10
global POST_MAX_BC_MAN  12


***********************************************************************
* 2.2) BUILD COMMON CONTROL VECTOR $X AND FE VECTOR $ABS
***********************************************************************

* Controls string used in all baseline regressions
global X ""
if $USE_HOURS==1    global X "$X tot_hrs_day"
if $USE_DAYFLAGS==1 global X "$X multi_contract_day"
if $USE_EXPER==1    global X "$X c.experience##c.experience"
if $USE_WEATHER==1  global X "$X c.temp##c.temp c.wind_spd##c.wind_spd precip0 c.precip_pos##c.precip_pos"

* Absorb string (fixed effects)
global ABS ""
if "$FE"=="worker_contract" global ABS "worker_id contract_id"
if "$FE"=="wc_only"         global ABS "wc_id"
if "$FE"=="worker_only"     global ABS "worker_id"
if "$FE"=="contract_only"   global ABS "contract_id"

* Inclusion condition used everywhere (NO dropping obs; this is an if-condition)
global SAMPLE_OK "1==1"
if $USE_FLOORS==1 {
    global SAMPLE_OK "$SAMPLE_OK & tot_hrs_day >= $MIN_HRS"
    global SAMPLE_OK "$SAMPLE_OK & piece > 0 & prod > 0"
    global SAMPLE_OK "$SAMPLE_OK & (piece*prod) >= $MIN_EARN"
    if $USE_MIN_PROD==1 global SAMPLE_OK "$SAMPLE_OK & prod >= $MIN_PROD"
}

display as text "[DEBUG] Y: $Y"
display as text "[DEBUG] Controls X:    $X"
display as text "[DEBUG] Absorb ABS:    $ABS"
display as text "[DEBUG] Sample floors: $SAMPLE_OK"


/**********************************************************************
  2A WINDOW SELECTION (TAU RULE; INLINE, ROBUST)
  Fix: Use province_s (short string) to avoid strL failures in collapse.
  - For ABON, fall back to manual because TAU rule is 2-group only.
**********************************************************************/

if "$CONTROL"=="ABON" global AUTO_WIN_2A 0

if $AUTO_WIN_2A==1 {

    preserve
        keep if sample_2A==1
        keep if inrange(event_week, -$MAXW, $MAXW)

        * Only BC + chosen control province
        keep if province_s=="BC" | province_s=="$CONTROL"

        * Counts by event_week × province_s
        collapse (count) N=event_week, by(event_week province_s)

        * Wide so we can apply thresholds to both groups
        reshape wide N, i(event_week) j(province_s) string

        * expects NBC and N<CONTROL>
        capture confirm variable NBC
        if _rc {
            di as error "Expected NBC after reshape, not found. Check province_s values."
            list in 1/15
            restore
            exit 198
        }

        local vctrl = "N$CONTROL"
        capture confirm variable `vctrl'
        if _rc {
            di as error "Expected `vctrl' after reshape, not found. CONTROL=$CONTROL. Check province_s values."
            list in 1/15
            restore
            exit 198
        }

        * Core-window average weekly counts
        quietly summarize NBC if inrange(event_week,$CORE0,$CORE1)
        scalar coreBC = r(mean)

        quietly summarize `vctrl' if inrange(event_week,$CORE0,$CORE1)
        scalar coreCT = r(mean)

        scalar thrBC = $TAU * coreBC
        scalar thrCT = $TAU * coreCT

        di as text "[DEBUG 2A] Core mean counts: BC=" coreBC " CTRL=" coreCT
        di as text "[DEBUG 2A] Thresholds: BC>=" thrBC " CTRL>=" thrCT

        * A week is admissible only if BOTH groups meet the thresholds
        gen byte pass = (NBC >= thrBC & `vctrl' >= thrCT)

        * PRE contiguous block ending at -1
        sort event_week
        gen byte okpre = pass & (event_week<=-1)
        gen byte cont_pre = .
        replace cont_pre = okpre if event_week==-1

        forvalues k=2/$MAXW {
            local wk = -`k'
            quietly replace cont_pre = okpre & cont_pre[_n+1] if event_week==`wk'
        }

        * POST contiguous block starting at 0
        gen byte okpost = pass & (event_week>=0)
        gen byte cont_post = .
        replace cont_post = okpost if event_week==0

        forvalues k=1/$MAXW {
            quietly replace cont_post = okpost & cont_post[_n-1] if event_week==`k'
        }

        quietly summarize event_week if cont_pre==1, meanonly
        scalar PRE_sug = r(min)

        quietly summarize event_week if cont_post==1, meanonly
        scalar POST_sug = r(max)

        di as text "[DEBUG 2A] Suggested window: PRE=" PRE_sug " POST=" POST_sug

    restore

    * Set globals used downstream
    global WPRE  = -PRE_sug
    global WPOST =  POST_sug
}
else {
    global WPRE  $WPRE_MAN
    global WPOST $WPOST_MAN
}

display as text "[DEBUG 2A] Using window [-$WPRE, $WPOST] (AUTO_WIN_2A=$AUTO_WIN_2A, TAU=$TAU, CORE=$CORE0..$CORE1)"

* shifted event week for plotting/factor vars
gen int event_week_s = event_week + 100 if sample_2A==1 & inrange(event_week, -$WPRE, $WPOST)

* label event time for plotting
cap label drop ew_lbl
local lo = 100 - $WPRE
local hi = 100 + $WPOST
forval i = `lo'/`hi' {
    local val = `i' - 100
    label define ew_lbl `i' "`val'", add
}
label values event_week_s ew_lbl
global EW_BASE 99   // week -1

* --- 2A DiD ---
reghdfe $Y ///
    i.treated_BC##i.post ///
    $X ///
    if (sample_2A==1) & inrange(event_week, -$WPRE, $WPOST) & !missing($Y) & ($SAMPLE_OK), ///
    absorb($ABS) ///
    vce(cluster contract_id)
estimates store did_2A

* --- 2A Event study ---
reghdfe $Y ///
    ib($EW_BASE).event_week_s##i.treated_BC ///
    $X ///
    if (sample_2A==1) & !missing(event_week_s) & !missing($Y) & ($SAMPLE_OK), ///
    absorb($ABS) ///
    vce(cluster contract_id)
estimates store es_2A

* --- 2A Plot ---
local clabel ""
forval i = `lo'/`hi' {
    local val = `i' - 100
    local clabel `clabel' `i'.event_week_s#1.treated_BC = "`val'" 1.treated_BC#`i'.event_week_s = "`val'"
}
local line_pos = $WPRE + 0.5

coefplot, ///
    keep(*.event_week_s#1.treated_BC 1.treated_BC#*.event_week_s) ///
    omitted baselevels vertical ///
    coeflabel(`clabel') ///
    yline(0, lpattern(dash)) ///
    xline(`line_pos', lwidth(medium)) ///
    title("Event study: $Y ($TY)", size(medium)) ///
    subtitle("BC vs $CONTROL | Window [-$WPRE, $WPOST] | TAU=$TAU | Floors=$USE_FLOORS", size(small)) ///
    ytitle("Effect relative to week -1") ///
    xtitle("Weeks relative to June 1 cutoff") ///
    graphregion(color(white)) xsize(12) ysize(7) ///
    xlabel(, labsize(small)) grid(none)

graph export "event_2A_${Y}_${TY}_${CONTROL}.png", replace


/**********************************************************************
  2B WINDOW SELECTION (TAU RULE; INLINE, ROBUST)
  Fix: Use group_s (short string) to avoid strL failures in collapse/reshape.
  Requires support in BOTH treated and placebo groups.
**********************************************************************/

* Short safe group label for collapse/reshape
cap drop group_s
gen str8 group_s = group_str if bc_placebo_sample==1

if $AUTO_WIN_2B==1 {

    preserve
        keep if bc_placebo_sample==1
        keep if inrange(event_week_y, -$MAXW, $MAXW)

        collapse (count) N=event_week_y, by(event_week_y group_s)
        reshape wide N, i(event_week_y) j(group_s) string

        capture confirm variable Ntreated
        if _rc {
            di as error "Expected Ntreated after reshape, not found. Check group_s."
            list in 1/15
            restore
            exit 198
        }
        capture confirm variable Nplacebo
        if _rc {
            di as error "Expected Nplacebo after reshape, not found. Check group_s."
            list in 1/15
            restore
            exit 198
        }

        quietly summarize Ntreated if inrange(event_week_y,$CORE0,$CORE1)
        scalar coreT = r(mean)

        quietly summarize Nplacebo if inrange(event_week_y,$CORE0,$CORE1)
        scalar coreP = r(mean)

        scalar thrT = $TAU * coreT
        scalar thrP = $TAU * coreP

        di as text "[DEBUG 2B] Core mean counts: treated=" coreT " placebo=" coreP
        di as text "[DEBUG 2B] Thresholds: treated>=" thrT " placebo>=" thrP

        gen byte pass = (Ntreated >= thrT & Nplacebo >= thrP)

        * PRE contiguous block ending at -1
        sort event_week_y
        gen byte okpre = pass & (event_week_y<=-1)
        gen byte cont_pre = .
        replace cont_pre = okpre if event_week_y==-1

        forvalues k=2/$MAXW {
            local wk = -`k'
            quietly replace cont_pre = okpre & cont_pre[_n+1] if event_week_y==`wk'
        }

        * POST contiguous block starting at 0
        gen byte okpost = pass & (event_week_y>=0)
        gen byte cont_post = .
        replace cont_post = okpost if event_week_y==0

        forvalues k=1/$MAXW {
            quietly replace cont_post = okpost & cont_post[_n-1] if event_week_y==`k'
        }

        quietly summarize event_week_y if cont_pre==1, meanonly
        scalar PRE_sug_BC = r(min)

        quietly summarize event_week_y if cont_post==1, meanonly
        scalar POST_sug_BC = r(max)

        di as text "[DEBUG 2B] Suggested window: PRE=" PRE_sug_BC " POST=" POST_sug_BC

    restore

    global PRE_MIN_BC  = PRE_sug_BC
    global POST_MAX_BC = POST_sug_BC
}
else {
    global PRE_MIN_BC  $PRE_MIN_BC_MAN
    global POST_MAX_BC $POST_MAX_BC_MAN
}

display as text "[DEBUG 2B] Using window [$PRE_MIN_BC, $POST_MAX_BC] (AUTO_WIN_2B=$AUTO_WIN_2B, TAU=$TAU, CORE=$CORE0..$CORE1)"

* shifted week var for factor vars / plotting (avoid negatives)
gen int event_week_sy = event_week_y + 100 if bc_placebo_sample==1 & inrange(event_week_y, $PRE_MIN_BC, $POST_MAX_BC)

* labels for plotting
cap label drop ew_lbl
local lo2 = 100 + $PRE_MIN_BC
local hi2 = 100 + $POST_MAX_BC
forval i = `lo2'/`hi2' {
    local val = `i' - 100
    label define ew_lbl `i' "`val'", add
}
label values event_week_sy ew_lbl
global EW_BASE_BC 99

* --- 2B DiD ---
reghdfe $Y ///
    i.treated_year##i.post_y ///
    $X ///
    if (bc_placebo_sample==1) & inrange(event_week_y, $PRE_MIN_BC, $POST_MAX_BC) & !missing($Y) & ($SAMPLE_OK), ///
    absorb($ABS) ///
    vce(cluster contract_id)
estimates store did_2B

* --- 2B Event study ---
reghdfe $Y ///
    ib($EW_BASE_BC).event_week_sy##i.treated_year ///
    $X ///
    if (bc_placebo_sample==1) & !missing(event_week_sy) & !missing($Y) & ($SAMPLE_OK), ///
    absorb($ABS) ///
    vce(cluster contract_id)
estimates store es_2B

* --- 2B Plot ---
local clabel2 ""
forval i = `lo2'/`hi2' {
    local val = `i' - 100
    local clabel2 `clabel2' `i'.event_week_sy#1.treated_year = "`val'" 1.treated_year#`i'.event_week_sy = "`val'"
}
local wpre_bc = -$PRE_MIN_BC
local line_pos2 = `wpre_bc' + 0.5

coefplot, ///
    keep(*.event_week_sy#1.treated_year 1.treated_year#*.event_week_sy) ///
    omitted baselevels vertical ///
    coeflabel(`clabel2') ///
    yline(0, lpattern(dash)) ///
    xline(`line_pos2', lwidth(medium)) ///
    title("Event study: $Y (BC only)", size(medium)) ///
    subtitle("Treated $TY_BC vs placebo ($BCPLACEBO) | Window [$PRE_MIN_BC, $POST_MAX_BC] | TAU=$TAU | Floors=$USE_FLOORS", size(small)) ///
    ytitle("Effect relative to week -1") ///
    xtitle("Weeks relative to June 1 cutoff") ///
    graphregion(color(white)) xsize(12) ysize(7) ///
    xlabel(, labsize(small)) grid(none)

graph export "event_2B_${Y}_${TY_BC}.png", replace


















/*
/**********************************************************************
  MASTER DO-FILE (Everything before the "additional analysis" block)
  This sets up:
  - import + cleaning
  - date conversion
  - destring numeric vars
  - core outcomes (ln_incentive, lpiece, lprod)
  - identifiers (worker_id, contract_id, wc_id)
  - day-level hours / multi-contract indicators (tot_hrs_day, n_contracts_day, etc.)
  - precipitation zero + positive components (precip0, precip_pos) used later
**********************************************************************/

clear all
set more off


/**********************************************************************
  0) LOAD DATA
**********************************************************************/

import delimited using "/Users/danaandersen/Downloads/planter_mw.csv", clear varnames(1)


/**********************************************************************
  1) SANITY CHECKS: REQUIRED RAW VARIABLES
**********************************************************************/

describe
confirm variable date year province contract name piece prod top_up tot_hrs
confirm variable temp total_precip wind_spd experience


/**********************************************************************
  2) DATE: ENSURE STATA DAILY DATE
**********************************************************************/

capture confirm numeric variable date
if _rc {
    gen double date_num = daily(date, "YMD")
    format date_num %td
    drop date
    rename date_num date
}
format date %td


/**********************************************************************
  3) NUMERIC STANDARDIZATION (SAFE DESTRING)
**********************************************************************/

foreach var in temp total_precip wind_spd experience tot_hrs piece prod top_up {
    capture confirm numeric variable `var'
    if _rc destring `var', replace ignore("NA" "." "" "null" "NULL")
}


/**********************************************************************
  4) CORE OUTCOMES
**********************************************************************/

cap drop ln_incentive lpiece lprod

* 1) Incentive-based hourly earnings (log of piece*prod)
gen double ln_incentive = ln(piece * prod) if piece>0 & prod>0

* 2) Productivity (log). Assumes prod is already trees/hour.
gen double lprod  = ln(prod)  if prod>0


/**********************************************************************
  5) STATIC IDENTIFIERS
**********************************************************************/

cap drop worker_id contract_id wc_id
encode name, gen(worker_id)
encode contract, gen(contract_id)
egen wc_id = group(worker_id contract_id)


/**********************************************************************
  6) DAY-LEVEL AGGREGATES / FLAGS (HANDLES MULTI-CONTRACT DAYS)
     tot_hrs is worker×contract-day hours in your data.
     tot_hrs_day sums across all contracts for the worker on that day.
**********************************************************************/

cap drop tot_hrs_day n_contracts_day multi_contract_day long_day extreme_day

bysort worker_id date: egen double tot_hrs_day = total(tot_hrs)
bysort worker_id date: egen int    n_contracts_day = count(contract_id)

gen byte multi_contract_day = (n_contracts_day > 1)
gen byte long_day           = (tot_hrs_day > 14)
gen byte extreme_day        = (tot_hrs_day > 18)

* quick checks
summ tot_hrs tot_hrs_day, detail
tab n_contracts_day
tab multi_contract_day


/**********************************************************************
  7) WEATHER TRANSFORMS USED LATER (NONLINEAR PRECIP)
**********************************************************************/

cap drop precip0 precip_pos

gen byte precip0 = (total_precip == 0) if !missing(total_precip)
gen double precip_pos = cond(total_precip > 0 & !missing(total_precip), total_precip, 0)

* (optional; if you want quadratic explicitly created)
cap drop temp2 wind2 precip2
gen double temp2   = temp^2       if !missing(temp)
gen double wind2   = wind_spd^2   if !missing(wind_spd)
gen double precip2 = precip_pos^2 if !missing(precip_pos)


/**********************************************************************
  8) OPTIONAL: QUICK CHECKS (SAFE TO DELETE)
**********************************************************************/

summ ln_incentive lpiece lprod
tab province
tab year



/**********************************************************************
  PART B: ADDITIONAL ANALYSIS CONTROL PANEL (RERUN FROM HERE)
  - BC placebo control (within BC)
  - AB control (BC vs AB)
  - Triple diff (optional)
  - Pre-treatment worker ability (AKM-style, pre-only)
  - Heterogeneous MW effects
  - Selection tests
  - Sorting tests
**********************************************************************/

/*** 0) FORCE CLEAR (safe rerun) ***/
foreach v in tdate_y post_y event_day_y event_week_y ///
            bc_placebo_sample ab_control_sample triplediff_sample treated_year ///
            ln_prod abil_sample fe_worker fe_contract ability_raw ability_z ability_hi ///
            treated_BC_y treated_BC_td ///
            any_post any_pre ln_piece {
    cap drop `v'
}
cap estimates clear


/**********************************************************************
  1) USER CHOICES (EDIT THESE)
**********************************************************************/

* Treated MW year in BC
global TY 2019

* BC placebo years with no MW change (numlist)
global BCPLACEBO "2016 2017"

* Province control choice for cross-province DiD
global PROVCTRL "AB"

* Window around June 1 (set based on support)
global PRE_MIN  -3
global POST_MAX  8


/**********************************************************************
  2) YEAR-SPECIFIC JUNE 1 CUTOFF + EVENT TIME (WORKS FOR ALL YEARS)
**********************************************************************/

* June 1 of each observation's year
gen double tdate_y = mdy(6,1,year)
format tdate_y %td

gen byte post_y = (date >= tdate_y) if !missing(date) & !missing(tdate_y)

gen int event_day_y  = date - tdate_y if !missing(date) & !missing(tdate_y)
gen int event_week_y = floor(event_day_y/7) if !missing(event_day_y)


/**********************************************************************
  3) SAMPLES (BC PLACEBO / AB CONTROL / TRIPLE DIFF)
**********************************************************************/

* helper: placebo-year indicator (avoids inlist() problems)
gen byte placebo_year = 0
foreach y of numlist $BCPLACEBO {
    replace placebo_year = 1 if year==`y'
}

* A) BC-only: treated year vs placebo years
gen byte bc_placebo_sample = (province=="BC") & (year==$TY | placebo_year==1)

* B) Treated year only: BC vs AB/ON
gen byte ab_control_sample = (year==$TY) & (province=="BC" | province=="$PROVCTRL")

* C) Triple diff sample: BC and control provinces in treated + placebo years
gen byte triplediff_sample = (province=="BC" | province=="$PROVCTRL") & (year==$TY | placebo_year==1)

* Treated-year indicator
gen byte treated_year = (year==$TY)

* quick checks
tab year if bc_placebo_sample
tab province if ab_control_sample
tab year province if triplediff_sample


/**********************************************************************
  4) PRE-TREATMENT WORKER ABILITY (AKM STYLE ON PRE DATA ONLY)
     - Ability = worker FE from pre-period ln(prod), net of contract FE
     - Uses BC treated year + BC placebo years, PRE-only window
**********************************************************************/

* log productivity (trees/hour) – assumes prod is already trees/hour
gen double ln_prod = ln(prod) if prod>0

gen byte abil_sample = bc_placebo_sample ///
    & post_y==0 ///
    & inrange(event_week_y, $PRE_MIN, -1) ///
    & !missing(ln_prod)

count if abil_sample
summ ln_prod if abil_sample

* Try to save FE with explicit names if your reghdfe supports it.
* If it errors, fall back to savefe and rename.

cap noi reghdfe ln_prod ///
    c.experience##c.experience ///
    c.temp##c.temp c.wind_spd##c.wind_spd ///
    precip0 c.precip_pos##c.precip_pos ///
    if abil_sample, ///
    absorb(worker_id=fe_worker contract_id=fe_contract) ///
    vce(cluster contract_id)

if _rc {
    * fallback: use savefe and then map __hdfe* to worker/contract by absorb order
    cap drop fe_worker fe_contract __hdfe*
    reghdfe ln_prod ///
        c.experience##c.experience ///
        c.temp##c.temp c.wind_spd##c.wind_spd ///
        precip0 c.precip_pos##c.precip_pos ///
        if abil_sample, ///
        absorb(worker_id contract_id, savefe) ///
        vce(cluster contract_id)

    * With absorb(worker_id contract_id), __hdfe1__ is worker FE, __hdfe2__ is contract FE
    gen double fe_worker   = __hdfe1__
    gen double fe_contract = __hdfe2__
}

* ability = worker FE from PRE data
gen double ability_raw = fe_worker

* standardize for interactions
summ ability_raw if !missing(ability_raw)
gen double ability_z = (ability_raw - r(mean)) / r(sd) if !missing(ability_raw)

* optional: low/high split
xtile ability_hi = ability_raw if !missing(ability_raw), nq(2)
replace ability_hi = ability_hi - 1   // 0=low, 1=high

summ ability_raw ability_z
tab ability_hi


/**********************************************************************
  5) HETEROGENEOUS MW EFFECTS — AB AS CONTROL (TREATED YEAR ONLY)
     Model: i.BC ## i.post ## ability_z
**********************************************************************/

gen byte treated_BC_y = (province=="BC") if ab_control_sample

reghdfe ln_incentive ///
    i.treated_BC_y##i.post_y##c.ability_z ///
    tot_hrs_day multi_contract_day ///
    c.experience##c.experience ///
    c.temp##c.temp c.wind_spd##c.wind_spd ///
    precip0 c.precip_pos##c.precip_pos ///
    if ab_control_sample ///
       & inrange(event_week_y, $PRE_MIN, $POST_MAX) ///
       & !missing(ln_incentive) & !missing(ability_z), ///
    absorb(contract_id) ///
    vce(cluster contract_id)


/**********************************************************************
  6) HETEROGENEOUS MW EFFECTS — BC PLACEBO YEARS AS CONTROL (WITHIN BC)
     Model: i.treated_year ## i.post ## ability_z   (BC only)
**********************************************************************/

reghdfe ln_incentive ///
    i.treated_year##i.post_y##c.ability_z ///
    tot_hrs_day multi_contract_day ///
    c.experience##c.experience ///
    c.temp##c.temp c.wind_spd##c.wind_spd ///
    precip0 c.precip_pos##c.precip_pos ///
    if bc_placebo_sample ///
       & inrange(event_week_y, $PRE_MIN, $POST_MAX) ///
       & !missing(ln_incentive) & !missing(ability_z), ///
    absorb(contract_id) ///
    vce(cluster contract_id)


/**********************************************************************
  7) OPTIONAL: TRIPLE DIFF (BC vs AB) × (treated year vs placebo) × (post)
     Model: i.BC ## i.treated_year ## i.post ## ability_z
**********************************************************************/

gen byte treated_BC_td = (province=="BC") if triplediff_sample

reghdfe ln_incentive ///
    i.treated_BC_td##i.treated_year##i.post_y##c.ability_z ///
    tot_hrs_day multi_contract_day ///
    c.experience##c.experience ///
    c.temp##c.temp c.wind_spd##c.wind_spd ///
    precip0 c.precip_pos##c.precip_pos ///
    if triplediff_sample ///
       & inrange(event_week_y, $PRE_MIN, $POST_MAX) ///
       & !missing(ln_incentive) & !missing(ability_z), ///
    absorb(contract_id) ///
    vce(cluster contract_id)


/**********************************************************************
  8) SELECTION TESTS
     A) Within BC: treated year vs placebo years
     B) BC vs AB in treated year
     Outcome: any_post (worked at least once post June 1)
**********************************************************************/

preserve
keep if bc_placebo_sample
keep worker_id year post_y ability_z

bys worker_id year: egen any_post = max(post_y==1)
bys worker_id year: egen any_pre  = max(post_y==0)
bys worker_id year: keep if _n==1
keep if any_pre==1

gen byte treated_year = (year==$TY)

reg any_post c.ability_z##i.treated_year, vce(cluster worker_id)
restore


preserve
keep if ab_control_sample
keep worker_id province post_y ability_z

bys worker_id province: egen any_post = max(post_y==1)
bys worker_id province: egen any_pre  = max(post_y==0)
bys worker_id province: keep if _n==1
keep if any_pre==1

gen byte treated_BC = (province=="BC")

reg any_post c.ability_z##i.treated_BC, vce(cluster worker_id)
restore


/**********************************************************************
  9) SORTING TESTS (ABILITY × PIECE RATE MATCHING)
     A) Within BC: treated year vs placebo years
     B) BC vs AB in treated year
**********************************************************************/

gen double ln_piece = ln(piece) if piece>0

* A) BC placebo design
reghdfe ln_piece ///
    i.treated_year##i.post_y##c.ability_z ///
    c.experience##c.experience ///
    c.temp##c.temp c.wind_spd##c.wind_spd ///
    precip0 c.precip_pos##c.precip_pos ///
    if bc_placebo_sample ///
       & inrange(event_week_y, $PRE_MIN, $POST_MAX) ///
       & !missing(ln_piece) & !missing(ability_z), ///
    absorb() ///
    vce(cluster contract_id)

* B) AB control design (treated year only)
reghdfe ln_piece ///
    i.treated_BC_y##i.post_y##c.ability_z ///
    c.experience##c.experience ///
    c.temp##c.temp c.wind_spd##c.wind_spd ///
    precip0 c.precip_pos##c.precip_pos ///
    if ab_control_sample ///
       & inrange(event_week_y, $PRE_MIN, $POST_MAX) ///
       & !missing(ln_piece) & !missing(ability_z), ///
    absorb() ///
    vce(cluster contract_id)





	
	
	
	
	


/*
/**********************************************************************
  PROJECT: Minimum Wage and Productivity
  PURPOSE: Persistent Section-by-Section Analysis (No Locals)
**********************************************************************/

clear all
set more off

/**********************************************************************
  0) LOAD DATA
**********************************************************************/
import delimited using "/Users/danaandersen/Downloads/planter_mw.csv", clear varnames(1)

/**********************************************************************
  1) CLEAN AND STANDARDIZE CORE VARIABLES
**********************************************************************/
capture confirm numeric variable date
if _rc {
    gen double date_num = daily(date, "YMD")
    format date_num %td
    drop date
    rename date_num date
}

foreach var in temp total_precip wind_spd experience tot_hrs piece prod top_up {
    capture confirm numeric variable `var'
    if _rc destring `var', replace ignore("NA" "." "" "null" "NULL")
}

gen double ln_incentive = ln(piece * prod) if piece > 0 & prod > 0

/**********************************************************************
  2) TOP-UP STATUS AND IDENTIFIERS
**********************************************************************/
bysort name: egen byte ever_topup = max(top_up)
gen byte never_topup = (ever_topup == 0)

encode name, gen(worker_id)
encode contract, gen(contract_id)
egen wc_id = group(worker_id contract_id)

/**********************************************************************
  6) WORKER-DAY AGGREGATION
**********************************************************************/
bysort worker_id date: egen tot_hrs_day = total(tot_hrs)
bysort worker_id date: egen n_contracts_day = count(contract_id)
gen byte multi_contract_day = (n_contracts_day > 1)

/**********************************************************************
  8) WEATHER CONTROLS (Persistent Globals)
**********************************************************************/
gen byte precip0 = (total_precip == 0) if !missing(total_precip)
gen double precip_pos = cond(total_precip > 0 & !missing(total_precip), total_precip, 0)

* We use a global so it stays in memory across different runs
global WEATHER "c.temp##c.temp c.wind_spd##c.wind_spd precip0 c.precip_pos##c.precip_pos"

/**********************************************************************
  9) ANALYSIS PREP (Set Year here manually if running section-by-section)
**********************************************************************/
global TY 2018
global CONTROL "AB"
global NEVERONLY 0
global WPRE  3   
global WPOST 8

global TDATE = td(01jun$TY)

* Create persistent filter
macro drop IFBASE
global IFBASE (year == $TY)
if "$CONTROL" == "AB"   global IFBASE `"$IFBASE & (province == "BC" | province == "AB")"'
if "$CONTROL" == "ON"   global IFBASE `"$IFBASE & (province == "BC" | province == "ON")"'
if "$CONTROL" == "ABON" global IFBASE `"$IFBASE & (province == "BC" | province == "AB" | province == "ON")"'
if $NEVERONLY == 1      global IFBASE `"$IFBASE & never_topup == 1"'

/**********************************************************************
  10) TREATMENT AND SHIFTED TIME (PERSISTENT VARIABLES)
**********************************************************************/
cap drop treated_BC post event_week event_week_s
gen byte treated_BC = (province == "BC")
gen byte post       = (date >= $TDATE)
gen int event_week  = floor((date - $TDATE)/7)

* Persistent Shift Variable (No negatives)
gen int event_week_s = event_week + 100
* Week -1 becomes 99
global EW_BASE 99

/**********************************************************************
  11) REGRESSIONS (Can be run independently)
**********************************************************************/

* BASELINE (No Weather)
reghdfe ln_incentive treated_BC##post tot_hrs tot_hrs_day multi_contract_day ///
    c.experience##c.experience if ($IFBASE), ///
    absorb(worker_id contract_id) vce(cluster contract_id)

* EVENT STUDY (With Weather)
reghdfe ln_incentive ib($EW_BASE).event_week_s##i.treated_BC ///
    tot_hrs tot_hrs_day multi_contract_day ///
    c.experience##c.experience $WEATHER ///
    if ($IFBASE) & inrange(event_week, -$WPRE, $WPOST), ///
    absorb(worker_id contract_id) vce(cluster contract_id)

/**********************************************************************
  13) PRE-TREND TEST (Fixed for Sectional Runs)
**********************************************************************/
* Since we can't use a loop to build a local list across sections easily, 
* we call the specific shifted weeks for our 3-week pre-window (-3, -2)
test 1.treated_BC#97.event_week_s 1.treated_BC#98.event_week_s

/**********************************************************************
  14) COEFPLOT (Manual Labeling for Sectional Runs)
**********************************************************************/
coefplot, keep(1.treated_BC#*.event_week_s) ///
    omitted baselevels vertical yline(0) ///
    rename(97.event_week_s="-3" 98.event_week_s="-2" 99.event_week_s="-1" ///
           100.event_week_s="0" 101.event_week_s="1" 102.event_week_s="2" ///
           103.event_week_s="3" 104.event_week_s="4" 105.event_week_s="5" ///
           106.event_week_s="6" 107.event_week_s="7" 108.event_week_s="8") ///
    title("Event Study: BC $TY") ///
    name(graph_$TY, replace)
	
	
	
	


/**********************************************************************
  PROJECT: Minimum Wage and Productivity
  PURPOSE: Load + clean data, build analysis sample, construct controls
           (including quadratic experience + quadratic weather with
           zero-precip dummy), and prepare for regressions.

  NOTE: This file stops RIGHT BEFORE the regression commands.
**********************************************************************/

clear all
set more off


/**********************************************************************
  0) LOAD DATA AND BASIC CHECKS
**********************************************************************/

import delimited using "/Users/danaandersen/Downloads/planter_mw.csv", clear varnames(1)

describe

* Core variables needed up to regressions
confirm variable date year province contract name piece prod top_up tot_hrs ///
               temp total_precip wind_spd experience


/**********************************************************************
  1) CLEAN AND STANDARDIZE CORE VARIABLES
**********************************************************************/

* --- date: ensure numeric Stata daily date ---
capture confirm numeric variable date
if _rc {
    capture drop date_num
    gen double date_num = daily(date, "YMD")
    format date_num %td
    drop date
    rename date_num date
}
else {
    format date %td
}

* --- weather controls: string -> numeric (safe destring) ---
capture confirm numeric variable temp
if _rc destring temp, replace ignore("NA" "." "" "null" "NULL")

capture confirm numeric variable total_precip
if _rc destring total_precip, replace ignore("NA" "." "" "null" "NULL")

capture confirm numeric variable wind_spd
if _rc destring wind_spd, replace ignore("NA" "." "" "null" "NULL")

* --- experience: ensure numeric ---
capture confirm numeric variable experience
if _rc destring experience, replace ignore("NA" "." "" "null" "NULL")

* --- tot_hrs: ensure numeric (hours on worker–contract–day) ---
capture confirm numeric variable tot_hrs
if _rc destring tot_hrs, replace ignore("NA" "." "" "null" "NULL")

* --- outcome: log incentive earnings (piece-rate based) ---
capture drop ln_incentive
gen double ln_incentive = ln(piece*prod) if piece>0 & prod>0

summ ln_incentive piece prod tot_hrs


/**********************************************************************
  2) TOP-UP STATUS (WORKER-LEVEL)
**********************************************************************/

capture drop ever_topup never_topup
bysort name: egen byte ever_topup = max(top_up)
gen byte never_topup = (ever_topup==0)

tab never_topup


/**********************************************************************
  PURPOSE: Analysis for 2018 & 2019 MW Changes
**********************************************************************/

* (Sections 1-2: Load and Clean as per your original script)
* [Assume data is loaded, Piece/Prod cleaned, and date_num created]

/**********************************************************************
  3) GLOBAL SETTINGS & WINDOW
**********************************************************************/
global CONTROL "AB"
global NEVERONLY 0
global WPRE  3   // Updated window
global WPOST 8

* Define weather robustness local (Comment this out to exclude from regressions)
local weather "c.temp##c.temp c.wind_spd##c.wind_spd precip0 c.precip_pos##c.precip_pos"

/**********************************************************************
  LOOP THROUGH TREATMENT YEARS
**********************************************************************/
foreach ty in 2018 2019 {
    
    display "************************** RUNNING YEAR: `ty' **************************"
    global TY `ty'
    global TDATE = td(01jun$TY)

    * --- 4) Setup IFBASE ---
    macro drop IFBASE
    global IFBASE (year == $TY)
    if "$CONTROL" == "AB"   global IFBASE `"$IFBASE & (province == "BC" | province == "AB")"'
    if "$CONTROL" == "ON"   global IFBASE `"$IFBASE & (province == "BC" | province == "ON")"'
    if "$CONTROL" == "ABON" global IFBASE `"$IFBASE & (province == "BC" | province == "AB" | province == "ON")"'
    if $NEVERONLY == 1      global IFBASE `"$IFBASE & never_topup == 1"'

    * --- 7) Treatment Variables ---
    cap drop treated_BC post event_week event_week_s
    gen byte treated_BC = (province == "BC")
    gen byte post       = (date >= $TDATE)
    gen int event_week  = floor((date - $TDATE)/7)
    
    * --- 11) Shift Event Time (Fixes "No Negatives" error) ---
    * Shift by 100 to stay safely positive
    gen int event_week_s = event_week + 100
    local ew_base = -1 + 100

    /**********************************************************************
      10) BASELINE DiD (Without Weather)
    **********************************************************************/
    reghdfe ln_incentive treated_BC##post tot_hrs tot_hrs_day multi_contract_day ///
        c.experience##c.experience if ($IFBASE), ///
        absorb(worker_id contract_id) vce(cluster contract_id)
    
    /**********************************************************************
      11) EVENT STUDY
    **********************************************************************/
    reghdfe ln_incentive ib(`ew_base').event_week_s##i.treated_BC ///
        tot_hrs tot_hrs_day multi_contract_day ///
        c.experience##c.experience `weather' ///  // Weather included here
        if ($IFBASE) & inrange(event_week, -$WPRE, $WPOST), ///
        absorb(worker_id contract_id) vce(cluster contract_id)

    /**********************************************************************
      13) DIAGNOSTICS: PRE-TREND TEST
    **********************************************************************/
    * Collect all interaction terms where original week < 0 (excluding the base)
    local pre_trends ""
    forval w = -$WPRE / -2 {
        local s_w = `w' + 100
        local pre_trends "`pre_trends' 1.treated_BC#`s_w'.event_week_s"
    }
    test `pre_trends'

    /**********************************************************************
      14) AUTOMATED COEFPLOT
    **********************************************************************/
    * This builds a mapping string to change "103" back to "3" on the axis
    local lab_map ""
    forval w = -$WPRE / $WPOST {
        local s_w = `w' + 100
        local lab_map `"`lab_map' `s_w'.event_week_s = "`w'""'
    }

    coefplot, keep(1.treated_BC#*.event_week_s) ///
        omitted baselevels ///
        vertical yline(0) ///
        rename(`lab_map') ///  // This applies the original week labels
        title("Event Study: BC $TY Change") ///
        xtitle("Weeks Relative to June 1st") ///
        name(graph_$TY, replace)

} // End of Year Loop






/*

/**********************************************************************
  3) GLOBAL ANALYSIS OPTIONS
**********************************************************************/

* Treated year (2018 or 2019)
global TY 2018

* Control group:
* "AB"   = Alberta
* "ON"   = Ontario
* "ABON" = Alberta + Ontario
global CONTROL "AB"

* Restrict to never-top-up workers? (1=yes, 0=no)
global NEVERONLY 0

* Event-study window (weeks)
global WPRE  8
global WPOST 12


/**********************************************************************
  4) ANALYSIS SAMPLE FILTER (IFBASE) — uses your working quoting style
**********************************************************************/

macro drop IFBASE
global IFBASE year==$TY

if "$CONTROL"=="AB"   ///
    global IFBASE `"$IFBASE & (province=="BC" | province=="AB")"'

if "$CONTROL"=="ON"   ///
    global IFBASE `"$IFBASE & (province=="BC" | province=="ON")"'

if "$CONTROL"=="ABON" ///
    global IFBASE `"$IFBASE & (province=="BC" | province=="AB" | province=="ON")"'

if $NEVERONLY==1 ///
    global IFBASE `"$IFBASE & never_topup==1"'

display `"[DEBUG] IFBASE: $IFBASE"'
count if $IFBASE
tab province if $IFBASE


/**********************************************************************
  5) IDENTIFIERS (WORKER, CONTRACT, WORKER×CONTRACT)
**********************************************************************/

* Worker ID
capture confirm numeric variable worker_id
if _rc {
    capture drop worker_id
    encode name, gen(worker_id)
}

* Contract ID
capture confirm numeric variable contract_id
if _rc {
    capture drop contract_id
    encode contract, gen(contract_id)
}

* Worker × Contract ID (for robustness FE)
capture confirm variable wc_id
if _rc {
    capture drop wc_id
    egen wc_id = group(worker_id contract_id)
}


/**********************************************************************
  6) WORKER-DAY AGGREGATION AND FLAGS
     tot_hrs = hours on worker–contract–day
     tot_hrs_day = total hours across all contracts for worker-day
**********************************************************************/

capture drop tot_hrs_day n_contracts_day ///
             multi_contract_day long_day extreme_day

bysort worker_id date: egen tot_hrs_day = total(tot_hrs)
bysort worker_id date: egen n_contracts_day = count(contract_id)

gen byte multi_contract_day = (n_contracts_day > 1)
gen byte long_day           = (tot_hrs_day > 14)
gen byte extreme_day        = (tot_hrs_day > 18)

summ tot_hrs tot_hrs_day, detail
tab n_contracts_day
tab multi_contract_day
tab extreme_day


/**********************************************************************
  7) TREATMENT TIMING AND EVENT TIME
**********************************************************************/

macro drop TDATE
global TDATE = td(01jun$TY)

capture drop treated_BC post event_day event_week
gen byte treated_BC = (province=="BC")

gen byte post = .
replace post = (date >= $TDATE) if year==$TY

gen int event_day  = .
replace event_day = date - $TDATE if year==$TY

gen int event_week = .
replace event_week = floor(event_day/7) if year==$TY

tab treated_BC if $IFBASE
tab post if $IFBASE
summ event_week if $IFBASE


/**********************************************************************
  8) FLEXIBLE CONTROLS (CONSISTENT QUADRATICS + PRECIP ZERO DUMMY)
     - Experience: quadratic
     - Temperature: quadratic
     - Wind: quadratic
     - Precipitation: zero dummy + quadratic in level (for precip>0)
       (precip level enters with precip and precip^2, with precip0 dummy)
**********************************************************************/

* Experience quadratic
capture drop experience2
gen double experience2 = experience^2 if !missing(experience)

* Temperature quadratic
capture drop temp2
gen double temp2 = temp^2 if !missing(temp)

* Wind quadratic
capture drop wind2
gen double wind2 = wind_spd^2 if !missing(wind_spd)

* Precip: zero dummy + quadratic in level
capture drop precip0 precip_pos precip2
gen byte precip0 = (total_precip==0) if !missing(total_precip)

* For quadratic intensity effects, define precip_pos = precip for precip>0, else 0
gen double precip_pos = 0
replace precip_pos = total_precip if total_precip>0 & !missing(total_precip)

gen double precip2 = precip_pos^2

* Quick checks
summ experience experience2 temp temp2 wind_spd wind2 total_precip precip0 precip_pos precip2
tab precip0 if $IFBASE



/**********************************************************************
  9) INSTALL REQUIRED PACKAGES (ONCE PER MACHINE)
**********************************************************************/

capture which reghdfe
if _rc {
    ssc install ftools, replace
    ssc install reghdfe, replace
}

capture which coefplot
if _rc {
    ssc install coefplot, replace
}


/**********************************************************************
  10) BASELINE DiD ESTIMATION
      - Outcome: ln(piece * prod)
      - FE: worker + contract
      - SEs: clustered at contract
**********************************************************************/

reghdfe ln_incentive ///
    treated_BC##post ///
    tot_hrs tot_hrs_day multi_contract_day ///
    experience experience2 ///
    temp temp2 ///
    wind_spd wind2 ///
    precip0 precip_pos precip2 ///
    if $IFBASE, ///
    absorb(worker_id contract_id) ///
    vce(cluster contract_id)


/**********************************************************************
  11) EVENT STUDY — BASELINE FE (worker + contract)
      FIX FOR NEGATIVE EVENT TIME:
      - Shift event_week so it is non-negative
      - Omit baseline corresponding to original week -1 manually
**********************************************************************/

* --- Shift event time ---
capture drop event_week_s
quietly summarize event_week if $IFBASE & !missing(event_week)
local ew_shift = -r(min)

gen int event_week_s = event_week + `ew_shift'

* Identify shifted baseline corresponding to event_week == -1
local ew_base = -1 + `ew_shift'

* Sanity check (optional)
* tab event_week event_week_s if $IFBASE & inrange(event_week,-$WPRE,$WPOST)

reghdfe ln_incentive ///
    i.event_week_s##i.treated_BC ///
    tot_hrs tot_hrs_day multi_contract_day ///
    experience experience2 ///
    temp temp2 ///
    wind_spd wind2 ///
    precip0 precip_pos precip2 ///
    if $IFBASE ///
       & inrange(event_week, -$WPRE, $WPOST) ///
       & event_week_s != `ew_base', ///
    absorb(worker_id contract_id) ///
    vce(cluster contract_id)


/**********************************************************************
  12) EVENT STUDY — ROBUSTNESS: WORKER × CONTRACT FE
**********************************************************************/

reghdfe ln_incentive ///
    i.event_week_s##i.treated_BC ///
    tot_hrs tot_hrs_day multi_contract_day ///
    experience experience2 ///
    temp temp2 ///
    wind_spd wind2 ///
    precip0 precip_pos precip2 ///
    if $IFBASE ///
       & inrange(event_week, -$WPRE, $WPOST) ///
       & event_week_s != `ew_base', ///
    absorb(wc_id) ///
    vce(cluster contract_id)


/**********************************************************************
  13) DIAGNOSTICS AND PRE-TREND TESTS
**********************************************************************/

* Joint test of all pre-treatment coefficients
testparm i.event_week_s#1.treated_BC if event_week < 0

* Support check: treated vs control by event week
tab event_week treated_BC if $IFBASE ///
    & inrange(event_week, -$WPRE, $WPOST)


/**********************************************************************
  14) EVENT-STUDY PLOT (TREATED INTERACTIONS ONLY)
      NOTE: x-axis is shifted event week; subtract ew_shift
            to recover true event week
**********************************************************************/

coefplot, ///
    keep(*1.treated_BC#*.event_week_s) ///
    vertical ///
    yline(0) ///
    title("Event Study: BC Minimum Wage Increase") ///
    subtitle("Baseline = week -1") ///
    xtitle("Event week (shifted; subtract ew_shift)") ///
    ytitle("Effect on ln(piece × prod)") ///
    ciopts(recast(rcap)) ///
    msymbol(O)















/*
/**********************************************************************
  4) Estimation: DiD (contract FE + worker FE)
     Outcome: ln_incentive = ln(piece * prod)
     Estimand: change in BC after June 1 relative to control provinces
**********************************************************************/

* Install reghdfe if needed
capture which reghdfe
if _rc {
    ssc install reghdfe, replace
    ssc install ftools, replace
}

reghdfe ln_incentive i.treated_BC##i.post ///
    if $IFBASE, ///
    absorb(contract name) vce(cluster contract)

* Coefficient of interest: 1.treated_BC#1.post


/**********************************************************************
  5) Estimation: Event study (weekly effects, baseline week = -1)
     - Window: [-$WPRE, +$WPOST]
     - Baseline week omitted: -1
     - Interpretation: each coefficient is (BC - control) in week k relative
       to week -1
**********************************************************************/

	
	
	
	
	
/**********************************************************************
  6) Optional: Plot event-study coefficients (coefplot)
**********************************************************************/

capture which coefplot
if _rc ssc install coefplot, replace

coefplot, ///
    keep(*.event_week#1.treated_BC) ///
    drop(*-1.event_week#1.treated_BC) ///
    vertical yline(0) ///
    title("Event study: BC vs $CONTROL ($TY)") ///
    xtitle("Event week (baseline = -1)") ///
    ytitle("Effect on ln(piece*prod)") ///
    ciopts(recast(rcap))

















/*

/**********************************************************************
  DERIVED GLOBALS (FULL, CLEAN RESET + REBUILD)
  Requires globals already set:
    $TY        (e.g., 2018 or 2019)
    $CONTROL   ("AB", "ON", or "ABON")
    $NEVERONLY (1 or 0)
**********************************************************************/

* --- 0) Hard reset so stale/broken globals cannot contaminate anything ---
macro drop IFBASE
macro drop TDATE
macro drop PROVCOND

* --- 1) Treatment date: June 1 of treated year ---
global TDATE = td(01jun$TY)

* --- 2) Province condition depends on CONTROL (build as its own global) ---
global PROVCOND 1==1
if "$CONTROL"=="AB"   global PROVCOND (province=="BC" | province=="AB")
if "$CONTROL"=="ON"   global PROVCOND (province=="BC" | province=="ON")
if "$CONTROL"=="ABON" global PROVCOND (province=="BC" | province=="AB" | province=="ON")

* Optional: fail loudly if CONTROL is mistyped
if !inlist("$CONTROL","AB","ON","ABON") {
    di as error "ERROR: CONTROL must be AB, ON, or ABON. You set CONTROL=$CONTROL"
    exit 198
}

* --- 3) Base sample restriction (year + province condition) ---
global IFBASE year==$TY & $PROVCOND

* --- 4) Optional never-top-up restriction ---
if $NEVERONLY==1 global IFBASE $IFBASE & never_topup==1

* --- 5) Sanity checks ---
di "TY      = $TY"
di "CONTROL = $CONTROL"
di "NEVER   = $NEVERONLY"
di "TDATE   = " %td $TDATE
di "IFBASE  = $IFBASE"

count if $IFBASE
tab province if $IFBASE


/*









/**********************************************************************
  DERIVED GLOBALS (robust version: no fragile quoting)
**********************************************************************/

* Treatment date (June 1 of treated year)
global TDATE = td(01jun$TY)

* Base if-condition: year and province set
global IFBASE "year==${TY}"

* Province restrictions depending on control group
if "$CONTROL"=="AB"   global IFBASE "$IFBASE & (province==""BC"" | province==""AB"")"
if "$CONTROL"=="ON"   global IFBASE "$IFBASE & (province==""BC"" | province==""ON"")"
if "$CONTROL"=="ABON" global IFBASE "$IFBASE & inlist(province,""BC"",""AB"",""ON"")"

* Optional never-top-up restriction
if $NEVERONLY==1 global IFBASE "$IFBASE & never_topup==1"

display "Running analysis with:"
display "  Treated year  = $TY"
display "  Control group = $CONTROL"
display "  Never-top-up  = $NEVERONLY"
display "  IF condition  = $IFBASE"


/**********************************************************************
  TREATMENT AND EVENT TIME FOR YEAR $TY
**********************************************************************/

capture drop treated_BC post event_day event_week
gen byte treated_BC = (province=="BC")
gen byte post = (date >= $TDATE) if year==$TY
gen int  event_day  = date - $TDATE if year==$TY
gen int  event_week = floor(event_day/7) if year==$TY


/**********************************************************************
  DiD: BC vs controls, year $TY
**********************************************************************/

reghdfe ln_incentive i.treated_BC##i.post ///
    if $IFBASE, ///
    absorb(contract name) vce(cluster contract
	
	



/*
import delimited using "/Users/danaandersen/Downloads/planter_mw.csv" , clear varnames(1)
 
g Date = date(date, "YMD")
 
egen id = group(name) 
 
duplicates drop id Date , force 
 
xtset id Date 

g bc = 1 if province == "BC"
	replace bc = 0 if bc == .
	
g min_wage_treat = mdy(6,1,2018) if bc == 1

g rel_date = Date - min_wage_treat 


eventdd lhrly_earn, timevar(rel_date) method(hdfe, absorb(id Date contract)) vce(cluster id )  leads(30) lags(15) accum



gen week = wofd(Date)
format week %tw

gen treat_w = wofd(mdy(6,1,2018)) if bc==1

* Relative WEEKS (treated only)
gen rel_w = week - treat_w if bc==1
label var rel_w "Weeks relative to minimum wage increase"

preserve
    collapse (mean) lhrly_earn ///
             (firstnm)  bc  rel_w contract, ///
             by(id week)

	
	* 3) Run event-study on weekly panel
    *-----------------------------------------
    * Choose window: about ±60 days -> use ±8 weeks
    local prew  = 6
    local postw = 6

	
    * Quick sanity checks
    bys id: assert bc == bc[1] if  rel_w>=`prew' & rel_w<=`postw'
    tab bc
	
    * eventdd requires leads/lags with inrange/accum/keepbal()
    eventdd lhrly_earn if rel_w>=`prew' & rel_w<=`postw' , ///
        timevar(rel_w) ///
        method(hdfe, absorb(id contract) vce(cluster id)) ///
        leads(`prew') lags(`postw') accum ///
        graph_op(xtitle("Weeks relative to treatment (t=0)") ///
                 ytitle("Effect on log hourly earnings") ///
                 yline(0, lpattern(dash)))
