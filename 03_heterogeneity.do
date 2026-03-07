/**********************************************************************
  03_HETEROGENEITY.DO
  Planter Minimum Wage Project

  PURPOSE:
    Section 3 heterogeneity analysis.

    Tests whether the productivity response to the BC minimum-wage
    reform scales with workers' pre-reform exposure to the wage floor.

  FRAMEWORK:
    See mw_exposure_framework_final.md for the full theoretical
    derivation. The core idea is:

      The MW floor binds whenever piece * prod < mw, i.e., hourly
      piece-rate earnings fall below the minimum wage. A worker's
      exposure is characterized by the distance between their expected
      earnings and the floor:

          d_ict = l r̂_ict − log(M_pre)

      where l r̂_ict is predicted log hourly piece-rate earnings under
      the pre-reform earnings process, and M_pre is the pre-reform
      provincial minimum wage.

      Workers with small d (near zero or negative) are at high risk of
      the floor binding. Workers with large d are rarely constrained.

    The main heterogeneity sorting variable is the pre-period mean of
    d_ict across a worker's pre-reform observations on each contract:

          d̄_ic = mean_pre(l r̂_ict) − log(M_pre)

    This incorporates persistent worker and contract earnings capacity
    plus average pre-period weather conditions (which affect actual
    earnings). It is entirely pre-determined and does not change
    mechanically after the reform.

    d̄_ic is standardized to d̃_ic (mean 0, SD 1) within the DiD
    estimation sample for interpretable coefficients.

  SECTIONS:
    3.0)  Clean rerun state
    3.1)  Settings / control panel (EDIT THESE)
    3.2)  Build controls + floor strings
    3.3)  Sample + event time
    3.4)  Earnings model (province-specific estimation)
    3.5)  Pre-period mean distance to floor (d̄_ic)
    3.6)  Residual variance (σ_resid)
    3.7)  Weather-driven variance (σ_weather)
    3.8)  Window flag (inwin_3), standardization, diagnostics
    3.9)  Parametric distance DiD (d̃_ic × treated × post)
    3.10) Nonparametric distance bins DiD
    3.11) Volatility heterogeneity DiD (σ_resid, σ_weather)
    3.12) Descriptive table: distance bins × pre/post (BC)
    3.13) Key estimates print (lincom)
    3.14) Binding probability (optional)

  OUTPUT:
    logs/03_heterogeneity_YYYYMMDD.log
    output/dbar_hist.png   (d̄_ic distribution diagnostic)

  REQUIRES: 00_data_prep.do must have been run first.
  REFERENCE: mw_exposure_framework_final.md
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
local logfile "logs/03_heterogeneity_`logdate'.log"
cap log close _all
log using "`logfile'", replace text
di as text "[03] Log opened: `logfile'"
di as text "[03] Started: $S_DATE $S_TIME"


/**********************************************************************
  LOAD
**********************************************************************/

use "planter_mw_prepped.dta", clear
di as text "[03] Loaded planter_mw_prepped.dta — N=" _N


/**********************************************************************
  3.0) CLEAN RERUN STATE
  Drop all variables and macros created by this do-file so it is
  safe to re-run without restarting Stata.
**********************************************************************/

cap estimates clear
cap program drop _trylincom

* Variables created by this do-file
foreach v in ///
    sample_3 treated_BC_3 post_3 event_week_3 event_week_3s ///
    sample_bc_earn sample_ab_earn ///
    resid_earn lr_hat lr_weather ///
    ln_Mpre d_ict use_dbar has_pre_ab_wc ///
    dbar_ic lrw_mean_c delta_weather ///
    sigma_resid_wc sigma_resid_w ///
    sigma_weather_c sigma_weather_i ///
    inwin_3 ///
    dtilde_ic sigma_resid_z sigma_weather_z ///
    dbar_bin ///
    z_ict pi_hat {
    cap drop `v'
}

* Macros set by this do-file
foreach m in ///
    TY_3 CONTROL_3 Y_HET ///
    CLUSTER_3 VCE_3 ///
    USE_WEATHER_3 USE_HOURS_3 USE_DAYFLAGS_3 USE_LIFECYCLE_3 ///
    USE_FLOORS_3 MIN_HRS_3 MIN_TREES_3 ///
    EARN_WEATHER EARN_LIFECYCLE EARN_FLOORS EARN_MIN_HRS EARN_MIN_TREES ///
    WPRE_3 WPOST_3 ///
    DBAR_BINS_N DO_BINDING_PROB ///
    DID_CONTRACTFE_3 DID_WEEKFE_3 DID_POST_MAIN_3 DID_POSTX_3 ///
    X_3 LIFECYCLE_3 SAMPLE_OK_3 EARN_W EARN_Z {
    cap macro drop `m'
}


/**********************************************************************
  3.1) CONTROL PANEL — EDIT THESE
**********************************************************************/

di as text "============================================================"
di as text "[3.1] Control panel"
di as text "============================================================"

* ----- Treatment year and control province -----
* Must match 02_main_analysis.do for consistency.
global TY_3      2018
global CONTROL_3 "AB"

* ----- Outcome for heterogeneity DiD regressions -----
* Options: ln_incentive  lprod  ln_piece  tot_hrs_day
* ln_incentive is the primary outcome and matches the earnings model.
global Y_HET "ln_incentive"

* ----- Clustering -----
global CLUSTER_3 "worker_id"
global VCE_3     "vce(cluster $CLUSTER_3)"

* ----- Controls for DiD regressions (Sections 3.9–3.11) -----
* These apply to the het DiD regressions, not the earnings model.
* Experience is excluded: cross-season experience is constant within
* a single year, so it is perfectly collinear with worker FE and would
* be silently dropped by reghdfe. See PROJECT_NOTES.md.
global USE_WEATHER_3   1    // temp, wind, precip (quadratic)
global USE_HOURS_3     0    // exclude hours — endogenous to MW response
global USE_DAYFLAGS_3  0    // exclude multi_contract_day flag
global USE_LIFECYCLE_3 1    // contract ramp + season ramp + ramp×cum_seasons

* ----- Floor filters for DiD regressions -----
* Set to match Section 2 (02_main_analysis.do) for sample comparability.
* Comparing het DiD results to the baseline 2A DiD requires the same
* sample; mismatched floors make the comparison ambiguous.
global USE_FLOORS_3  1
global MIN_HRS_3     2     // minimum total hours/day (tot_hrs_day)
global MIN_TREES_3   100   // minimum trees this observation (prod * tot_hrs)

* ----- Earnings model specification (Section 3.4) -----
* Controls for the earnings model may differ from the DiD controls.
* The earnings model is estimated on the pre-period (BC) or pre+post (AB)
* and uses reghdfe with absorbed worker and contract FEs.
global EARN_WEATHER   1    // include weather controls in earnings model
global EARN_LIFECYCLE 1    // include lifecycle controls in earnings model

* Floors for the earnings model (separate from DiD floors).
* Turning floors off (EARN_FLOORS=0) maximizes the estimation sample
* and improves FE precision, since outlier influence is absorbed by FEs.
* Floors can be turned on for robustness.
global EARN_FLOORS    0
global EARN_MIN_HRS   2
global EARN_MIN_TREES 100

* ----- Event window for Section 3 DiD regressions -----
* Fixed window (not auto-selected). Symmetry with Section 2A window
* is desirable for comparability, but the window here may be narrower
* to maximize sample balance around June 1.
* WPRE_3: weeks of pre-period to include (e.g. 3 = weeks -3 to -1)
* WPOST_3: weeks of post-period to include (e.g. 8 = weeks 0 to 8)
global WPRE_3  3
global WPOST_3 8

* ----- Distance bins for nonparametric het DiD (Section 3.10) -----
global DBAR_BINS_N 4    // 4 = quartiles; increase for finer splits

* ----- Binding probability (Section 3.14) -----
* Setting to 0 skips the optional binding probability computation.
* Computationally inexpensive but adds output that may be noisy with
* small samples; run after main results are stable.
global DO_BINDING_PROB 0

* ----- DiD regression spec toggles (Sections 3.9–3.11) -----
* DID_CONTRACTFE_3: absorb contract_id in all het DiD regressions.
*   Recommended = 1: removes time-invariant contract-level confounders.
global DID_CONTRACTFE_3 1

* DID_WEEKFE_3: absorb event_week_3s in all het DiD regressions.
*   Recommended = 1: removes common within-season time trends.
global DID_WEEKFE_3 1

* DID_POST_MAIN_3: include 1.post_3 as a main effect in het DiD.
*   With worker FE absorbed, the level effect of post is absorbed;
*   but when DID_WEEKFE_3=0, a post main effect may be needed.
*   With DID_WEEKFE_3=1 this is redundant; set 0 to avoid collinearity.
global DID_POST_MAIN_3 0

* DID_POSTX_3: include common post × exposure slopes.
*   1 = include c.dtilde_ic#1.post_3 (and equivalent for σ terms).
*   This identifies whether the exposure gradient exists in the control
*   group (AB) post-period as well, beyond the BC treatment effect.
*   Recommended = 1 for transparency; set 0 for simpler spec.
global DID_POSTX_3 1


/**********************************************************************
  3.2) BUILD CONTROLS + FLOOR STRINGS
**********************************************************************/

di as text "------------------------------------------------------------"
di as text "[3.2] Building globals from control panel..."
di as text "------------------------------------------------------------"

* ---- Controls for DiD regressions (X_3, LIFECYCLE_3) ----
global X_3 ""
if $USE_HOURS_3==1    global X_3 "$X_3 tot_hrs_day"
if $USE_DAYFLAGS_3==1 global X_3 "$X_3 multi_contract_day"
if $USE_WEATHER_3==1  global X_3 "$X_3 c.temp##c.temp c.wind_spd##c.wind_spd precip0 c.precip_pos##c.precip_pos"

global LIFECYCLE_3 ""
if $USE_LIFECYCLE_3==1 {
    global LIFECYCLE_3 "$LIFECYCLE_3 c.contract_day##c.contract_day"
    global LIFECYCLE_3 "$LIFECYCLE_3 c.season_days_worked##c.season_days_worked"
    global LIFECYCLE_3 "$LIFECYCLE_3 c.season_days_worked##c.cum_seasons"
}

* ---- Controls for the earnings model (EARN_W, EARN_Z) ----
* The earnings model absorbs worker and contract FEs (reghdfe absorb).
* Weather and lifecycle controls are the non-absorbed regressors.
* Note: experience is also excluded here — constant within worker in
* a single year, collinear with worker FE.
global EARN_W ""
if $EARN_WEATHER==1 ///
    global EARN_W "c.temp##c.temp c.wind_spd##c.wind_spd precip0 c.precip_pos##c.precip_pos"

global EARN_Z ""
if $EARN_LIFECYCLE==1 {
    global EARN_Z "$EARN_Z c.contract_day##c.contract_day"
    global EARN_Z "$EARN_Z c.season_days_worked##c.season_days_worked"
    global EARN_Z "$EARN_Z c.season_days_worked##c.cum_seasons"
}

* ---- Floor filter for DiD regressions ----
global SAMPLE_OK_3 "1==1"
if $USE_FLOORS_3==1 {
    global SAMPLE_OK_3 "$SAMPLE_OK_3 & tot_hrs_day >= $MIN_HRS_3"
    global SAMPLE_OK_3 "$SAMPLE_OK_3 & piece > 0 & prod > 0"
    global SAMPLE_OK_3 "$SAMPLE_OK_3 & (prod * tot_hrs) >= $MIN_TREES_3"
}

* ---- Safe numeric locals from globals (for use in absorb lists) ----
* reghdfe absorb() lists cannot reference globals with ==1 syntax.
local did_contractfe = real(trim("$DID_CONTRACTFE_3"))
local did_weekfe     = real(trim("$DID_WEEKFE_3"))
local did_post_main  = real(trim("$DID_POST_MAIN_3"))
local did_postx      = real(trim("$DID_POSTX_3"))

* ---- Debug output ----
di as text "[3.2] Y_HET:         $Y_HET"
di as text "[3.2] TY_3:          $TY_3   CONTROL_3: $CONTROL_3"
di as text "[3.2] X_3:           $X_3"
di as text "[3.2] LIFECYCLE_3:   $LIFECYCLE_3"
di as text "[3.2] EARN_W:        $EARN_W"
di as text "[3.2] EARN_Z:        $EARN_Z"
di as text "[3.2] SAMPLE_OK_3:   $SAMPLE_OK_3"
di as text "[3.2] Window:        [-$WPRE_3, $WPOST_3]"
di as text "[3.2] DBAR_BINS_N:   $DBAR_BINS_N"
di as text "[3.2] DID_CONTRACTFE: $DID_CONTRACTFE_3   DID_WEEKFE: $DID_WEEKFE_3"
di as text "[3.2] DID_POST_MAIN:  $DID_POST_MAIN_3   DID_POSTX:  $DID_POSTX_3"


/**********************************************************************
  3.3) DEFINE SAMPLE + EVENT TIME

  sample_3:     BC and AB observations in the treated year (TY_3).
  treated_BC_3: =1 for BC, =0 for AB within sample_3.
  post_3:       =1 from June 1 onward (within sample_3).
  event_week_3: weeks relative to June 1 (0 = week of June 1).
  event_week_3s: event_week_3 + 100 to avoid negative factor values.

  ABS_DID_3: the absorbed FE list for all Section 3 DiD regressions,
  built from the DID_* toggles. Always includes worker_id (baseline).
  Optionally adds contract_id and event_week_3s.

  Note: event_week_3s is constructed for all sample_3 obs (even
  outside the window) so that it is non-missing for any obs that might
  enter the DiD. The DiD restricts to inwin_3, so only window-obs
  are used in estimation.
**********************************************************************/

di as text "------------------------------------------------------------"
di as text "[3.3] Sample + event time"
di as text "------------------------------------------------------------"

gen byte sample_3 = (year == $TY_3) & inlist(province_s, "BC", "$CONTROL_3")

gen byte treated_BC_3 = .
replace treated_BC_3 = (province_s == "BC") if sample_3==1

local tdate_3 = mdy(6, 1, $TY_3)
gen byte post_3       = (date >= `tdate_3')            if sample_3==1
gen int  event_week_3 = floor((date - `tdate_3') / 7) if sample_3==1
gen int  event_week_3s = event_week_3 + 100            if sample_3==1

di as text "[3.3] Sample counts:"
tab province_s if sample_3==1
tab post_3     if sample_3==1

* ---- Unified absorb list for ALL Section 3 DiD regressions ----
* Built from local toggles (numeric) not globals, to avoid macro
* evaluation issues inside absorb() specifications.
local ABS_DID_3 "worker_id"
if `did_contractfe'==1 local ABS_DID_3 "`ABS_DID_3' contract_id"
if `did_weekfe'==1     local ABS_DID_3 "`ABS_DID_3' event_week_3s"
di as text "[3.3] ABS_DID_3: `ABS_DID_3'"


/**********************************************************************
  3.4) EARNINGS MODEL

  Estimate log hourly piece-rate earnings as a function of worker FE,
  contract FE, weather, and lifecycle controls. Province-specific
  estimation samples (see mw_exposure_framework_final.md, Section 4):

    BC: pre-period only (post_3==0 & province=="BC")
        Rationale: avoids contamination from the MW reform response.
        After June 1 in BC, both piece rates and productivity change
        as firms and workers respond to the new wage floor. Including
        post-period BC data would bias the pre-reform earnings process.

    AB: pre + post (province=="AB")
        Rationale: AB has only 1 of 6 contracts active before June 1,
        2018 (contract B-18-0408, N=331 obs). Using pre-period BC data
        alone would leave 5 AB contracts (N≈3,320 obs) with unidentified
        contract FEs. Because AB is the control province (no MW reform),
        post-period AB observations do not introduce policy contamination.

  Outcome: ln_incentive = ln(piece * prod)
    Using incentive rather than productivity alone means the worker FE
    α̂_i absorbs baseline earnings capacity (ability × typical piece
    rate) and the contract FE γ̂_c absorbs each contract's piece rate
    level and inherent difficulty. Predicted earnings l r̂_ict therefore
    represent expected hourly earnings for worker i on contract c on
    day t under the pre-reform earnings process — the exact quantity
    needed for distance-to-floor calculations (see framework Section 3).

  Outputs (combined BC + AB):
    resid_earn   û_ict: idiosyncratic earnings shocks
    lr_hat       l r̂_ict = ln_incentive − û_ict (in-sample only)
    lr_weather   W_ict β̂: weather component only, computed from stored
                 regression coefficients immediately after each
                 province-specific regression (before coefficients are
                 overwritten by the next regression). Used for σ_weather.

  Note on experience: cross-season experience is constant within a
  single year (2018). It would be collinear with worker_id FE and
  silently dropped. It is excluded from EARN_W and EARN_Z intentionally.

  Note on separate β̂ by province: estimating separately for BC and
  AB yields province-specific weather and lifecycle coefficients.
  This is intentional — we do not impose that weather affects earnings
  identically across provinces. The cost is that predicted earnings
  are not on a perfectly comparable scale across provinces, but this
  does not affect het-sort rankings within each province.
**********************************************************************/

di as text "------------------------------------------------------------"
di as text "[3.4] Earnings model"
di as text "------------------------------------------------------------"

* ---- Define estimation samples ----
cap drop sample_bc_earn sample_ab_earn

* BC: pre-period observations only
gen byte sample_bc_earn = (sample_3==1) & (province_s=="BC") & (post_3==0) ///
    & !missing(ln_incentive)
if $EARN_FLOORS==1 {
    replace sample_bc_earn = 0 if tot_hrs_day < $EARN_MIN_HRS
    replace sample_bc_earn = 0 if piece <= 0 | prod <= 0
    replace sample_bc_earn = 0 if (prod * tot_hrs) < $EARN_MIN_TREES
}

* AB: full 2018 season (pre + post, untreated)
gen byte sample_ab_earn = (sample_3==1) & (province_s=="AB") ///
    & !missing(ln_incentive)
if $EARN_FLOORS==1 {
    replace sample_ab_earn = 0 if tot_hrs_day < $EARN_MIN_HRS
    replace sample_ab_earn = 0 if piece <= 0 | prod <= 0
    replace sample_ab_earn = 0 if (prod * tot_hrs) < $EARN_MIN_TREES
}

di as text "[3.4] Estimation sample sizes:"
di as text "[3.4]   BC pre-period: " %7.0f `=r(N)' _skip(2) "(will update below)"
count if sample_bc_earn==1
di as text "[3.4]   BC pre-period: " r(N)
count if sample_ab_earn==1
di as text "[3.4]   AB pre+post:   " r(N)

* ---- Initialize output variables ----
cap drop resid_earn lr_hat lr_weather
gen double resid_earn = .
gen double lr_hat     = .
gen double lr_weather = .

* ============================================================
* BC EARNINGS MODEL
* Estimated on BC pre-period only.
* ============================================================

di as text "[3.4] --- BC earnings model ---"
di as text "[3.4]   Outcome:  ln_incentive"
di as text "[3.4]   Absorb:   worker_id contract_id"
di as text "[3.4]   Controls: $EARN_W $EARN_Z"

reghdfe ln_incentive $EARN_W $EARN_Z if sample_bc_earn==1, ///
    absorb(worker_id contract_id) ///
    resid(resid_bc_tmp) ///
    vce(cluster $CLUSTER_3)

di as text "[3.4] BC model fit: N=" e(N) "  R2_within=" %6.4f e(r2_within)

* Store residuals and fitted values for BC observations.
* l r̂_ict = ln_incentive − û_ict (equivalently: actual minus residual).
replace resid_earn = resid_bc_tmp      if sample_bc_earn==1
replace lr_hat     = ln_incentive - resid_bc_tmp if sample_bc_earn==1

* Compute the weather component W_ict β̂_BC using the stored BC coefficients.
* This must be done NOW, before running the AB regression which overwrites _b[].
*
* The weather component is the linear combination of weather predictors
* at the BC-estimated coefficients. Coefficient names follow Stata's
* factor-variable notation: _b[c.temp] = linear temp coefficient,
* _b[c.temp#c.temp] = quadratic temp coefficient, etc.
*
* If EARN_WEATHER=0, no weather predictors were in the model; set to 0.
if $EARN_WEATHER==1 {
    replace lr_weather = ///
        _b[c.temp]                      * (temp)           + ///
        _b[c.temp#c.temp]               * (temp * temp)    + ///
        _b[c.wind_spd]                  * (wind_spd)       + ///
        _b[c.wind_spd#c.wind_spd]       * (wind_spd * wind_spd) + ///
        _b[precip0]                     * (precip0)        + ///
        _b[c.precip_pos]                * (precip_pos)     + ///
        _b[c.precip_pos#c.precip_pos]   * (precip_pos * precip_pos) ///
        if sample_bc_earn==1
    di as text "[3.4] BC weather component stored."
}
else {
    replace lr_weather = 0 if sample_bc_earn==1
    di as text "[3.4] EARN_WEATHER=0 — weather component set to 0 for BC."
}

cap drop resid_bc_tmp

* ============================================================
* AB EARNINGS MODEL
* Estimated on AB pre + post (untreated control province).
* ============================================================

di as text "[3.4] --- AB earnings model ---"
di as text "[3.4]   Outcome:  ln_incentive"
di as text "[3.4]   Absorb:   worker_id contract_id"
di as text "[3.4]   Controls: $EARN_W $EARN_Z"
di as text "[3.4]   Note: AB has 1/6 contracts active pre-June 1."
di as text "[3.4]         Post-period AB data used to identify all 6 contract FEs."

reghdfe ln_incentive $EARN_W $EARN_Z if sample_ab_earn==1, ///
    absorb(worker_id contract_id) ///
    resid(resid_ab_tmp) ///
    vce(cluster $CLUSTER_3)

di as text "[3.4] AB model fit: N=" e(N) "  R2_within=" %6.4f e(r2_within)

* Store residuals and fitted values for AB observations.
replace resid_earn = resid_ab_tmp      if sample_ab_earn==1
replace lr_hat     = ln_incentive - resid_ab_tmp if sample_ab_earn==1

* Compute weather component W_ict β̂_AB using AB-specific coefficients.
if $EARN_WEATHER==1 {
    replace lr_weather = ///
        _b[c.temp]                      * (temp)           + ///
        _b[c.temp#c.temp]               * (temp * temp)    + ///
        _b[c.wind_spd]                  * (wind_spd)       + ///
        _b[c.wind_spd#c.wind_spd]       * (wind_spd * wind_spd) + ///
        _b[precip0]                     * (precip0)        + ///
        _b[c.precip_pos]                * (precip_pos)     + ///
        _b[c.precip_pos#c.precip_pos]   * (precip_pos * precip_pos) ///
        if sample_ab_earn==1
    di as text "[3.4] AB weather component stored."
}
else {
    replace lr_weather = 0 if sample_ab_earn==1
    di as text "[3.4] EARN_WEATHER=0 — weather component set to 0 for AB."
}

cap drop resid_ab_tmp

* ---- Summary diagnostics ----
di as text "[3.4] lr_hat summary (BC pre + AB pre+post):"
summ lr_hat if (sample_bc_earn==1 | sample_ab_earn==1)

di as text "[3.4] resid_earn summary:"
summ resid_earn if (sample_bc_earn==1 | sample_ab_earn==1)


/**********************************************************************
  3.5) PRE-PERIOD MEAN DISTANCE TO FLOOR (d̄_ic)

  Construct the main heterogeneity sorting variable:

      d̄_ic = mean_pre(l r̂_ict) − log(M_pre)
             = mean_pre(lr_hat) − log(M_pre)

  This is the pre-period mean of the realized daily distance d_ict,
  where d_ict = lr_hat − log(M_pre).

  KEY DESIGN DECISIONS (see framework Sections 2 and 7.1):

  1. M_pre, not M_t:
     We use the PRE-REFORM minimum wage (M_pre), not the contemporaneous
     wage. This ensures d̄_ic does not mechanically shift downward for
     BC workers after June 1 (when M_t increases). d̄_ic is a pure
     pre-reform characteristic — it measures how close a worker was to
     the floor BEFORE the policy changed.

  2. Include weather in d̄_ic:
     d̄_ic includes the average weather component W̄_c,pre β̂ from the
     pre-period. Because weather is measured at the contract level (one
     Env Canada station per camp), the weather average is the same for
     all workers on the same contract. Including weather means d̄_ic
     reflects the worker's full earnings environment pre-reform —
     not just their ability and the contract's baseline piece rate,
     but also whether they were on a contract in systematically good or
     bad weather conditions.

  3. Aggregation level: worker × contract (wc_id):
     d̄_ic varies across wc pairs. Workers who work multiple contracts
     within the season have multiple d̄_ic values. The DiD regressions
     use the wc-level variable and absorb worker + contract FEs, so
     het effects are identified from within-province cross-wc variation.

  AGGREGATION FOR AB POST-ONLY CONTRACTS:
     AB has 5 contracts with no pre-June-1 observations (see 3.4 note).
     For workers on those contracts, there are no pre-period lr_hat
     values to average. Because AB is untreated, we fall back to using
     post-period lr_hat values for those workers only. This avoids
     losing the majority of the AB control group from the het analysis
     while maintaining the pre-determined nature of d̄_ic for all BC
     workers and for the AB workers with pre-period coverage.
**********************************************************************/

di as text "------------------------------------------------------------"
di as text "[3.5] Pre-period mean distance d̄_ic"
di as text "------------------------------------------------------------"

* ---- Step 1: Compute log(M_pre) by province ----
* M_pre = the minimum wage before June 1 in each province.
* In the pre-period (post_3==0), mw is constant within province
* (no within-season changes before the June 1 reform date). We take
* the mean of mw over pre-period observations as the M_pre scalar.
* For BC: M_pre_bc = BC minimum wage before the June 1, 2018 increase.
* For AB: M_pre_ab = AB minimum wage (stable throughout 2018).
cap drop ln_Mpre

quietly summ mw if sample_3==1 & province_s=="BC" & post_3==0, meanonly
scalar M_pre_bc = r(mean)
quietly summ mw if sample_3==1 & province_s=="AB" & post_3==0, meanonly
scalar M_pre_ab = r(mean)

di as text "[3.5] M_pre_BC = " M_pre_bc "  log = " %6.4f log(M_pre_bc)
di as text "[3.5] M_pre_AB = " M_pre_ab "  log = " %6.4f log(M_pre_ab)

* ln_Mpre: log of the pre-reform provincial MW for each observation.
* Applied UNIFORMLY to pre AND post period observations —
* the reform does not change this variable for BC post-period obs.
* This is the core of the "pre-determined" exposure design.
gen double ln_Mpre = .
replace ln_Mpre = log(M_pre_bc) if sample_3==1 & province_s=="BC"
replace ln_Mpre = log(M_pre_ab) if sample_3==1 & province_s=="AB"

* ---- Step 2: Realized daily distance d_ict ----
* d_ict = l r̂_ict − log(M_pre): the predicted log gap between a
* worker's expected earnings and the pre-reform floor on each day.
* Positive d_ict → expected earnings exceed floor (low binding risk).
* Negative d_ict → expected earnings fall below floor (high risk).
* Near-zero d_ict → on the margin (highest uncertainty about binding).
cap drop d_ict
gen double d_ict = lr_hat - ln_Mpre if !missing(lr_hat) & !missing(ln_Mpre)

di as text "[3.5] d_ict distribution (all obs with lr_hat):"
summ d_ict if sample_3==1, detail

* ---- Step 3: Flag which observations to use for d̄_ic ----
* use_dbar = 1 for observations that contribute to the pre-period mean.
cap drop use_dbar has_pre_ab_wc
gen byte use_dbar = 0

* BC: use pre-period obs within the BC earnings model sample.
* (post_3==0 & sample_bc_earn ensures: BC, pre-period, valid d_ict)
replace use_dbar = 1 if sample_bc_earn==1 & post_3==0 & !missing(d_ict)

* AB: prefer pre-period obs (post_3==0) where available.
replace use_dbar = 1 if sample_ab_earn==1 & post_3==0 & !missing(d_ict)

* AB fallback: for wc pairs with NO pre-period observations (post-only
* AB contracts), use post-period lr_hat as a substitute.
* Step 1: identify which wc pairs have at least one pre-period obs.
bysort wc_id: egen byte has_pre_ab_wc = ///
    max(cond(sample_ab_earn==1 & post_3==0 & !missing(d_ict), 1, 0))
* Step 2: for post-only wc pairs, allow post-period obs to contribute.
replace use_dbar = 1 if sample_ab_earn==1 & post_3==1 ///
    & !missing(d_ict) & has_pre_ab_wc==0

di as text "[3.5] use_dbar counts by province × period:"
tab province_s post_3 if use_dbar==1

* ---- Step 4: Compute d̄_ic = pre-period mean of d_ict by wc_id ----
cap drop dbar_ic
bysort wc_id: egen double dbar_ic = mean(cond(use_dbar==1, d_ict, .))

di as text "[3.5] d̄_ic coverage:"
count if sample_3==1 & province_s=="BC" & !missing(dbar_ic)
di as text "[3.5]   BC observations with d̄_ic: " r(N)
count if sample_3==1 & province_s=="AB" & !missing(dbar_ic)
di as text "[3.5]   AB observations with d̄_ic: " r(N)

di as text "[3.5] d̄_ic summary by province:"
summ dbar_ic if sample_3==1 & province_s=="BC"
summ dbar_ic if sample_3==1 & province_s=="AB"


/**********************************************************************
  3.6) RESIDUAL VARIANCE (σ_resid)

  Construct worker-contract-level earnings volatility from the
  idiosyncratic residuals û_ict from the earnings model (Section 3.4).

  σ_resid_wc = SD(û_ict) for worker i on contract c, over the pre-
               period observations used to compute d̄_ic (use_dbar==1).

  Interpretation: σ_resid_wc captures the day-to-day variability in
  a worker's earnings on a given contract that is unexplained by
  persistent worker ability, contract characteristics, weather, or
  lifecycle controls. This is the idiosyncratic shock component.

  Workers with high σ_resid face more random variation in earnings
  and are therefore at greater risk of the floor binding on any given
  day, even holding expected earnings constant.

  AGGREGATION NOTE:
  - Worker-contract level is preferred (see framework Section 8) if
    the number of pre-period observations per wc pair is sufficient to
    estimate a meaningful SD. Pairs with fewer than ~5 obs will have
    very noisy σ_resid_wc.
  - Worker-level (σ_resid_w) is a fallback that pools across contracts
    for more stable estimation but loses within-worker variation across
    contracts.
  Both are constructed; the DiD regression uses σ_resid_wc as baseline
  and σ_resid_w as robustness.
**********************************************************************/

di as text "------------------------------------------------------------"
di as text "[3.6] Residual variance σ_resid"
di as text "------------------------------------------------------------"

cap drop sigma_resid_wc sigma_resid_w

* Worker-contract level: SD of residuals over pre-period (use_dbar==1) obs.
bysort wc_id: egen double sigma_resid_wc = ///
    sd(cond(use_dbar==1 & !missing(resid_earn), resid_earn, .))

* Worker level: SD of residuals over pre-period obs, pooled across contracts.
bysort worker_id: egen double sigma_resid_w = ///
    sd(cond(use_dbar==1 & !missing(resid_earn), resid_earn, .))

* Diagnostic: how many wc pairs have a valid σ_resid_wc?
* SD requires at least 2 obs; wc pairs with only 1 pre-period obs
* will have σ_resid_wc = missing.
di as text "[3.6] σ_resid_wc coverage:"
count if sample_3==1 & province_s=="BC" & !missing(sigma_resid_wc)
di as text "[3.6]   BC obs with σ_resid_wc: " r(N)
count if sample_3==1 & province_s=="AB" & !missing(sigma_resid_wc)
di as text "[3.6]   AB obs with σ_resid_wc: " r(N)

di as text "[3.6] σ_resid_wc summary (all obs in sample_3):"
summ sigma_resid_wc if sample_3==1


/**********************************************************************
  3.7) WEATHER-DRIVEN VARIANCE (σ_weather)

  Construct a measure of how much a worker's earnings vary due to
  weather fluctuations on their contract.

  CONSTRUCTION:
  1. Weather deviation from contract-level pre-period mean:
       Δ_weather_ict = W_ict β̂ − W̄_c,pre β̂
     where W_ict β̂ = lr_weather (stored in Section 3.4) and
     W̄_c,pre β̂ is the pre-period average weather component for
     contract c.

  2. Contract-level weather volatility:
       σ_weather_c = SD(Δ_weather_ict) for contract c, pre-period.
     Since weather is measured at the contract level (one station per
     camp), σ_weather_c is a contract-level scalar — constant across
     all workers on contract c.

  3. Worker-level weather volatility:
       σ_weather_i = mean(σ_weather_c) across contracts worked by
                     worker i in the pre-period.
     This averages the contract-level volatility across the contracts
     a worker appeared on, weighted by days worked (since the mean is
     taken over all observations for that worker).

  IDENTIFICATION NOTE (see framework Section 9):
  In DiD regressions with contract FE absorbed, σ_weather_c is
  constant within contract and its level is absorbed by contract FE.
  The interaction σ_weather_c × treated × post is formally identified
  from cross-contract differences in weather volatility interacted with
  the BC reform. However, identification relies on cross-contract
  variation in σ_weather_c, which may be limited with ~25 BC contracts.

  σ_weather_i (worker-level) varies across workers even within the same
  contract if they worked different sets of contracts across the season.
  It is the preferred variable for DiD regressions (more stable).
**********************************************************************/

di as text "------------------------------------------------------------"
di as text "[3.7] Weather-driven variance σ_weather"
di as text "------------------------------------------------------------"

cap drop lrw_mean_c delta_weather sigma_weather_c sigma_weather_i

* ---- Step 1: Contract-level pre-period mean weather component ----
* Average weather effect (W̄_c,pre β̂) for each contract, pre-period only.
bysort contract_id: egen double lrw_mean_c = ///
    mean(cond(sample_3==1 & post_3==0 & !missing(lr_weather), lr_weather, .))

* ---- Step 2: Daily weather deviation from contract mean ----
* Δ_weather_ict = W_ict β̂ − W̄_c,pre β̂
gen double delta_weather = lr_weather - lrw_mean_c ///
    if sample_3==1 & !missing(lr_weather) & !missing(lrw_mean_c)

* ---- Step 3: Contract-level SD of weather deviations (pre-period) ----
bysort contract_id: egen double sigma_weather_c = ///
    sd(cond(sample_3==1 & post_3==0 & !missing(delta_weather), delta_weather, .))

* ---- Step 4: Worker-level mean of contract weather volatility ----
* For each worker, average the contract weather volatility across their
* pre-period contracts. Since sigma_weather_c is constant within a
* contract, this averages across the distinct contracts the worker
* appeared on, weighted by days worked on each contract.
bysort worker_id: egen double sigma_weather_i = ///
    mean(cond(sample_3==1 & post_3==0 & !missing(sigma_weather_c), sigma_weather_c, .))

di as text "[3.7] σ_weather_c summary (contract level, pre-period):"
summ sigma_weather_c if sample_3==1 & post_3==0

di as text "[3.7] σ_weather_i summary (worker level):"
summ sigma_weather_i if sample_3==1 & province_s=="BC"
summ sigma_weather_i if sample_3==1 & province_s=="AB"


/**********************************************************************
  3.8) WINDOW FLAG (inwin_3), STANDARDIZATION (d̃_ic), DIAGNOSTICS

  inwin_3: the analysis sample for all Section 3 DiD regressions.
    Applies: event window, floor filters, and non-missing d̄_ic.

  d̃_ic: standardized d̄_ic (mean 0, SD 1) within the inwin_3 sample.
    Standardizing within the DiD sample rather than the full dataset
    ensures the coefficient on d̃_ic × treated × post has a consistent
    interpretation: the difference in treatment effects between workers
    one SD apart in pre-reform floor exposure, among the workers
    actually used in the DiD regression.

  dbar_bin: quartile of d̄_ic within the inwin_3 sample.
    Workers are sorted into DBAR_BINS_N equal groups by pre-reform
    exposure. Bin 1 = closest to the floor (highest exposure risk).
    Bin DBAR_BINS_N = furthest from the floor (lowest risk).
    The bin ordering follows the economic prediction: effects should
    be largest in bin 1 if the MW floor drives productivity changes.

  Distribution diagnostic: histogram of d̄_ic in the inwin_3 sample.
**********************************************************************/

di as text "------------------------------------------------------------"
di as text "[3.8] Window flag, standardization, diagnostics"
di as text "------------------------------------------------------------"

* ---- Define inwin_3: analysis sample for DiD regressions ----
cap drop inwin_3
gen byte inwin_3 = (sample_3==1) ///
    & inrange(event_week_3, -$WPRE_3, $WPOST_3) ///   // within event window
    & ($SAMPLE_OK_3) ///                               // floor filters
    & !missing($Y_HET) ///                             // non-missing outcome
    & !missing(dbar_ic)                                // non-missing exposure

di as text "[3.8] inwin_3 sample:"
tab province_s post_3 if inwin_3==1

* ---- Standardize d̄_ic within inwin_3 ----
* Compute mean and SD of d̄_ic in the DiD sample.
* The standardized version d̃_ic has the same value for all observations
* belonging to the same wc pair (since d̄_ic is wc-level constant).
cap drop dtilde_ic

quietly summ dbar_ic if inwin_3==1
scalar dbar_mu = r(mean)
scalar dbar_sd = r(sd)

di as text "[3.8] d̄_ic in inwin_3: mean=" %6.3f dbar_mu "  SD=" %6.3f dbar_sd

if dbar_sd > 0 {
    gen double dtilde_ic = (dbar_ic - dbar_mu) / dbar_sd if !missing(dbar_ic)
    di as text "[3.8] d̃_ic constructed (mean 0, SD 1 within inwin_3)."
}
else {
    di as error "[3.8] WARN: dbar_sd = 0. d̃_ic not constructed."
    gen double dtilde_ic = 0 if !missing(dbar_ic)
}

* ---- Create distance bins (dbar_bin) ----
* Quartile assignment from the inwin_3 sample.
* Bin 1 = lowest d̄_ic = closest to floor = highest binding risk.
* Post-assignment: propagate bin to all observations for the same wc_id
* (since the bin is a wc-level characteristic, workers appearing outside
* inwin_3 in some obs still get the bin from their inwin_3 obs).
cap drop dbar_bin
xtile dbar_bin_tmp = dbar_ic if inwin_3==1, nq($DBAR_BINS_N)
bysort wc_id: egen byte dbar_bin = max(dbar_bin_tmp)
cap drop dbar_bin_tmp

di as text "[3.8] d̄_ic bin distribution in inwin_3:"
tab dbar_bin if inwin_3==1

* ---- Distribution diagnostic: histogram of d̄_ic ----
* Exported to output/dbar_hist.png for visual inspection.
* A roughly bell-shaped distribution is expected. A spike at 0 or
* a very left-skewed distribution would warrant investigation.
histogram dbar_ic if inwin_3==1, ///
    width(0.1) frequency ///
    title("Pre-period mean distance d̄_ic (inwin_3 sample)", size(medium)) ///
    subtitle("Positive = above floor; Negative = at/below floor", size(small)) ///
    xtitle("d̄_ic = mean_pre(l r̂_ict) − log(M_pre)") ///
    ytitle("Frequency") ///
    xline(0, lpattern(dash) lcolor(red)) ///
    graphregion(color(white))
graph export "output/dbar_hist.png", replace
di as text "[3.8] Histogram saved: output/dbar_hist.png"


/**********************************************************************
  3.9) PARAMETRIC DISTANCE DiD

  Tests whether the treatment effect on ln_incentive varies
  continuously with d̃_ic (standardized pre-reform floor exposure).

  Model:
    Y = β_0 (treated × post)
      + [β_1 (post)]                       if DID_POST_MAIN_3=1
      + β_2 (d̃_ic × treated × post)       ← key het slope
      + [β_3 (d̃_ic × post)]               if DID_POSTX_3=1
      + controls + absorbed FEs + ε

  Interpretation of β_2:
    The difference in treatment effect between workers one SD apart in
    pre-reform floor exposure. Under the theory that effects are larger
    for workers closer to the floor (lower d̄_ic), β_2 > 0 (more exposed
    workers = lower d̃_ic = larger negative effect, so β_2 coefficient
    on d̃_ic × treated × post should be positive if being FURTHER from
    the floor reduces the effect).

  β_3 (if included) tests whether d̃_ic predicts a different pre/post
  change in the CONTROL group (AB). A significant β_3 would indicate
  that d̃_ic captures something that trends differently across exposure
  levels even in the absence of a reform — a parallel-trends concern.

  Collinearity note: d̃_ic as a level is effectively absorbed by worker
  + contract FEs (it is constant within wc and approximately spans the
  worker × contract space). The interactions with post and treated×post
  are identified because they vary within wc pairs across time.
**********************************************************************/

di as text "------------------------------------------------------------"
di as text "[3.9] Parametric distance DiD"
di as text "------------------------------------------------------------"

* Build RHS for the parametric het DiD
local rhs9 "1.treated_BC_3#1.post_3"
if `did_post_main'==1 local rhs9 "`rhs9' 1.post_3"

* Triple interaction: d̃_ic × treated × post (the key het slope)
local rhs9 "`rhs9' c.dtilde_ic#1.treated_BC_3#1.post_3"

* Common post × d̃_ic slope (identifies whether exposure predicts
* pre/post change in control group)
if `did_postx'==1 local rhs9 "`rhs9' c.dtilde_ic#1.post_3"

di as text "[3.9] RHS: `rhs9'"
di as text "[3.9] ABS: `ABS_DID_3'"

reghdfe $Y_HET ///
    `rhs9' ///
    $X_3 $LIFECYCLE_3 ///
    if inwin_3==1 & !missing(dtilde_ic), ///
    absorb(`ABS_DID_3') ///
    $VCE_3
estimates store het_dist_param
di as text "[3.9] Parametric distance DiD: DONE"


/**********************************************************************
  3.10) NONPARAMETRIC DISTANCE BINS DiD

  Tests whether the treatment effect differs across DBAR_BINS_N groups
  of workers sorted by pre-reform floor exposure (d̄_ic quartiles by
  default).

  Model:
    Y = β_0 (treated × post)                             ← bin 1 baseline
      + [β_1 (post)]                       if DID_POST_MAIN_3=1
      + Σ_{k=2}^{K} β_k (bin_k × treated × post)       ← bin diffs from bin 1
      + [Σ_{k=2}^{K} γ_k (bin_k × post)]  if DID_POSTX_3=1
      + controls + absorbed FEs + ε

  Bin 1 = lowest d̄_ic (closest to floor, highest risk). The baseline
  treatment effect β_0 applies to bin 1 workers. β_k = β_0 + Δ_k
  where Δ_k is the incremental difference for bin k vs bin 1.

  The prediction from theory: β_0 should be the most negative (largest
  productivity decline) and β_k should increase toward zero (or become
  positive) as k increases (workers further from the floor).

  This specification is more flexible than the parametric version but
  requires sufficient sample in each bin for precise estimation.
**********************************************************************/

di as text "------------------------------------------------------------"
di as text "[3.10] Nonparametric distance bins DiD"
di as text "------------------------------------------------------------"

di as text "[3.10] Bin distribution in inwin_3:"
tab dbar_bin if inwin_3==1

local rhs10 "1.treated_BC_3#1.post_3"
if `did_post_main'==1 local rhs10 "`rhs10' 1.post_3"

* Bin × treated × post interactions (bin 1 is absorbed as baseline)
local rhs10 "`rhs10' i.dbar_bin#1.treated_BC_3#1.post_3"

* Common bin × post interactions
if `did_postx'==1 local rhs10 "`rhs10' i.dbar_bin#1.post_3"

di as text "[3.10] RHS: `rhs10'"
di as text "[3.10] ABS: `ABS_DID_3'"

reghdfe $Y_HET ///
    `rhs10' ///
    $X_3 $LIFECYCLE_3 ///
    if inwin_3==1 & !missing(dbar_bin), ///
    absorb(`ABS_DID_3') ///
    $VCE_3
estimates store het_dist_bins
di as text "[3.10] Distance bins DiD: DONE"


/**********************************************************************
  3.11) VOLATILITY HETEROGENEITY DiD

  Tests whether the treatment effect varies with earnings VOLATILITY
  (σ_resid and σ_weather), separate from the mean distance.

  Two volatility dimensions:
    σ_resid_wc: idiosyncratic residual volatility (worker-contract).
                Workers with high σ_resid have unpredictable daily
                earnings and face binding risk even with high expected
                earnings.
    σ_weather_i: worker-level average weather volatility.
                 Workers on contracts with high weather variability
                 face more earnings uncertainty from external shocks.

  Both are standardized to mean 0 SD 1 within inwin_3 for interpretable
  coefficients (σ_resid_z, σ_weather_z).

  Model:
    Y = β_0 (treated × post)
      + [β_1 (post)]                        if DID_POST_MAIN_3=1
      + β_2 (σ_resid_z × treated × post)
      + β_3 (σ_weather_z × treated × post)
      + [β_4 (σ_resid_z × post)]            if DID_POSTX_3=1
      + [β_5 (σ_weather_z × post)]          if DID_POSTX_3=1
      + controls + absorbed FEs + ε

  Identification note for σ_weather: σ_weather_c is a contract-level
  scalar. With contract FE absorbed, only the interaction with post
  (and treated × post) is identified — from contracts with different
  weather volatility showing different pre/post changes. Identification
  may be weak with ~25 BC contracts. σ_weather_i (used here) varies at
  the worker level and is more stable.
**********************************************************************/

di as text "------------------------------------------------------------"
di as text "[3.11] Volatility heterogeneity DiD"
di as text "------------------------------------------------------------"

* ---- Standardize volatility measures within inwin_3 ----
cap drop sigma_resid_z sigma_weather_z

quietly summ sigma_resid_wc if inwin_3==1 & !missing(sigma_resid_wc)
if r(N) > 0 & r(sd) > 0 {
    scalar sr_mu = r(mean)
    scalar sr_sd = r(sd)
    gen double sigma_resid_z = (sigma_resid_wc - sr_mu) / sr_sd ///
        if !missing(sigma_resid_wc)
    di as text "[3.11] σ_resid_z: mean=" %6.3f sr_mu "  SD=" %6.3f sr_sd
}
else {
    di as error "[3.11] WARN: insufficient σ_resid_wc obs. Falling back to worker-level."
    quietly summ sigma_resid_w if inwin_3==1
    scalar sr_mu = r(mean)
    scalar sr_sd = r(sd)
    gen double sigma_resid_z = (sigma_resid_w - sr_mu) / sr_sd ///
        if !missing(sigma_resid_w)
}

quietly summ sigma_weather_i if inwin_3==1 & !missing(sigma_weather_i)
if r(N) > 0 & r(sd) > 0 {
    scalar sw_mu = r(mean)
    scalar sw_sd = r(sd)
    gen double sigma_weather_z = (sigma_weather_i - sw_mu) / sw_sd ///
        if !missing(sigma_weather_i)
    di as text "[3.11] σ_weather_z: mean=" %6.3f sw_mu "  SD=" %6.3f sw_sd
}
else {
    di as error "[3.11] WARN: insufficient σ_weather_i obs in inwin_3."
    gen double sigma_weather_z = .
}

* ---- Volatility het DiD regression ----
local rhs11 "1.treated_BC_3#1.post_3"
if `did_post_main'==1 local rhs11 "`rhs11' 1.post_3"
local rhs11 "`rhs11' c.sigma_resid_z#1.treated_BC_3#1.post_3"
local rhs11 "`rhs11' c.sigma_weather_z#1.treated_BC_3#1.post_3"
if `did_postx'==1 {
    local rhs11 "`rhs11' c.sigma_resid_z#1.post_3"
    local rhs11 "`rhs11' c.sigma_weather_z#1.post_3"
}

di as text "[3.11] RHS: `rhs11'"
di as text "[3.11] ABS: `ABS_DID_3'"

count if inwin_3==1 & !missing(sigma_resid_z) & !missing(sigma_weather_z)
di as text "[3.11] Estimation N = " r(N)

reghdfe $Y_HET ///
    `rhs11' ///
    $X_3 $LIFECYCLE_3 ///
    if inwin_3==1 & !missing(sigma_resid_z) & !missing(sigma_weather_z), ///
    absorb(`ABS_DID_3') ///
    $VCE_3
estimates store het_volatility
di as text "[3.11] Volatility DiD: DONE"


/**********************************************************************
  3.12) DESCRIPTIVE TABLE: d̄_ic BINS × PRE/POST (BC ONLY)

  Tabulates mean piece rate, productivity, and hours by distance bin
  and pre/post period within BC. Useful for:
  - Confirming that bin assignment correlates with raw outcome patterns.
  - Checking whether pre-period differences across bins are plausible.
  - Verifying that the DiD results have intuitive raw counterparts.

  Note: this uses the Stata 17+ table syntax. In Stata 17+, the
  table command produces a formatted output to the results window
  and log. For paper-ready tables, use collect export or putexcel.
**********************************************************************/

di as text "------------------------------------------------------------"
di as text "[3.12] Descriptive table: d̄_ic bins × pre/post (BC)"
di as text "------------------------------------------------------------"

table dbar_bin post_3 ///
    if (year==$TY_3) & (province_s=="BC") ///
       & inrange(event_week_3, -$WPRE_3, $WPOST_3) ///
       & !missing(dbar_bin), ///
    statistic(mean piece) statistic(mean prod) statistic(mean tot_hrs_day) ///
    nformat(%9.3f)


/**********************************************************************
  3.13) KEY ESTIMATES PRINT

  Extracts and prints lincom results for the key het DiD coefficients
  from each stored estimate. Uses _trylincom for graceful handling of
  omitted/unavailable coefficients.

  Reading the output:
  - het_dist_param: β on d̃_ic × treated × post. Positive = workers
    FURTHER from the floor (higher d̃_ic) show a LARGER response.
    Under the theory, we expect the opposite (β < 0): workers closer
    to the floor (lower d̃_ic) should respond more.
  - het_dist_bins: β for each bin vs bin 1 (closest to floor). If
    theory holds, these should be positive and increasing (bins 2-4
    respond less negatively than bin 1).
  - het_volatility: β on σ_resid_z × treated × post and
    σ_weather_z × treated × post. Positive = more volatile workers
    show a larger treatment response.
**********************************************************************/

di as text "------------------------------------------------------------"
di as text "[3.13] Key estimates"
di as text "------------------------------------------------------------"

cap program drop _trylincom
program define _trylincom
    syntax , COEF(string) LABEL(string)
    capture noisily lincom `coef'
    if _rc==0 di as text "   -> `label'"
    else      di as text "   -> `label' (not in model / omitted / name mismatch)"
end

di as text "============================================================"
di as text "SECTION 3: KEY RESULTS"
di as text "============================================================"

* (A) Parametric distance het DiD
capture quietly estimates restore het_dist_param
if _rc {
    di as error "het_dist_param not found — Section 3.9 may have failed."
}
else {
    di as text "---- (A) Parametric distance: het_dist_param ----"
    _trylincom, coef("1.treated_BC_3#1.post_3") ///
        label("Treat×Post (d̃_ic=0, i.e., at mean exposure)")
    _trylincom, coef("c.dtilde_ic#1.treated_BC_3#1.post_3") ///
        label("Treat×Post × d̃_ic slope (1 SD further from floor)")
    if `did_postx'==1 {
        _trylincom, coef("c.dtilde_ic#1.post_3") ///
            label("Post × d̃_ic slope (control group — pre-trend check)")
    }
}

di as text ""

* (B) Nonparametric distance bins het DiD
capture quietly estimates restore het_dist_bins
if _rc {
    di as error "het_dist_bins not found — Section 3.10 may have failed."
}
else {
    di as text "---- (B) Distance bins: het_dist_bins ----"
    _trylincom, coef("1.treated_BC_3#1.post_3") ///
        label("Treat×Post in bin 1 (closest to floor, highest risk)")
    forvalues k = 2/$DBAR_BINS_N {
        _trylincom, coef("`k'.dbar_bin#1.treated_BC_3#1.post_3") ///
            label("Bin `k' vs bin 1 differential effect")
    }
}

di as text ""

* (C) Volatility het DiD
capture quietly estimates restore het_volatility
if _rc {
    di as error "het_volatility not found — Section 3.11 may have failed."
}
else {
    di as text "---- (C) Volatility: het_volatility ----"
    _trylincom, coef("1.treated_BC_3#1.post_3") ///
        label("Treat×Post (at mean σ_resid and σ_weather)")
    _trylincom, coef("c.sigma_resid_z#1.treated_BC_3#1.post_3") ///
        label("Treat×Post × σ_resid_z (1 SD more idiosyncratic volatility)")
    _trylincom, coef("c.sigma_weather_z#1.treated_BC_3#1.post_3") ///
        label("Treat×Post × σ_weather_z (1 SD more weather volatility)")
}

di as text "============================================================"


/**********************************************************************
  3.14) BINDING PROBABILITY (OPTIONAL)

  Construct π̂_ict = Pr(lr_ict < log(M_pre)) using the empirical CDF
  of pre-period residuals (see framework Section 10).

  For each observation, the argument to the CDF is:
      z_ict = log(M_pre) − l r̂_ict = −d_ict

  A more negative z_ict (large positive d_ict, far above floor) maps
  to a lower probability. A positive z_ict (below floor on average)
  maps to a probability above 0.5.

  DEFAULT: pooled empirical CDF across all pre-period residuals.
    F̂(z) = share of pre-period residuals ≤ z.
  Advantage: stable with moderate sample sizes.
  Disadvantage: assumes homoskedastic residuals across workers.

  Set DO_BINDING_PROB=1 in the control panel (Section 3.1) to run.
  Skipped by default to keep the baseline run fast.
**********************************************************************/

di as text "------------------------------------------------------------"
di as text "[3.14] Binding probability (optional)"
di as text "------------------------------------------------------------"

* ---- Helper program: Mata empirical CDF ----
* Defining the Mata block inside a 'program define ... end' wrapper
* isolates it from Stata's if { } brace counter. This avoids a known
* Stata batch-mode issue where mata: { } inside an if { } block causes
* the if-block's closing } to be consumed by the Mata parser, producing
* "} is not a valid command name" (r(199)) or Mata parse errors.
*
* The program is always DEFINED here (unconditionally), but is only
* CALLED when DO_BINDING_PROB=1 (below). When not called, it costs
* nothing. The program is also dropped + redefined on rerun to ensure
* the latest version is always used.
*
* Strategy inside the program:
*   (a) Load all resid_earn values and the use_dbar flag into Mata.
*   (b) Select pre-period residuals (use_dbar==1, non-missing); sort once.
*   (c) For each obs in sample_3 with a valid z_ict, evaluate
*       F_hat(z) = #{resid_sorted <= z} / N_r
*       using Mata's vectorised compiled :<=  operator — vastly faster
*       than a Stata count-in-loop (O(N*N_r) interpreter passes).
*   (d) Write all results back to Stata in one st_store() call.
*
* z_ict = -d_ict = log(M_pre) - lr_hat (the CDF argument).
* pi_hat = F_hat(z_ict) = Pr(resid <= z_ict) = Pr(lr_ict < log(M_pre)).

cap program drop _run_pi_cdf
program define _run_pi_cdf
    mata: {

        // (a) Load residuals and pre-period flag
        resid_all    = st_data(., "resid_earn")
        dbar_flag    = st_data(., "use_dbar")

        // (b) Select pre-period, non-missing residuals; sort once
        nonmiss_r    = (resid_all :< .)          // 1 if not Mata missing
        resid_keep   = (dbar_flag :== 1) :& nonmiss_r
        resid_sorted = sort(select(resid_all, resid_keep), 1)
        N_r          = rows(resid_sorted)

        if (N_r == 0) {
            errprintf("[3.14 Mata] ERROR: No pre-period residuals found. " +
                      "Check use_dbar and resid_earn. pi_hat not computed.\n")
        }
        else {

            printf("[3.14 Mata] N_r = %g pre-period residuals in CDF.\n", N_r)

            // (c) Load z_ict, sample flag, and pi_hat column
            z_all   = st_data(., "z_ict")
            s3_all  = st_data(., "sample_3")
            pi_all  = st_data(., "pi_hat")    // initialised to missing
            N_obs   = rows(z_all)

            // Evaluate F_hat(z_ict) for each valid obs.
            // sum(resid_sorted :<= z_all[i]) is a single compiled
            // vector operation: compare N_r values, sum the 0/1 result.
            for (i = 1; i <= N_obs; i++) {
                if (s3_all[i] == 1 & z_all[i] < .) {
                    pi_all[i] = sum(resid_sorted :<= z_all[i]) / N_r
                }
            }

            // (d) Write results back to Stata
            st_store(., "pi_hat", pi_all)
            printf("[3.14 Mata] pi_hat stored for sample_3 obs with valid z_ict.\n")
        }
    }
end    // end program _run_pi_cdf

* ---- Main execution: conditional on DO_BINDING_PROB ----
if $DO_BINDING_PROB==1 {

    di as text "[3.14] Constructing pi_hat via Mata empirical CDF..."

    * z_ict = -d_ict = log(M_pre) - lr_hat (the CDF argument).
    * pi_hat will be filled by the Mata program below.
    cap drop z_ict pi_hat
    gen double z_ict  = -d_ict if sample_3==1 & !missing(d_ict)
    gen double pi_hat = .

    _run_pi_cdf

    di as text "[3.14] pi_hat summary (sample_3):"
    summ pi_hat if sample_3==1, detail

    di as text "[3.14] Mean pi_hat by province and pre/post:"
    tabstat pi_hat if sample_3==1, by(post_3) s(mean sd) nototal

    * Sanity checks:
    * - pi_hat should be in [0, 1] by construction
    * - Workers with d_ict >> 0 (far above floor) should have pi_hat ~ 0
    * - Workers with d_ict < 0 (below floor on average) should have pi_hat > 0.5
    count if (pi_hat < 0 | pi_hat > 1) & sample_3==1
    if r(N) > 0 di as error "[3.14] WARN: " r(N) " obs with pi_hat outside [0,1]."

    di as text "[3.14] Correlation between dbar_ic and -mean(pi_hat|wc) (should be ~ -1):"
    bysort wc_id: egen double pi_hat_mean_wc = mean(pi_hat) if sample_3==1
    corr dbar_ic pi_hat_mean_wc if sample_3==1 & !missing(dbar_ic)
    cap drop pi_hat_mean_wc

}
else {
    * Create empty variables so the Section 3.0 drop list does not error on rerun.
    cap drop z_ict pi_hat
    gen double z_ict  = .
    gen double pi_hat = .
    di as text "[3.14] Skipped (DO_BINDING_PROB=0). Set to 1 to compute."
}


/**********************************************************************
  LOG: close
**********************************************************************/

di as text ""
di as text "[03] Finished: $S_DATE $S_TIME"
log close
