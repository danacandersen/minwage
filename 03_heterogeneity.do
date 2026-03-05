/**********************************************************************
  03_HETEROGENEITY.DO
  Planter Minimum Wage Project

  PURPOSE:
    Section 3 heterogeneity analysis.
    Constructs ability, weather shock, contract difficulty, and
    top-up probability (phat), then runs heterogeneity DiD regressions.

  SECTIONS:
    3.0)  Clean rerun state
    3.1)  Settings / control panel (EDIT THESE)
    3.2)  Build controls + floor strings
    3.3)  Sample + event time
    3.4)  Pre-treatment ability (ability_z)
    3.5)  Daily weather shock (wx_shock)
    3.6)  Contract difficulty (contract_diff_z)
    3.7)  Top-up risk model -> phat_topup
    3.8)  Phat moments within window + distribution diagnostics
    3.9)  Parametric phat heterogeneity DiD
    3.10) Nonparametric phat bins DiD
    3.11) Contract difficulty + weather shock heterogeneity DiD
    3.12) Ability bins heterogeneity DiD
    3.13) Descriptive table: ability bins × pre/post (BC only)
    3.14) Key estimates print (lincom)

  OUTPUT:
    logs/03_heterogeneity_YYYYMMDD.log

  REQUIRES: 00_data_prep.do to have been run first.

  KEY OPEN QUESTION:
    phat_topup is constructed from ability_z + weather + controls.
    Because ability_z enters both the residualization (Section 3.4)
    and the phat model (Section 3.7), the two may be collinear.
    This makes it hard to separately identify ability vs MW-binding
    channels in parametric heterogeneity regressions (Section 3.9).
**********************************************************************/

clear all
set more off
version 19


/**********************************************************************
  LOG
**********************************************************************/

cap mkdir "logs"

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
**********************************************************************/

cap estimates clear

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
    PHAT_VAR_LEVEL PHAT_BINS_N PHAT_SD_BINS ///
    DID_CONTRACTFE_3 DID_WEEKFE_3 DID_POST_MAIN_3 DID_POSTX_3 ///
    X_3 LIFECYCLE_3 SAMPLE_OK_3 {
    cap macro drop `m'
}


/**********************************************************************
  3.1) SETTINGS — EDIT THESE
**********************************************************************/

di as text "============================================================"
di as text "[3.1] Settings"
di as text "============================================================"

* Treatment year and control province
global TY_3      2018
global CONTROL_3 "AB"

* Outcome for heterogeneity regressions
global Y_HET "ln_incentive"

* Clustering
global CLUSTER_3 "worker_id"
global VCE_3     "vce(cluster $CLUSTER_3)"

* Control toggles (Section 3 may differ from Section 2)
global USE_WEATHER_3   1
global USE_EXPER_3     1
global USE_HOURS_3     0    // hours typically excluded from het specs
global USE_DAYFLAGS_3  0
global USE_LIFECYCLE_3 1

* Floor toggles
global USE_FLOORS_3   0    // typically no floors in het analysis
global MIN_HRS_3      2
global MIN_EARN_3    30
global USE_MIN_PROD_3 1
global MIN_PROD_3    50

* Event window for Section 3 (fixed, not auto-selected)
global WPRE_3  3
global WPOST_3 8

* Contract difficulty timing:
* 0 = PRE only (AB+BC)
* 1 = AB pre+post allowed; BC pre only   [recommended: avoids post contamination in BC]
* 2 = POST allowed for all (AB+BC)
global CD_POST_MODE 1

* Week FE in contract difficulty residualization
global CD_TIMEFE 0

* Phat moments level: "worker" or "wc" (worker×contract)
global PHAT_VAR_LEVEL "wc"

* Phat bins for nonparametric heterogeneity
global PHAT_BINS_N  4
global PHAT_SD_BINS 0    // 1 = also run sd-bins model

* ============================================================
* GLOBAL SPEC TOGGLES (apply to ALL Section 3 DiD regressions)
* ============================================================
* Include contract FE in all DiD absorb lists?
global DID_CONTRACTFE_3 1

* Include event-week FE in all DiD absorb lists?
global DID_WEEKFE_3 1

* Include post_3 main effect explicitly (or rely on worker FE for level)?
global DID_POST_MAIN_3 1

* Include common post × heterogeneity covariate terms (not just treated × post × X)?
global DID_POSTX_3 1


/**********************************************************************
  3.2) BUILD CONTROLS + FLOOR STRINGS
**********************************************************************/

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

* ---- Safe numeric locals from globals (avoids "==1 invalid name" error) ----
* reghdfe/if conditions cannot use $GLOBAL==1 syntax directly in absorb lists.
local did_contractfe = real(trim("$DID_CONTRACTFE_3"))
local did_weekfe     = real(trim("$DID_WEEKFE_3"))
local did_post_main  = real(trim("$DID_POST_MAIN_3"))
local did_postx      = real(trim("$DID_POSTX_3"))
local cdmode         = real(trim("$CD_POST_MODE"))
local cd_timefe      = real(trim("$CD_TIMEFE"))

di as text "------------------------------------------------------------"
di as text "[3] TY_3:           $TY_3   CONTROL_3: $CONTROL_3   Y_HET: $Y_HET"
di as text "[3] Window:         [-$WPRE_3,$WPOST_3]"
di as text "[3] X_3:            $X_3"
di as text "[3] LIFECYCLE_3:    $LIFECYCLE_3"
di as text "[3] SAMPLE_OK_3:    $SAMPLE_OK_3"
di as text "[3] CD_POST_MODE:   $CD_POST_MODE   CD_TIMEFE: $CD_TIMEFE"
di as text "[3] PHAT_VAR_LEVEL: $PHAT_VAR_LEVEL   PHAT_BINS_N: $PHAT_BINS_N   PHAT_SD_BINS: $PHAT_SD_BINS"
di as text "[3] DID_CONTRACTFE: $DID_CONTRACTFE_3   DID_WEEKFE: $DID_WEEKFE_3"
di as text "[3] DID_POST_MAIN:  $DID_POST_MAIN_3   DID_POSTX: $DID_POSTX_3"
di as text "------------------------------------------------------------"


/**********************************************************************
  3.3) DEFINE SAMPLE + EVENT TIME
**********************************************************************/

di as text "------------------------------------------------------------"
di as text "[3.3] Sample + event time"
di as text "------------------------------------------------------------"

gen byte sample_3 = (year == $TY_3) & inlist(province_s, "BC", "$CONTROL_3")

gen byte treated_BC_3 = .
replace treated_BC_3 = (province_s == "BC") if sample_3==1

local tdate_3 = mdy(6, 1, $TY_3)
gen byte post_3        = (date >= `tdate_3')              if sample_3==1
gen int  event_week_3  = floor((date - `tdate_3') / 7)   if sample_3==1
gen int  event_week_3s = event_week_3 + 100               if sample_3==1

di as text "[3.3] Sample counts:"
tab province_s if sample_3==1
tab post_3     if sample_3==1

* ---- Unified absorb list for ALL Section 3 DiD regressions ----
* Built from local toggles to avoid $GLOBAL==1 issues.
local ABS_DID_3 "worker_id"
if `did_contractfe'==1 local ABS_DID_3 "`ABS_DID_3' contract_id"
if `did_weekfe'==1     local ABS_DID_3 "`ABS_DID_3' event_week_3s"
di as text "[3.3] ABS_DID_3: `ABS_DID_3'"


/**********************************************************************
  3.4) PRE-TREATMENT ABILITY (ability_z)

  Approach: residualize log productivity on controls + contract FE,
  using only pre-June-1 observations in the sample. Take worker-level
  mean of the residual; standardize to mean 0, SD 1.

  Interpretation: ability_z is a measure of time-invariant worker
  productivity net of contract characteristics and controls.
  Higher ability_z = more productive worker (relative to contract peers).
**********************************************************************/

di as text "------------------------------------------------------------"
di as text "[3.4] Pre-treatment ability"
di as text "------------------------------------------------------------"

cap drop ln_prod
gen double ln_prod = ln(prod) if prod > 0

* Sample for ability: pre-period only, non-missing outcomes
cap drop abil_sample
gen byte abil_sample = (sample_3==1) & (post_3==0) & !missing(ln_prod) & ($SAMPLE_OK_3)

di as text "[3.4] abil_sample N = "
count if abil_sample==1

* Residualize: absorb contract FE to get within-contract productivity residual
cap drop resid_preprod
reghdfe ln_prod $X_3 $LIFECYCLE_3 if abil_sample==1, absorb(contract_id) resid(resid_preprod)

* Worker-level mean residual = pre-treatment ability
cap drop ability_raw ability_z
bysort worker_id: egen double ability_raw = mean(resid_preprod)

* Standardize ability
quietly summ ability_raw if !missing(ability_raw)
scalar abil_mu = r(mean)
scalar abil_sd = r(sd)
gen double ability_z = (ability_raw - abil_mu) / abil_sd ///
    if !missing(ability_raw) & abil_sd > 0

di as text "[3.4] ability_z:"
summ ability_z


/**********************************************************************
  3.5) DAILY WEATHER SHOCK (wx_shock)

  Fitted values from a weather-only OLS (no worker FE), estimated on
  the full sample_3 sample. This captures expected productivity given
  weather conditions, not individual ability.

  NOTE: Using OLS (not reghdfe) intentionally — we want a single
  fitted value per day, not within-worker deviation.
**********************************************************************/

di as text "------------------------------------------------------------"
di as text "[3.5] Daily weather shock"
di as text "------------------------------------------------------------"

cap drop wx_shock
quietly reg ln_prod ///
    c.temp##c.temp c.wind_spd##c.wind_spd precip0 c.precip_pos##c.precip_pos ///
    if sample_3==1 & !missing(ln_prod), vce(robust)
predict double wx_shock if sample_3==1, xb

di as text "[3.5] wx_shock:"
summ wx_shock if sample_3==1


/**********************************************************************
  3.6) CONTRACT DIFFICULTY (contract_diff_z)

  Goal: measure inherent difficulty of each contract, net of worker
  ability and controls.

  Approach:
    - Residualize ln_prod on ability_z + controls + worker FE separately
      for AB and BC (using CD_POST_MODE to determine which period).
    - Merge residuals; take contract-level mean.
    - Standardize.

  CD_POST_MODE=1 (recommended): AB uses pre+post data (more power);
  BC uses only pre data (avoids contamination from MW response in post).
**********************************************************************/

di as text "------------------------------------------------------------"
di as text "[3.6] Contract difficulty"
di as text "------------------------------------------------------------"

cap drop cd_sample_ab cd_sample_bc cd_sample ///
         resid_cd_ab resid_cd_bc resid_cd ///
         contract_diff_raw contract_diff_z

gen byte cd_sample_ab = (year==$TY_3) & (province_s=="AB") ///
    & !missing(ln_prod) & !missing(ability_z) & ($SAMPLE_OK_3)
gen byte cd_sample_bc = (year==$TY_3) & (province_s=="BC") ///
    & !missing(ln_prod) & !missing(ability_z) & ($SAMPLE_OK_3)

* Apply timing restriction per CD_POST_MODE
if `cdmode'==0 {
    * PRE only for both provinces
    replace cd_sample_ab = cd_sample_ab & (post_3==0)
    replace cd_sample_bc = cd_sample_bc & (post_3==0)
}
else if `cdmode'==1 {
    * AB: pre+post; BC: pre only
    replace cd_sample_bc = cd_sample_bc & (post_3==0)
}
* cdmode==2: allow post for both — do nothing

gen byte cd_sample = cd_sample_ab | cd_sample_bc

* Absorb list for difficulty residualization
local ABS_CD_3 "worker_id"
if `cd_timefe'==1 local ABS_CD_3 "`ABS_CD_3' event_week_3s"

di as text "[3.6] cd_sample_ab N = "
count if cd_sample_ab==1
di as text "[3.6] cd_sample_bc N = "
count if cd_sample_bc==1

reghdfe ln_prod c.ability_z $X_3 $LIFECYCLE_3 if cd_sample_ab==1, absorb(`ABS_CD_3') resid(resid_cd_ab)
reghdfe ln_prod c.ability_z $X_3 $LIFECYCLE_3 if cd_sample_bc==1, absorb(`ABS_CD_3') resid(resid_cd_bc)

gen double resid_cd = .
replace resid_cd = resid_cd_ab if !missing(resid_cd_ab)
replace resid_cd = resid_cd_bc if !missing(resid_cd_bc)

bysort contract_id: egen double contract_diff_raw = mean(resid_cd)

quietly summ contract_diff_raw if !missing(contract_diff_raw)
scalar cd_mu = r(mean)
scalar cd_sd = r(sd)
gen double contract_diff_z = (contract_diff_raw - cd_mu) / cd_sd ///
    if !missing(contract_diff_raw) & cd_sd > 0

di as text "[3.6] Contract difficulty coverage:"
count if sample_3==1 & province_s=="AB" & !missing(contract_diff_z)
count if sample_3==1 & province_s=="BC" & !missing(contract_diff_z)
summ contract_diff_z


/**********************************************************************
  3.7) TOP-UP RISK MODEL: phat_topup

  LPM: topup_ind ~ mw + ability_z + wx_shock + controls + contract FE
  Estimated on PRE-period only.
  Predict phat for entire sample_3 (pre + post).
  Clamp to [0,1].

  phat_topup = predicted probability that a worker×day observation
  would receive a top-up (piece earnings below MW floor).
**********************************************************************/

di as text "------------------------------------------------------------"
di as text "[3.7] Top-up risk model"
di as text "------------------------------------------------------------"

cap drop topup_ind topup_pre phat_topup
gen byte topup_ind = (top_up > 0) if !missing(top_up)

gen byte topup_pre = (sample_3==1) & (post_3==0) ///
    & !missing(topup_ind) & !missing(ability_z) ///
    & !missing(mw) & !missing(wx_shock) ///
    & ($SAMPLE_OK_3)

di as text "[3.7] topup_pre N = "
count if topup_pre==1
di as text "[3.7] topup_ind mean in topup_pre sample:"
summ topup_ind if topup_pre==1

* LPM with contract FE (i.contract_id) for within-contract variation
reg topup_ind c.mw c.ability_z c.wx_shock $X_3 i.contract_id ///
    if topup_pre==1, vce(cluster $CLUSTER_3)

* Predict for full sample_3 window
predict double phat_topup if sample_3==1, xb

* Clamp to valid probability range
replace phat_topup = max(0, min(1, phat_topup)) if sample_3==1 & !missing(phat_topup)

di as text "[3.7] phat_topup in sample_3:"
summ phat_topup if sample_3==1


/**********************************************************************
  3.8) PHAT MOMENTS + ANALYSIS WINDOW FLAG (inwin_3)

  inwin_3: observations used for heterogeneity DiD regressions.
           Applies event window, floor filters, and non-missing on
           outcome + phat.

  phat moments computed at the level set by $PHAT_VAR_LEVEL:
    "worker": worker-level mean/SD of phat within window
    "wc":     worker×contract mean/SD of phat within window
**********************************************************************/

di as text "------------------------------------------------------------"
di as text "[3.8] Phat moments + window flag"
di as text "------------------------------------------------------------"

cap drop inwin_3
gen byte inwin_3 = (sample_3==1) ///
    & inrange(event_week_3, -$WPRE_3, $WPOST_3) ///
    & ($SAMPLE_OK_3) ///
    & !missing($Y_HET) ///
    & !missing(phat_topup)

di as text "[3.8] inwin_3 N = "
count if inwin_3==1

cap drop phat_mu_w phat_sd_w phat_mu_wc phat_sd_wc

if "$PHAT_VAR_LEVEL"=="worker" {
    bysort worker_id: egen double phat_mu_w = mean(cond(inwin_3==1, phat_topup, .))
    bysort worker_id: egen double phat_sd_w = sd(  cond(inwin_3==1, phat_topup, .))
}
else {
    capture confirm variable wc_id
    if _rc {
        di as error "[3.8] wc_id not found — falling back to worker-level moments"
        bysort worker_id: egen double phat_mu_w = mean(cond(inwin_3==1, phat_topup, .))
        bysort worker_id: egen double phat_sd_w = sd(  cond(inwin_3==1, phat_topup, .))
        global PHAT_VAR_LEVEL "worker"
    }
    else {
        bysort wc_id: egen double phat_mu_wc = mean(cond(inwin_3==1, phat_topup, .))
        bysort wc_id: egen double phat_sd_wc = sd(  cond(inwin_3==1, phat_topup, .))
    }
}

* Short label for phat variable names used in regressions
local mu = cond("$PHAT_VAR_LEVEL"=="worker", "phat_mu_w", "phat_mu_wc")
local sd = cond("$PHAT_VAR_LEVEL"=="worker", "phat_sd_w", "phat_sd_wc")


* --- 3.8B: Distribution diagnostics ---
di as text "------------------------------------------------------------"
di as text "[3.8B] phat distribution in window (inwin_3==1)"
di as text "------------------------------------------------------------"

if "$PHAT_VAR_LEVEL"=="worker" {
    summ phat_mu_w phat_sd_w if inwin_3==1, detail
    histogram phat_mu_w if inwin_3==1, ///
        width(.02) frequency ///
        title("phat mean (worker), in window") ///
        graphregion(color(white))
}
else {
    summ phat_mu_wc phat_sd_wc if inwin_3==1, detail
    histogram phat_mu_wc if inwin_3==1, ///
        width(.02) frequency ///
        title("phat mean (wc), in window") ///
        graphregion(color(white))
}


/**********************************************************************
  3.9) PARAMETRIC phat HETEROGENEITY DiD

  RHS: treated×post + (optional) post main + phat_mu slopes ×
       treated×post + (optional) phat_mu × post common slope.
  Purpose: test if workers with higher pre-treatment top-up risk
  show different productivity response to the MW increase.

  NOTE: collinearity concern — phat_mu correlates with ability_z.
  Interpret slopes with caution; see PROJECT_NOTES.md.
**********************************************************************/

di as text "------------------------------------------------------------"
di as text "[3.9] Parametric phat heterogeneity DiD"
di as text "------------------------------------------------------------"

local rhs9 "1.treated_BC_3#1.post_3"
if `did_post_main'==1 local rhs9 "`rhs9' 1.post_3"

* Always: treated×post×phat slopes
local rhs9 "`rhs9' c.`mu'#1.treated_BC_3#1.post_3"
local rhs9 "`rhs9' c.`sd'#1.treated_BC_3#1.post_3"

* Optional: common post×phat slopes (identifies level shift vs heterogeneous effect)
if `did_postx'==1 {
    local rhs9 "`rhs9' c.`mu'#1.post_3"
    local rhs9 "`rhs9' c.`sd'#1.post_3"
}

di as text "[3.9] RHS: `rhs9'"
di as text "[3.9] ABS: `ABS_DID_3'"

reghdfe $Y_HET ///
    `rhs9' ///
    $X_3 $LIFECYCLE_3 ///
    if inwin_3==1 & !missing(`mu') & !missing(`sd'), ///
    absorb(`ABS_DID_3') ///
    $VCE_3
estimates store het_phat_param


/**********************************************************************
  3.10) NONPARAMETRIC phat BINS DiD

  Assign observations to phat quartiles (or $PHAT_BINS_N bins).
  Run DiD with bin × treated × post interactions.
  Baseline bin = lowest phat (least MW-binding risk).

  Optional: also run sd-bins model (if PHAT_SD_BINS==1).
**********************************************************************/

di as text "------------------------------------------------------------"
di as text "[3.10] Nonparametric phat bins DiD"
di as text "------------------------------------------------------------"

cap drop phat_mu_bin phat_sd_bin

if "$PHAT_VAR_LEVEL"=="worker" {
    xtile phat_mu_bin = phat_mu_w if inwin_3==1, nq($PHAT_BINS_N)
    if $PHAT_SD_BINS==1 xtile phat_sd_bin = phat_sd_w if inwin_3==1, nq($PHAT_BINS_N)

    * Propagate bin assignment to all obs for the same worker
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

    * Propagate bin to all obs for the same wc_id
    bysort wc_id: egen byte phat_mu_bin2 = max(phat_mu_bin)
    cap drop phat_mu_bin
    rename phat_mu_bin2 phat_mu_bin

    if $PHAT_SD_BINS==1 {
        bysort wc_id: egen byte phat_sd_bin2 = max(phat_sd_bin)
        cap drop phat_sd_bin
        rename phat_sd_bin2 phat_sd_bin
    }
}

di as text "[3.10] phat_mu_bin distribution:"
tab phat_mu_bin if inwin_3==1

* mu-bins model
local rhs10 "1.treated_BC_3#1.post_3"
if `did_post_main'==1 local rhs10 "`rhs10' 1.post_3"
local rhs10 "`rhs10' i.phat_mu_bin#1.treated_BC_3#1.post_3"
if `did_postx'==1 local rhs10 "`rhs10' i.phat_mu_bin#1.post_3"

di as text "[3.10] RHS (mu bins): `rhs10'"

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

    di as text "[3.10] RHS (sd bins): `rhs10s'"

    reghdfe $Y_HET ///
        `rhs10s' ///
        $X_3 $LIFECYCLE_3 ///
        if inwin_3==1 & !missing(phat_sd_bin), ///
        absorb(`ABS_DID_3') ///
        $VCE_3
    estimates store het_phat_bins_sd
}


/**********************************************************************
  3.11) CONTRACT DIFFICULTY + WEATHER SHOCK HETEROGENEITY DiD

  NOTE: identification of contract_diff_z requires spanning contracts
  (workers appearing in multiple contracts) when DID_CONTRACTFE_3==1.
  Check coverage with count if ... above.
**********************************************************************/

di as text "------------------------------------------------------------"
di as text "[3.11] Difficulty + weather shock heterogeneity DiD"
di as text "------------------------------------------------------------"

cap drop est11
gen byte est11 = inwin_3==1 & !missing(contract_diff_z) & !missing(wx_shock)

di as text "[3.11] est11 N = "
count if est11==1

local rhs11 "1.treated_BC_3#1.post_3"
if `did_post_main'==1 local rhs11 "`rhs11' 1.post_3"
local rhs11 "`rhs11' c.contract_diff_z#1.treated_BC_3#1.post_3"
local rhs11 "`rhs11' c.wx_shock#1.treated_BC_3#1.post_3"
if `did_postx'==1 {
    local rhs11 "`rhs11' c.contract_diff_z#1.post_3"
    local rhs11 "`rhs11' c.wx_shock#1.post_3"
}

di as text "[3.11] RHS: `rhs11'"

reghdfe $Y_HET ///
    `rhs11' ///
    $X_3 $LIFECYCLE_3 ///
    if est11==1, ///
    absorb(`ABS_DID_3') ///
    $VCE_3
estimates store het_difficulty


/**********************************************************************
  3.12) ABILITY BINS HETEROGENEITY DiD

  Assign workers to ability quartiles (within sample_3).
  Run DiD with ability_bin × treated × post interactions.
  Captures heterogeneity by baseline productivity.
**********************************************************************/

di as text "------------------------------------------------------------"
di as text "[3.12] Ability bins heterogeneity DiD"
di as text "------------------------------------------------------------"

cap drop abil_bin_all
xtile abil_bin_all = ability_z if sample_3==1 & !missing(ability_z), nq(4)

di as text "[3.12] abil_bin_all distribution:"
tab abil_bin_all if sample_3==1

local rhs12 "1.treated_BC_3#1.post_3"
if `did_post_main'==1 local rhs12 "`rhs12' 1.post_3"
local rhs12 "`rhs12' i.abil_bin_all#1.treated_BC_3#1.post_3"
if `did_postx'==1 local rhs12 "`rhs12' i.abil_bin_all#1.post_3"

di as text "[3.12] RHS: `rhs12'"

reghdfe $Y_HET ///
    `rhs12' ///
    $X_3 $LIFECYCLE_3 ///
    if inwin_3==1 & !missing(ability_z) & !missing(abil_bin_all), ///
    absorb(`ABS_DID_3') ///
    $VCE_3
estimates store het_ability_bins


/**********************************************************************
  3.13) DESCRIPTIVE TABLE: ability bins × pre/post (BC only)

  Tabulates piece rate, productivity, and hours by ability quartile
  and pre/post period within BC. Useful for checking raw patterns.
**********************************************************************/

di as text "------------------------------------------------------------"
di as text "[3.13] Descriptive table: ability bins x pre/post (BC)"
di as text "------------------------------------------------------------"

table abil_bin_all post_3 ///
    if (year==$TY_3) & (province_s=="BC") ///
       & inrange(event_week_3, -$WPRE_3, $WPOST_3), ///
    statistic(mean piece) statistic(mean prod) statistic(mean tot_hrs_day) ///
    nformat(%9.3f)


/**********************************************************************
  3.14) KEY ESTIMATES PRINT (copy/paste friendly)

  Uses a helper program _trylincom to gracefully handle cases where
  a coefficient is omitted or not in the model.
**********************************************************************/

di as text "------------------------------------------------------------"
di as text "[3.14] Key estimates"
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

* (A) Parametric phat
capture quietly estimates restore het_phat_param
if _rc {
    di as error "het_phat_param not found"
}
else {
    di as text "---- (A) Parametric phat: het_phat_param ----"
    _trylincom, coef("1.treated_BC_3#1.post_3")              label("Treat×Post (baseline)")
    _trylincom, coef("c.`mu'#1.treated_BC_3#1.post_3")       label("Treat×Post × phat_mu slope")
    _trylincom, coef("c.`sd'#1.treated_BC_3#1.post_3")       label("Treat×Post × phat_sd slope")
}

di as text ""

* (B) Nonparametric phat bins
capture quietly estimates restore het_phat_bins
if _rc {
    di as error "het_phat_bins not found"
}
else {
    di as text "---- (B) phat mu-bins: het_phat_bins ----"
    _trylincom, coef("1.treated_BC_3#1.post_3") label("Treat×Post (bin 1 baseline)")
    foreach k of numlist 2/$PHAT_BINS_N {
        _trylincom, coef("`k'.phat_mu_bin#1.treated_BC_3#1.post_3") label("Bin `k' vs bin 1")
    }
}

di as text ""

* (C) Difficulty
capture quietly estimates restore het_difficulty
if _rc {
    di as error "het_difficulty not found"
}
else {
    di as text "---- (C) Difficulty: het_difficulty ----"
    _trylincom, coef("1.treated_BC_3#1.post_3")                     label("Treat×Post (baseline)")
    _trylincom, coef("c.contract_diff_z#1.treated_BC_3#1.post_3")   label("Treat×Post × contract_diff slope")
    _trylincom, coef("c.wx_shock#1.treated_BC_3#1.post_3")          label("Treat×Post × wx_shock slope")
}

di as text ""

* (D) Ability bins
capture quietly estimates restore het_ability_bins
if _rc {
    di as error "het_ability_bins not found"
}
else {
    di as text "---- (D) Ability bins: het_ability_bins ----"
    _trylincom, coef("1.treated_BC_3#1.post_3") label("Treat×Post (bin 1 baseline)")
    foreach k in 2 3 4 {
        _trylincom, coef("`k'.abil_bin_all#1.treated_BC_3#1.post_3") label("Ability bin `k' vs bin 1")
    }
}

di as text "============================================================"


/**********************************************************************
  LOG: close
**********************************************************************/

di as text ""
di as text "[03] Finished: $S_DATE $S_TIME"
log close
