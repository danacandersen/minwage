# Planter Minimum Wage Project — Notes

## Research Question

Does a minimum wage increase affect tree planter worker outcomes? Specifically, the analysis exploits
BC's June 1 minimum wage cutoff to study effects on incentive pay (piece-rate × productivity),
productivity (trees/hour), piece rates, and hours worked.

---

## Data

**File:** `planter_mw.csv` (from `~/Downloads/`; built by Liyuan from payroll + Env Canada weather)

**Unit of observation:** Worker × contract × day

**N = 116,742 obs**, years 2016–2019, provinces BC (49%), ON (37%), AB (14%)

**Key variables:**

| Variable | Description |
|---|---|
| `name` | Worker identifier (string) |
| `contract` | Contract identifier (string) |
| `date`, `year`, `month` | Date fields (date imported as string YMD, converted to Stata daily) |
| `province` | BC, AB, or ON |
| `piece` | Piece rate ($ per tree) |
| `prod` | Productivity (trees per hour) |
| `tot_hrs` | Total hours worked |
| `top_up` | **0/1 dummy** — whether worker received MW top-up (not dollar amount) |
| `mw` | Minimum wage applicable |
| `experience` | Cross-season experience (years) |
| `temp`, `total_precip`, `wind_spd` | Weather controls (from Env Canada) |
| `earnings`, `hrly_earn` | Total and hourly earnings (pre-computed in raw data) |
| `lpiece`, `lhrly_earn` | Pre-computed log variables (redundant with constructed outcomes) |

**Constructed outcomes (in `00_data_prep.do`):**

| Variable | Formula |
|---|---|
| `ln_incentive` | ln(piece × prod) |
| `lprod` | ln(prod) — trees/hour |
| `ln_piece` | ln(piece) |
| `topup_amount` | max(0, (mw − piece×prod) × tot\_hrs) — constructed dollar top-up |

**Missingness (after correct import):**
- All outcomes, piece, prod, mw, experience: **0% missing**
- `temp`: 1.2% missing (1,404 obs)
- `wind_spd`: 0.6% missing (723 obs)
- `total_precip`: 8.0% missing (9,352 obs)

---

## Empirical Strategy

### Treatment

BC received a minimum wage increase effective **June 1** (year varies by analysis).
The cutoff creates a sharp within-season treatment date.

### Section 2A — Cross-province DiD + Event Study

- **Treated:** BC workers in the treated year (default: 2018)
- **Control:** AB or ON workers in the same year
- **Comparison:** BC vs. AB (or ON), before vs. after June 1
- **Regression:** `reghdfe Y i.treated_BC_2A##i.post_2A [controls]`, absorbed FEs + event-week FE
- **Event study:** Week-by-week coefficients relative to week −1 baseline
- **Preferred control: AB** (see Ontario note below)

### Section 2B — Within-BC Placebo Years

- **Treated group:** BC workers in the treated year (2018)
- **Placebo group:** BC workers in non-reform years (2016, 2017)
- **Purpose:** Test whether effects are specific to the reform year vs. generic seasonal patterns
- **Event time:** Year-specific June 1 for each year (so seasons are aligned)

### Section 3 — Heterogeneity by Pre-Determined Exposure (DDD)

See full framework in `mw_exposure_framework_final.md`.

**Key measure: d̄_ic (pre-determined exposure to the MW floor)**

Constructed as:
```
d̄_ic = α̂_i + γ̂_c + W̄_c,pre β̂ − log(M_pre)
```
where:
- `α̂_i` = worker FE from a province-specific earnings model estimated on **pre-period untreated data**
- `γ̂_c` = contract FE from the same model
- `W̄_c,pre β̂` = fitted weather index from pre-period average weather for contract c
- `log(M_pre)` = log of the pre-reform minimum wage (treated year)

`d̄_ic` is the expected log gap between predicted incentive earnings and the MW floor,
**computed entirely from pre-reform information**. It is strictly pre-determined with
respect to the reform and serves as a proxy for the probability of the floor binding.

- **High d̄_ic** → worker/contract combination predicted far above the floor → low binding risk
- **Low d̄_ic** → worker/contract close to or below predicted floor → high binding risk

**Province-specific estimation:**
- BC: estimated on pre-period BC obs only (N=5,582, R²_within=0.228)
- AB: estimated on all AB obs pre+post (AB is untreated, so all years are pre-reform; only 1 of 6
  AB contracts has pre-period BC-window data, so full AB sample needed to identify all contract FEs)
- ON: not used as control (ON 2018 had its own large MW increase)

**Standardization:** d̄_ic standardized to mean 0, SD 1 within the heterogeneity sample.
Quartile bins constructed on the pre-period BC+AB distribution.

**Binding probability (Section 3.14):** `π̂_ict = F̂(−d̄_ict)` where F̂ is the empirical CDF
of pre-period residual earnings. Computed via Mata for speed; stored as `pi_hat`.

---

## Key Design Choices

### Fixed Effects (global `$FE`)

| Option | Absorbed |
|---|---|
| `worker_contract` *(default)* | `worker_id`, `contract_id` |
| `wc_only` | `wc_id` (worker × contract) |
| `worker_only` | `worker_id` |
| `contract_only` | `contract_id` |

### Clustering
- Default: `worker_id`
- Alternatives: `contract_id`, `wc_id`

### Sample Floors (toggleable)
- Min hours per day: 2
- Min productivity: 50 trees/hour
- Min earnings: $30
- Only observations with positive piece rate and production

### Controls (all toggleable)
- Weather: temp, temp², wind, wind², precip indicators and quadratics
- Hours: `tot_hrs_day`
- Day flags: `multi_contract_day`
- Lifecycle: contract ramp (days since contract start, quadratic), season ramp (days worked in season, quadratic), ramp × cumulative seasons interaction
- Note: experience dropped from Section 3 (collinear with worker FE when FEs are absorbed)

### Auto Event-Window Selection (TAU rule)
- Finds the widest contiguous window around June 1 where both groups have ≥ TAU × core-week observations
- TAU = 0.25 (25% of average core-week count)
- Core weeks: 0 and 1 (just after cutoff)
- Max window: ±30 weeks
- Falls back to manual window if auto-selection fails

---

## Key Data Facts (from 01_data_checks.do, 2026-03-04)

### Coverage
| Province | 2016 | 2017 | 2018 | 2019 |
|---|---|---|---|---|
| BC obs | 14,618 | 14,587 | 13,946 | 13,618 |
| AB obs | 3,806 | 3,840 | 4,144 | 4,993 |
| ON obs | 11,822 | 12,771 | 10,857 | 7,740 |

| Province | Workers/yr | Contracts/yr |
|---|---|---|
| BC | ~377–427 | ~23–27 |
| AB | ~122–157 | 6 (fixed) |
| ON | ~258–424 | ~11–15 |

### top_up Validity
`top_up` is a 0/1 dummy for whether MW floor was binding on a given day.
Constructed binding condition: `piece*prod < mw` (hourly comparison).
Cross-tab shows near-perfect alignment: only 37 misclassified obs out of 116,742.

### MW Binding Rates by Province × Year
| Province | 2016 | 2017 | 2018 | 2019 |
|---|---|---|---|---|
| BC | 2.1% | 3.1% | 3.4% | **5.9%** |
| AB | 4.0% | 3.7% | 4.3% | 3.4% |
| ON | 13.2% | 17.6% | **23.9%** | 13.4% |

- **BC 2019** jump to 5.9% reflects the second BC MW increase.
- **ON 2018** spikes to 23.9% because ON had its own large MW increase ($11.40→$14.00) in 2018.
  **Ontario is not a clean control for 2018.** Use AB as the preferred control.
- **AB** is flat and low throughout — cleanest control province.

### Imputed Top-Up Amounts (when MW binding)
- Mean: $32.29/day; Median: $27.72/day; P75: $49.50; Max: $112/day
- AB binding obs: mean earn_hourly ≈ $9.23/hr vs MW ≈ $13.08/hr
- BC binding obs: mean earn_hourly ≈ $8.26/hr vs MW ≈ $11.88/hr

---

## Section 3 Key Results (from 03_heterogeneity.do, 2026-03-06)

**Sample:** BC vs. AB, 2018, in-window (N=11,606 obs, inwin_3). Outcome: `ln_incentive`.
Worker FE + Contract FE, clustered by `worker_id`.

### d̄_ic descriptives (Section 3.8)
- Mean: 0.766, SD: 0.387 (in DiD sample before standardization)
- Pre-period BC+AB sample: N≈9,181 obs used to construct quartile bins
- Distribution plotted in `output/dbar_hist.png`

### 3.9 Parametric DDD
| Term | Coeff | p-value |
|---|---|---|
| treat × post | −0.155 | 0.004 |
| d̃_ic × treat × post | −0.143 | 0.048 |
| d̃_ic × post (AB placebo) | — | p=0.317 ✅ |

Interpretation: a one-SD increase in d̄_ic (worker further from floor) is associated with
a 14.3 pp larger reduction in incentive pay post-reform in BC relative to AB.

### 3.10 Quartile Bins DDD
| Bin | Coeff vs. Bin 1 | p-value |
|---|---|---|
| Bin 1 (closest to floor) | −0.008 | 0.861 |
| Bin 2 vs. Bin 1 | −0.100 | 0.024 |
| Bin 3 vs. Bin 1 | −0.145 | 0.028 |
| Bin 4 (furthest from floor) vs. Bin 1 | −0.227 | 0.014 |

Pattern: **monotonically increasing** productivity decline as d̄_ic increases.
Workers furthest from the floor show the largest negative response — paradoxical relative
to the naïve binding-probability channel.

### 3.11 Volatility Heterogeneity
| Term | Coeff | p-value |
|---|---|---|
| σ_resid_z × treat × post | +0.069 | 0.009 |
| σ_weather_z × treat × post | −0.138 | 0.384 |

Workers with higher idiosyncratic earnings volatility show a *smaller* productivity
decline — consistent with the idea that volatile workers face higher effective binding risk
and may substitute toward hours or effort on high-productivity days.

### Economic Interpretation (working hypotheses)
1. **Income effect / backward-bending labor supply** (leading candidate): The piece-rate
   increase post-reform benefits high-d̄_ic workers most (they plant the most trees, earn
   more per tree). If income effects are strong, these workers reduce effort to hit a target
   income, producing larger output declines.
2. **Threat effects stabilizing Bin 1**: Low-d̄_ic workers face genuine binding risk;
   threat of firing or soft pressure from foremen may keep effort near the minimum,
   attenuating any response.
3. **Firm plot reallocation**: Firms may assign high-ability workers to harder plots
   post-reform to manage costs, mechanically lowering observed productivity.
4. **Mean reversion in α̂_i** (key validity threat): If high estimated α̂_i partly reflects
   transient luck rather than ability, d̄_ic is mismeasured → effects attributed to
   "high exposure" may partly be mean reversion. Check: event study by bin pre-trends.

---

## Bugs Found and Fixed

### 1. `destring` decimal point corruption (CRITICAL — fixed 2026-03-04)
**Bug:** `destring` was called with `ignore("NA" "." ...)`. The `"."` option strips the
literal period character from all string values before conversion. So `"14.7778"` became
`"147778"`, and `"14.7777777777778"` became `"147777777777778"` ≈ **1.48e+14**.
This corrupted `temp` and `wind_spd` entirely (88%+ appeared missing after a sentinel
replacement attempt).

**Root cause:** `temp` and `wind_spd` are imported as strings (because some rows contain
"NA"), so `destring` runs on them. The `"."` in the ignore list corrupts decimal values.

**Fix:** Removed `"."` from the `ignore()` list. Now use `ignore("NA" "null" "NULL")` only.

**Lesson:** Never include `"."` in `destring ignore()` for variables with decimal values.
The raw CSV data (from Env Canada via Liyuan) is clean — this was entirely our bug.

### 2. Space in project path breaks Stata batch mode (fixed 2026-03-04)
**Bug:** `run_analysis.py` passed the full do-file path to Stata's `-b do` argument.
Stata truncates paths at spaces, so `Research New/01_data_checks.do` became `Research.do`.
**Fix:** Pass only the filename (not full path); `cwd=project_dir` in subprocess handles the rest.

### 3. Python 3.9 type hint syntax (fixed 2026-03-04)
**Bug:** `str | None` union type hint requires Python 3.10+. System Python is 3.9.6.
**Fix:** Changed to untyped `def find_stata(override=None)`.

### 4. `top_up` is a dummy, not a dollar amount (clarified 2026-03-04)
**Clarification:** `top_up` in the raw data is 0/1, not a dollar amount.
Imputed top-up amount = `max(0, (mw - piece*prod) * tot_hrs)`. Fixed in `01_data_checks.do`.

### 5. Mata `{ }` inside Stata `if { }` block corrupts brace counting (fixed 2026-03-06)
**Bug:** Embedding `mata: { ... }` inside a Stata `if $TOGGLE==1 { ... }` block in batch
mode (`-b` flag) causes `r(199)` ("} is not a valid command name"). Stata's brace counter
for the outer `if { }` block is confused: the Mata closing `}` is consumed, leaving the
outer if-block's closing `}` as a stray command.

**First fix attempt:** Moved `mata: { }` outside the if-block, using `mata` ... `end`
delimiter syntax with a condition check (`st_global()`) inside Mata. This produced
`r(3000)` ("unexpected end of line") — `mata` + `end` is unreliable in Stata `-b` mode.

**Definitive fix:** Wrap `mata: { }` inside a `program define ... end` block. The program
definition fully isolates Mata's brace counting from any outer Stata flow-control context.
Define the program unconditionally; call it conditionally.

```stata
cap program drop _run_pi_cdf
program define _run_pi_cdf
    mata: {
        // ... Mata code here ...
    }
end  // brace counting is isolated inside the program

if $DO_BINDING_PROB==1 {
    _run_pi_cdf  // safe — no nested mata: { } here
}
```

**Lesson:** Never embed `mata: { }` or `mata` ... `end` inside a Stata `if { }` block in
batch mode. Always use the `program define` wrapper pattern.

---

## Do-File Structure

**Legacy monolith:** `dana_minwage_2026_feb.do` (do not edit — reference only)

**Modular files:**

| File | Purpose | Status |
|---|---|---|
| `00_data_prep.do` | Import CSV → `planter_mw_prepped.dta` | ✅ Runs clean |
| `01_data_checks.do` | Coverage, missingness, top-up validity, weather, experience | ✅ Runs clean |
| `02_main_analysis.do` | 2A + 2B DiD, auto event-window, coefplot → `output/` | ✅ Runs clean |
| `03_heterogeneity.do` | d̄_ic construction, ability, DDD by bin/parametric/volatility | ✅ Runs clean |
| `run_analysis.py` | Python master runner | ✅ Working |

**Supporting documents:**
- `binding_risk_measure.md` — initial theoretical framework notes
- `mw_exposure_framework_final.md` — authoritative framework: d̄_ic definition, province
  estimation choices, M_pre design, standardization, and binding probability

**How to run:**
```bash
python3 run_analysis.py              # all steps
python3 run_analysis.py --step 02   # single step
python3 run_analysis.py --from 01   # from step 01 onward
```
Say "run 02" to Claude and it will execute, read the log, and report findings.

**Output locations:**
- Logs: `logs/NN_stepname_YYYYMMDD.log`
- Figures: `output/event_*.png`, `output/dbar_hist.png`
- Prepared data: `planter_mw_prepped.dta` (gitignored — rebuilt by 00)
- Worktree data symlink: `planter_mw_prepped.dta` → `../../../planter_mw_prepped.dta`
  (symlink created manually; not tracked by git)

---

## Next Steps

1. Run event study split by d̄_ic quartile — check for pre-existing trends by bin (key
   validity check for the heterogeneity result)
2. Assess mean reversion threat: regress Δα̂_i (post − pre within-worker FE change) on d̄_ic;
   if large, heterogeneity result is contaminated
3. Consider robustness: use weather-only version of d̄_ic (drop α̂_i) to reduce endogeneity
   in the exposure measure
4. Check 2019 replication: run Section 3 for 2019 BC reform to test heterogeneity stability
5. Refine figures for paper quality (currently exploratory/monochrome)
6. Draft Section 3 of paper using working hypotheses above

---

## Open Questions

### Heterogeneity paradox
The DDD shows the *largest* productivity declines among workers furthest from the MW floor
(highest d̄_ic / highest ability). This is the opposite of what a naïve binding-probability
model would predict. Leading explanations: income effects from piece-rate increases,
threat effects stabilizing Bin 1, firm plot reallocation, or mean reversion in α̂_i.
Pre-trend check by bin is the most important next diagnostic.

### AB estimation sample for d̄_ic
Only 1 of 6 AB contracts (B-18-0408, 331 obs) falls within the pre-period event window
for the 2018 treated year. The other 5 AB contracts are post-only. To identify all 6 AB
contract FEs, the AB earnings model uses the full AB sample (pre+post, all years). This is
valid because AB is never treated, so all AB observations are effectively pre-reform in
the counterfactual sense. d̄_ic for post-only AB wc pairs uses their post-period `lr_hat`
fitted values (Section 3.6 fallback), which is internally consistent.

---

## Software Requirements

- Stata 19 MP (`/Applications/Stata/StataMP.app/Contents/MacOS/StataMP`)
- Packages: `reghdfe`, `ftools`, `coefplot` (auto-installed by `00_data_prep.do`)
- Python 3.9+ for `run_analysis.py`
- Raw data: `~/Downloads/planter_mw.csv`
