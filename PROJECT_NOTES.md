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
- Experience: experience and experience²
- Hours: `tot_hrs_day`
- Day flags: `multi_contract_day`
- Lifecycle: contract ramp (days since contract start, quadratic), season ramp (days worked in season, quadratic), ramp × cumulative seasons interaction

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

## Bugs Found and Fixed (2026-03-04 session)

### 1. `destring` decimal point corruption (CRITICAL — now fixed)
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

### 2. Space in project path breaks Stata batch mode
**Bug:** `run_analysis.py` passed the full do-file path to Stata's `-b do` argument.
Stata truncates paths at spaces, so `Research New/01_data_checks.do` became `Research.do`.
**Fix:** Pass only the filename (not full path); `cwd=project_dir` in subprocess handles the rest.

### 3. Python 3.9 type hint syntax
**Bug:** `str | None` union type hint requires Python 3.10+. System Python is 3.9.6.
**Fix:** Changed to untyped `def find_stata(override=None)`.

### 4. `top_up` is a dummy, not a dollar amount
**Clarification:** `top_up` in the raw data is 0/1, not a dollar amount.
Imputed top-up amount = `max(0, (mw - piece*prod) * tot_hrs)`. Fixed in `01_data_checks.do`.

---

## Do-File Structure

**Legacy monolith:** `dana_minwage_2026_feb.do` (do not edit — reference only)

**Modular files (created 2026-03-04):**

| File | Purpose | Status |
|---|---|---|
| `00_data_prep.do` | Import CSV → `planter_mw_prepped.dta` | ✅ Runs clean |
| `01_data_checks.do` | Coverage, missingness, top-up validity, weather, experience | ✅ Runs clean |
| `02_main_analysis.do` | 2A + 2B DiD, auto event-window, coefplot → `output/` | Not yet run |
| `03_heterogeneity.do` | Ability, phat, contract diff, het DiD | Not yet run |
| `run_analysis.py` | Python master runner | ✅ Working |

**How to run:**
```bash
python3 run_analysis.py              # all steps
python3 run_analysis.py --step 02   # single step
python3 run_analysis.py --from 01   # from step 01 onward
```
Say "run 02" to Claude and it will execute, read the log, and report findings.

**Output locations:**
- Logs: `logs/NN_stepname_YYYYMMDD.log`
- Figures: `output/event_*.png`
- Prepared data: `planter_mw_prepped.dta` (gitignored — rebuilt by 00)

---

## Next Steps

1. Run `02_main_analysis.do` — produces 2A and 2B DiD + event study plots
2. Check event study plots against existing PNGs from the monolith for consistency
3. Run `03_heterogeneity.do` — ability, phat, contract difficulty het DiD
4. Address phat collinearity question (see open questions below)

---

## Open Questions

### phat / ability collinearity
`phat_topup` is constructed from `ability_z` + weather + controls (Section 3.7).
Because `ability_z` enters both the residualization (3.4) and the phat model (3.7),
the two are likely collinear. This makes it hard to separately identify the ability
channel vs. the MW-binding channel in parametric heterogeneity regressions (3.9).
Paradox from earlier runs: high-phat (low ability) workers show smallest response;
low-phat (high ability) workers show largest productivity declines.

---

## Software Requirements

- Stata 19 MP (`/Applications/Stata/StataMP.app/Contents/MacOS/StataMP`)
- Packages: `reghdfe`, `ftools`, `coefplot` (auto-installed by `00_data_prep.do`)
- Python 3.9+ for `run_analysis.py`
- Raw data: `~/Downloads/planter_mw.csv`
