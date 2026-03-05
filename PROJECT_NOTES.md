# Planter Minimum Wage Project — Notes

## Research Question

Does a minimum wage increase affect tree planter worker outcomes? Specifically, the analysis exploits
BC's June 1 minimum wage cutoff to study effects on incentive pay (piece-rate × productivity),
productivity (trees/hour), piece rates, and hours worked.

---

## Data

**File:** `planter_mw.csv`

**Unit of observation:** Worker × contract × day

**Key variables:**

| Variable | Description |
|---|---|
| `name` | Worker identifier (string) |
| `contract` | Contract identifier (string) |
| `date`, `year`, `month` | Date fields |
| `province` | BC, AB, or ON |
| `piece` | Piece rate ($ per tree) |
| `prod` | Productivity (trees per hour) |
| `tot_hrs` | Total hours worked |
| `top_up` | Top-up pay (MW supplement) |
| `experience` | Cross-season experience (years) |
| `temp`, `total_precip`, `wind_spd` | Weather controls |
| `earnings`, `hrly_earn` | Total and hourly earnings |
| `mw` | Minimum wage applicable |

**Constructed outcomes:**

| Variable | Formula |
|---|---|
| `ln_incentive` | ln(piece × prod) |
| `lprod` | ln(prod) — trees/hour |
| `ln_piece` | ln(piece) |

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

## Output Files

Event study plots saved as PNG in the working directory:

| File | Contents |
|---|---|
| `event_2A_ln_incentive_2018_AB.png` | Cross-province event study, ln_incentive, 2018, AB control |
| `event_2A_ln_incentive_2019_AB.png` | Same for 2019 |
| `event_2A_ln_piece_2018_AB.png` | Cross-province event study, ln_piece |
| `event_2A_ln_prod_2018_AB.png` | Cross-province event study, ln_prod |
| `event_2A_lprod_2018_AB.png` | Cross-province event study, lprod |
| `event_2A_lprod_2019_AB.png` | Same for 2019 |
| `event_2B_ln_incentive_2018.png` | Within-BC placebo event study, ln_incentive |
| `event_2B_ln_incentive_2019.png` | Same for 2019 |
| `event_2B_ln_prod_2018.png` | Within-BC placebo event study |
| `event_2B_lprod_2018.png` | Within-BC placebo event study |
| `es_bc_placebo_by_ability_ln_incentive_2019.png` | Placebo by ability group |

---

## Data Issues Found (2026-03-04)

### Weather Sentinel Values — CRITICAL
`temp` and `wind_spd` contain large numeric sentinel codes (~1e+14) for missing weather
readings. These were not caught by `destring` because they are already numeric.
After replacing sentinels with missing in `00_data_prep.do`:
- `temp` missing: **88.7%** (103,497 obs)
- `wind_spd` missing: **87.5%** (102,177 obs)
- `total_precip` missing: 9.3% (acceptable)

**Implication:** Including `temp` or `wind_spd` as controls drops ~88% of the sample.
**Decision needed:** Drop weather controls from baseline, or investigate coverage in raw CSV.
Thresholds used: temp outside [-200, 500], precip > 5000, wind > 500 → set to missing.

### top_up is a 0/1 Dummy
`top_up` records whether a worker received a top-up, not the dollar amount.
Constructed top-up amount: `max(0, (mw - piece*prod) * tot_hrs)`
Validity check: only 37 misclassified obs out of 116,742 (essentially error-free).

### Ontario as Control Province
ON had its own large MW increase in 2018 ($11.40 → $14.00), causing ON binding rate
to spike to 23.9% in 2018 vs ~15% in other years. **ON is not a clean control for 2018.**
AB binding rate is flat (3.5–4.3%) — AB is the preferred control.

### MW Binding Rates by Province × Year
| Province | 2016 | 2017 | 2018 | 2019 |
|---|---|---|---|---|
| BC | 2.1% | 3.1% | 3.4% | 5.9% |
| AB | 4.0% | 3.7% | 4.3% | 3.4% |
| ON | 13.2% | 17.6% | 23.9% | 13.4% |

BC 2019 jump to 5.9% reflects the second BC MW increase.

---

## Do-File Structure

**Legacy monolith:** `dana_minwage_2026_feb.do` (do not edit)

**Modular files (created 2026-03-04):**
- `00_data_prep.do` → `01_data_checks.do` → `02_main_analysis.do` → `03_heterogeneity.do`
- Runner: `python3 run_analysis.py [--step NN] [--from NN]`
- Logs: `logs/NN_name_YYYYMMDD.log`; Figures: `output/`

**Legacy monolith — original `dana_minwage_2026_feb.do`:**

| Section | Description |
|---|---|
| (0) Packages | Install `reghdfe`, `ftools`, `coefplot` if missing |
| (1A) Data prep | Import CSV, type conversions, ID encoding, outcomes |
| (1B) String copies | Safe `str` copies to avoid `strL` issues in collapse/reshape |
| (1C) Event time | Year-specific June 1 event time (`post_y`, `event_week_y`) |
| (1D) Lifecycle | Contract lifecycle, season ramp, cross-season experience |
| (2.1) Control panel | All globals for outcomes, controls, FE, floors, windows |
| (2.2) Build globals | Construct `$X`, `$LIFECYCLE`, `$ABS`, `$SAMPLE_OK` from toggles |
| (2A) Cross-province | DiD + event study, BC vs AB/ON, treated year |
| (2B) Placebo years | DiD + event study, within BC, treated vs placebo years |

---

## Software Requirements

- Stata with: `reghdfe`, `ftools`, `coefplot` (auto-installed if missing)
- Data path hardcoded: `/Users/danaandersen/Downloads/planter_mw.csv`
