# Project Overview: Minimum Wage and Productivity in Piece-Rate Work

## Research Question

How does a minimum wage increase affect worker productivity, and what are the mechanisms?
The project focuses on understanding both the average effect and the heterogeneous responses
across workers with different pre-determined exposure to the minimum wage floor (`d̄_ic`).

## Setting

Tree planting in Canada, using payroll data covering **2016–2019**. Workers are paid piece-rate
(per tree planted), with a mandatory top-up if piece-rate earnings fall below the hourly
minimum wage. The key reforms are BC minimum wage increases effective **June 1, 2018** and
**June 1, 2019**, announced well in advance as part of a multi-year schedule published in
early 2018.

The tree planting context is well-suited for this question:
- Rich payroll data with daily worker-level observations on output, hours, and pay
- Clean piece-rate structure makes productivity directly observable
- The June 1 cutoff creates a sharp within-season treatment date
- The setting is illustrative of broader piece-rate and gig economy workers

## Mechanisms of Interest

1. **Piece-rate adjustment**: Firms may raise the piece rate ($ per tree) post-reform to
   attract and retain workers under a higher floor — all workers see higher piece rates,
   which raises the opportunity cost of reduced effort
2. **Income effects / backward-bending labor supply**: Higher piece rates mean high-ability
   workers can hit income targets with less effort. Workers far from the floor benefit most
   from piece-rate increases → largest effort reductions among high-d̄_ic workers
3. **Threat effects**: Workers close to the floor face binding risk and potential firing
   pressure, which may hold effort stable even post-reform — attenuating the response
   among low-d̄_ic workers
4. **Pre-determined exposure (`d̄_ic`)**: The expected log gap between predicted incentive
   earnings and the MW floor, constructed entirely from pre-reform data. Workers differ
   in `d̄_ic` by ability and contract difficulty. Theory and evidence suggest effects
   scale with this measure — but in the opposite direction to simple binding-probability
   predictions (see Current Findings)

## Identification Strategy

**Section 2A — Cross-province DiD (primary):**
BC workers vs. AB workers in the treated year, before and after June 1.
Event study around the June 1 cutoff with auto-selected event window (TAU rule).
AB is the preferred control (ON had its own large MW increase in 2018).

**Section 2B — Within-BC placebo years (secondary):**
BC workers in the treated year (2018) vs. BC workers in placebo years (2016, 2017).
Identification is weaker here; 2A is the primary focus.

**Heterogeneity (Section 3) — DDD:**
Triple-difference interacting BC × post with quartile bins and a continuous parametric
interaction of `d̄_ic` (pre-determined exposure measure). Tests whether productivity
responses scale with pre-reform distance to the MW floor.

## Baseline Specification

- **Estimator**: `reghdfe`
- **Fixed effects**: Worker FE + Contract FE (TWFE)
- **Clustering**: Worker-level
- **Controls**: Weather (temp, wind, precip, quadratics), hours, multi-contract day flag,
  contract lifecycle ramp, season ramp, ramp × cumulative seasons interaction
- **Sample floors**: ≥2 hrs/day, ≥50 trees/hr, ≥$30 earnings, positive piece rate and output

## Current Status & Key Findings (2026-03-06)

**Status**: All four do-files run clean. Sections 2 and 3 complete; results in hand.
Full modular pipeline: `00_data_prep` → `01_data_checks` → `02_main_analysis` → `03_heterogeneity`.

**Average effect (Section 2A):** MW increase reduces `ln_incentive` in BC relative to AB
after June 1, 2018. Consistent across DiD and event study specifications.

**Heterogeneity (Section 3) — Headline finding:**
- Parametric DDD: treat × post = **−0.155** (p=0.004); d̄_ic × treat × post = **−0.143** (p=0.048)
- Quartile bins: Bin 1 (closest to floor) = −0.008 (n.s.); Bin 4 (furthest) = **−0.227** (p=0.014)
- Pattern is **monotonically increasing** — workers furthest from the floor show the largest
  productivity declines, the opposite of what a naïve binding-probability story predicts
- Volatility DDD: workers with higher idiosyncratic earnings volatility show *smaller* declines

**Leading interpretation:** Income effects from piece-rate increases dominate among high-ability
workers (high d̄_ic), who can reduce effort while still earning more than before the reform.
Low-d̄_ic workers may be stabilized by threat effects or genuine binding risk.

**Key validity concern:** Mean reversion in α̂_i — pre-trend check by bin is the most
important pending diagnostic.

## Next Steps

- [ ] Event study split by d̄_ic quartile — check for pre-existing trends by bin
- [ ] Assess mean reversion: regress within-worker FE change on d̄_ic
- [ ] Robustness: weather-only version of d̄_ic (drop α̂_i from exposure measure)
- [ ] 2019 replication: run Section 3 for the second BC reform year
- [ ] Refine figures for paper quality (currently exploratory/monochrome)
- [ ] Draft Section 3 of paper

## Target Output

Journal paper targeting **Journal of Human Resources (JHR)** or
**Journal of Labor Economics (JLE)**.

---

## Project Preferences

### Workflow

- **Primary tool**: Stata, everything runs from do-files (no interactive work)
- **IDE**: VS Code, with Stata called via `run_analysis.py`
- **Version control**: Git. Commits made on explicit user instruction only
- **Log files**: All Stata runs save log files to `logs/`; logs used for interpretation,
  sanity checks, and debugging without re-running

### File Organization

- Analysis split into modular sub-do-files, orchestrated by a master Python script
- Framework documented in `mw_exposure_framework_final.md` (authoritative reference)
- Figures exported as PNG to `output/`
- Log files saved to `logs/` alongside do-files

### Stata Coding Style

- **Heavy commenting** — all sections, decisions, and non-obvious choices explained
- **Globals for all toggles** — outcomes, FE choice, clustering, floors, windows
  all controlled from a single control panel at the top of each do-file
- **Section headers** with clear dividers (`/***...***/` style)
- **Debug output** — `di as text "[DEBUG] ..."` lines to trace macro values
- **Defensive coding** — `cap drop`, `cap macro drop`, `capture confirm` before
  creating variables/macros to allow safe reruns
- **No interactive work** — everything reproducible from the do-file

### Analysis Defaults

- Baseline: Worker FE + Contract FE, clustered by worker
- Event window: Auto-selected via TAU rule (TAU=0.25), falls back to manual
- Outcome: `ln_incentive` as headline; `lprod`, `ln_piece`, `tot_hrs_day` as secondary
- Exploration phase — figure style is minimal, monochrome; will refine for paper

### Figures

- Minimal style, minimal color
- Currently exploratory — not production-ready
- Exported as PNG via `coefplot`
