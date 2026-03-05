# Project Overview: Minimum Wage and Productivity in Piece-Rate Work

## Research Question

How does a minimum wage increase affect worker productivity, and what are the mechanisms?
The project focuses on understanding both the average effect and the heterogeneous responses
across workers with different probabilities of the minimum wage binding (`phat`).

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

1. **Piece-rate adjustment**: Firms may lower the piece rate ($ per tree) to offset the wage
   floor, passing costs back to workers — consistent with firm-level incentives
2. **Threat effects**: Workers who risk falling below the floor may face firing threats,
   affecting effort even before the floor binds
3. **Binding probability (`phat`)**: Only ~3% of workers receive top-ups, but output is
   stochastic, so all workers face some probability of the floor binding. Theory predicts
   effects should scale with this probability

## Identification Strategy

**Section 2A — Cross-province DiD (primary):**
BC workers vs. AB or ON workers in the treated year, before and after June 1.
Event study around the June 1 cutoff with auto-selected event window (TAU rule).

**Section 2B — Within-BC placebo years (secondary):**
BC workers in the treated year (2018) vs. BC workers in placebo years (2016, 2017).
Identification is weaker here; 2A is the primary focus.

**Heterogeneity (Section 3) — DDD:**
Triple-difference interacting BC × post with phat quartiles or ability bins to test
whether productivity responses scale with binding probability.

## Baseline Specification

- **Estimator**: `reghdfe`
- **Fixed effects**: Worker FE + Contract FE (TWFE)
- **Clustering**: Worker-level
- **Controls**: Weather (temp, wind, precip, quadratics), experience (quadratic),
  hours, multi-contract day flag, contract lifecycle ramp, season ramp,
  ramp × cumulative seasons interaction
- **Sample floors**: ≥2 hrs/day, ≥50 trees/hr, ≥$30 earnings, positive piece rate and output

## Current Status & Key Findings (Preliminary)

**Status**: Early stage. Inherited from student dissertation; now taking the project forward.

**Headline finding**: The MW increase reduces productivity on average, with heterogeneous
effects by `phat` that are not yet well understood. Paradoxically, the largest productivity
declines appear among highly productive workers with near-zero top-up probability, while
low-productivity workers (high `phat`) respond very little. All workers see a piece-rate
increase post-reform, consistent with firm-level adjustment.

**Open question**: Whether the heterogeneity result reflects a genuine mechanism or a
specification/measurement issue with `phat` construction — in particular, whether `phat`
is too collinear with worker ability to separately identify their effects.

## Next Steps

- [ ] Split `dana_minwage_2026_feb.do` into modular sub-do-files
- [ ] Build `run_analysis.py` master script to run Stata from VS Code
- [ ] Set up log file output for all Stata runs (saved to project directory)
- [ ] Investigate `phat` construction — check R² and variable importance,
      consider weather+MW-only version for cleaner exogenous variation
- [ ] Run event study split by phat quartile to check for pre-existing trends
- [ ] Assess whether ability and phat effects can be separately identified

## Target Output

Journal paper targeting **Journal of Human Resources (JHR)** or
**Journal of Labor Economics (JLE)**.

---

## Project Preferences

### Workflow

- **Primary tool**: Stata, everything runs from do-files (no interactive work)
- **IDE**: VS Code, with Stata called via `run_analysis.py`
- **Version control**: Git. Commits made on explicit user instruction only
- **Log files**: All Stata runs save log files to the project directory;
  logs used for interpretation, sanity checks, and debugging without re-running

### File Organization

- Analysis split into modular sub-do-files, orchestrated by a master Python script
- Date-stamped do-files for major versions (e.g., `dana_minwage_2026_feb.do`)
- Figures exported as PNG to the project directory
- Log files saved to project directory alongside do-files

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
