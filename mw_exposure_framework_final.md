# Pre‑Determined Exposure to the Minimum‑Wage Floor in Piece‑Rate Work

*(Updated specification incorporating all design decisions and
implementation clarifications)*

------------------------------------------------------------------------

# Purpose

This document defines a **robust empirical framework** for constructing
a *pre‑determined exposure measure* to a minimum‑wage floor in a
stochastic piece‑rate environment.

The framework must satisfy:

1.  **No post‑treatment contamination**
2.  **Predictability using only pre‑policy information**
3.  **Compatibility with fixed‑effects earnings models**
4.  **Applicability to both BC (treated) and AB (control)** despite
    differences in contract support
5.  **Explicit handling of missing worker or contract FE support**
6.  **Clear separation between persistent exposure and realized daily
    exposure**
7.  **Transparent treatment of residual distributions and variance**
8.  **Interpretability of heterogeneity coefficients via
    standardization**

The exposure variable measures **risk of falling below the minimum wage
under the pre‑policy regime**, not realized binding after the reform.

------------------------------------------------------------------------

# 1 Economic Environment

Hourly piece‑rate earnings:

    r_ict = p_ct * q_ict

where

-   `p_ct` = piece rate\
-   `q_ict` = productivity (trees/hour)

Observed hourly compensation:

    w_ict = max(r_ict , M_t)

The minimum‑wage floor binds when

    r_ict < M_t

The relevant economic quantity is the **probability that piece‑rate
earnings fall below the minimum wage**.

------------------------------------------------------------------------

# 2 Key Design Decision: Use the Pre‑Policy Minimum Wage

Define

    M_pre = minimum wage before the reform

All exposure calculations use `M_pre`, including post‑period
observations.

Distance to the floor therefore becomes

    d = log(r̂) − log(M_pre)

Using the pre‑policy minimum wage ensures that the exposure variable:

-   does **not mechanically change after the reform**
-   remains a **pre‑determined characteristic**
-   reflects **risk under the original compensation regime**

------------------------------------------------------------------------

# 3 Pre‑Period Earnings Model

Define log hourly piece‑rate earnings

    lr_ict = log(piece_ct * prod_ict)

Estimate

    lr_ict = α_i + γ_c + W_ict β + Z_ict δ + u_ict

where

-   `α_i` = worker earnings component
-   `γ_c` = contract earnings component
-   `W_ict` = weather variables
-   `Z_ict` = optional lifecycle controls

------------------------------------------------------------------------

# 4 Estimation Samples by Province

### BC (treated province)

Estimate earnings model using **BC pre‑policy observations only**.

Reason: avoids contamination from the minimum‑wage reform.

### AB (control province)

Estimate earnings model using **AB pre + post observations**.

Reason: AB has only one contract in the pre‑period, so contract
heterogeneity cannot otherwise be identified. Because AB is untreated,
post observations do not introduce policy contamination.

------------------------------------------------------------------------

# 5 Recovering Fixed Effects for Prediction

Prediction requires recovering absorbed fixed effects.

Predicted log earnings:

    l r̂_ict = α̂_i + γ̂_c + W_ict β̂ + Z_ict δ̂

Implementation must therefore recover

    α̂_i
    γ̂_c

from the estimation procedure (e.g. FE recovery or savefe).

------------------------------------------------------------------------

# 6 Support Requirements for Prediction

Prediction requires **pre‑period support** for both worker and contract
components.

Three cases arise.

### Case 1: Contract appears in pre and post

    γ̂_c identified

Prediction valid.

### Case 2: Contract appears only in pre

No issue. The contract simply disappears after the reform.

### Case 3: Contract appears only in post

    γ̂_c not identified

### Default strategy

The exclusion rule applies **within each province's estimation
sample**.

-   For BC, the estimation sample is BC pre‑period only. Contracts
    that appear only after June 1 in BC have no support in this
    sample and their observations are **excluded from exposure
    construction**.

-   For AB, the estimation sample is AB pre + post. All AB 2018
    contracts appear in this sample and receive identified γ̂_c.
    **No AB observations are excluded under this rule.**

In the 2018 data, AB has only one contract active before June 1
(contract B‑18‑0408, 331 observations). The remaining five AB
contracts are post‑only and are identified from the untreated
post‑period. This confirms the need for the province‑specific
estimation samples defined in Section 4.

This preserves interpretation of exposure as a **pre‑policy
characteristic** while retaining the full AB sample.

### Optional robustness

Contract effects may be imputed from contract observables, but this
introduces model dependence and is recommended only as a robustness
exercise.

------------------------------------------------------------------------

# 7 Structural vs Realized Distance to the Floor

Two distinct exposure measures are constructed.

## 7.1 Pre‑period mean distance (heterogeneity sorting variable)

The heterogeneity sorting variable is the pre‑period mean of realized
daily distance, averaged across all pre‑policy observations for
worker i on contract c:

    d̄_ic = α̂_i + γ̂_c + W̄_ic,pre β̂ − log(M_pre)

where W̄_ic,pre β̂ is the average weather component over the worker's
pre‑period observations on contract c.

Because weather is measured at the contract level (one station per
camp), W̄_ic,pre β̂ = W̄_c,pre β̂ and is constant across workers on
the same contract.

Properties:

-   time‑invariant (pre‑period average)
-   purely pre‑determined (uses only pre‑period observations)
-   incorporates average weather conditions experienced pre‑reform

Rationale: expected earnings — and therefore binding risk — depend on
the worker's full earnings environment, including the weather
conditions typical of their contract. Excluding weather would
understate exposure for workers on contracts in systematically adverse
conditions.

This variable is used for **heterogeneity sorting**.

------------------------------------------------------------------------

## 7.2 Realized daily distance

Replaces the pre‑period average weather with the actual weather on
day t:

    d_ict = α̂_i + γ̂_c + W_ict β̂ − log(M_pre)

Properties:

-   varies day‑to‑day
-   reflects realized environmental conditions rather than pre‑period
    averages

The difference between d_ict and d̄_ic is purely the weather
deviation from its pre‑period mean: W_ict β̂ − W̄_ic,pre β̂.

Used for **probability calculations**.

------------------------------------------------------------------------

# 8 Residual Variance

Residuals from the earnings model:

    û_ict = lr_ict − l r̂_ict

Construct residual volatility

    σ_resid

Possible aggregation levels:

-   worker
-   contract
-   worker‑contract

Worker‑contract is preferred if data allow sufficient observations.

------------------------------------------------------------------------

# 9 Weather‑Driven Variance

Weather component of earnings:

    l r̂_weather = W_ict β̂

Weather deviation relative to average conditions:

    Δ_weather = W_ict β̂ − W̄ β̂

Weather volatility:

    σ_weather

### Identification note

Weather is measured at the **contract level**.

Thus

    σ_weather_c

is constant within contracts.

In DiD specifications with contract FE:

-   the level of `σ_weather_c` is absorbed
-   interactions with treatment and post remain identified

However identification relies entirely on **cross‑contract variation**,
which may be weak if the number of contracts is small.

### Recommended alternative

Construct worker‑level weather volatility

    σ_weather_i

defined as the worker's average exposure to contract weather volatility
across pre‑period contracts.

This provides a **more stable heterogeneity dimension**.

------------------------------------------------------------------------

# 10 Binding Probability

Define

    π_ict = Pr(lr_ict < log(M_pre))

Approximate using empirical residual distributions.

### Default: pooled empirical CDF

    F̂(z) = (1/N_pre) Σ 1{û_n ≤ z}

Binding probability

    π̂_ict = F̂(log(M_pre) − l r̂_ict)

Advantages:

-   stable estimate
-   avoids small‑sample noise

### Robustness: worker‑specific CDF

    F̂_i(z) = (1/N_i) Σ_t 1{û_it ≤ z}

Allows worker‑specific residual volatility.

### Not recommended as baseline

Worker‑contract CDFs are often unstable because many worker‑contract
pairs have few observations.

------------------------------------------------------------------------

# 11 Exposure Objects

The framework produces the following core variables.

### Expected log earnings

    l r̂_ict

### Pre‑period mean distance (het‑sort variable)

    d̄_ic = α̂_i + γ̂_c + W̄_ic,pre β̂ − log(M_pre)

### Realized daily distance

    d_ict = α̂_i + γ̂_c + W_ict β̂ − log(M_pre)

### Volatility measures

    σ_resid
    σ_weather

### Binding probability (optional)

    π̂_ict

------------------------------------------------------------------------

# 12 Standardization of Structural Exposure

The pre‑period mean distance variable is measured in **log earnings
units**, which makes raw coefficients difficult to interpret.

Define standardized exposure

    d̃_ic = (d̄_ic − mean(d̄_ic)) / sd(d̄_ic)

Mean and standard deviation are computed **within the estimation sample
used in the DiD regression**.

Continuous heterogeneity regressions should use

    d̃_ic × treated × post

Interpretation:

The coefficient represents the difference in treatment effects between
workers **one standard deviation apart in pre‑reform exposure to the
minimum‑wage floor** (reflecting both persistent earning capacity and
average pre‑period weather conditions).

------------------------------------------------------------------------

# 13 Worker vs Contract Heterogeneity

Recovered fixed effects allow two dimensions of heterogeneity:

Worker component

    α̂_i

Contract component

    γ̂_c

These represent **reduced‑form earnings environments**, not structural
productivity parameters.

The main exposure object remains

    α̂_i + γ̂_c

------------------------------------------------------------------------

# 14 Implementation Workflow

1.  Estimate BC earnings model using BC pre‑policy data.
2.  Estimate AB earnings model using AB pre + post data.
3.  Recover worker and contract fixed effects.
4.  Restrict exposure construction to observations with support in
    that province's estimation sample (BC pre‑period contracts for
    BC; all 2018 AB contracts for AB).
5.  Construct predicted earnings `l r̂_ict`.
6.  Compute pre‑period mean distance `d̄_ic` (including average
    pre‑period weather component).
7.  Standardize exposure to obtain `d̃_ic`.
8.  Construct residual and weather volatility measures.
9.  Use `d̃_ic` as the main heterogeneity sorting variable.
10. Optionally construct binding probabilities using empirical CDFs.
11. Use worker‑level weather volatility for main DiD heterogeneity
    specifications.

------------------------------------------------------------------------

# Bottom Line

The exposure variable measures **pre‑policy risk of falling below the
minimum wage floor**, not realized binding after the reform.

Achieving this requires:

-   pre‑period earnings models estimated separately by province
-   recovery of worker and contract fixed effects
-   use of the pre‑policy minimum wage
-   province‑specific handling of contract support (BC pre‑only;
    AB pre + post)
-   separation between pre‑period mean distance (het sorting) and
    realized daily distance (probability calculations)
-   inclusion of average pre‑period weather in the het‑sort variable
-   careful specification of residual distributions
-   standardization for interpretable heterogeneity coefficients.

The resulting measure provides a **clean, pre‑determined indicator of
workers' incentive exposure to the minimum‑wage floor**.
