
# Construction of a Pre‑Determined Binding‑Risk Measure in Piece‑Rate Settings

## Purpose

This document formalizes an approach for constructing a measure of exposure to a minimum‑wage floor in a piece‑rate compensation environment with stochastic productivity.

The goal is to define an empirical object that is:

1. Economically meaningful and directly related to the worker's incentive environment.
2. Implementable using pre‑policy information.
3. Not mechanically contaminated by post‑policy outcomes.
4. Modular so that both simple and structural variants can be explored.

This document is written primarily to guide downstream implementation by an AI agent.

---

# 1. Economic Environment

## 1.1 Piece‑Rate Earnings and the Minimum Wage Floor

Let worker `i` on contract `c` on day `t` produce hourly output:

```
q_ict
```

Let the piece rate be:

```
p_ct
```

Hourly piece‑rate earnings are:

```
r_ict = p_ct * q_ict
```

In the data this corresponds to:

```
r_ict = piece_ct * prod_ict
```

Let

```
M_t
```

denote the applicable minimum wage.

Observed hourly earnings are therefore

```
w_ict = max(r_ict, M_t)
```

The minimum wage binds whenever

```
r_ict < M_t
```

---

# 2. Target Object

## 2.1 Expected Piece‑Rate Earnings

The first key object is expected hourly piece‑rate earnings conditional on predetermined information:

```
r̂_ict = E[r_ict | I_ict]
```

It is preferable to work in logs:

```
lr_ict = ln(r_ict)
```

and

```
l r̂_ict = E[lr_ict | I_ict]
```

---

## 2.2 Distance to the Floor

The simplest exposure measure is the distance between expected earnings and the minimum wage floor.

In levels:

```
d_L = r̂_ict − M_t
```

In logs:

```
d = l r̂_ict − ln(M_t)
```

Interpretation:

* `d > 0` expected earnings exceed the floor
* `d < 0` expected earnings fall below the floor
* `d ≈ 0` worker is close to the threshold

---

## 2.3 Binding Probability

Define

```
π_ict = Pr(r_ict < M_t | I_ict)
```

In logs:

```
π_ict = Pr(lr_ict < ln(M_t) | I_ict)
```

Let

```
lr_ict = l r̂_ict + ε_ict
```

Then

```
π_ict = Pr(ε_ict < ln(M_t) − l r̂_ict)
```

Thus binding probability depends on

1. expected earnings relative to the floor
2. dispersion of earnings shocks

---

# 3. Empirical Decomposition

The proposal is to estimate the following structure in the pre‑period.

```
lr_ict = α_i + γ_c + W_ict β + Z_ict δ + u_ict
```

Where

* `α_i` = worker fixed effect
* `γ_c` = contract fixed effect
* `W_ict` = weather variables
* `Z_ict` = optional lifecycle controls
* `u_ict` = residual

Interpretation:

* worker ability component
* contract difficulty component
* systematic weather component
* unexplained shock component

This regression must be estimated using **pre‑period observations only**.

---

# 4. Predicted Earnings

Using pre‑period estimates construct predicted earnings for all observations:

```
l r̂_ict = α̂_i + γ̂_c + W_ict β̂ + Z_ict δ̂
```

This provides predicted piece‑rate earnings under the **pre‑policy environment**.

---

# 5. Variance Decomposition

Two distinct sources of volatility should be constructed.

## 5.1 Residual Variance

Residuals:

```
û_ict = lr_ict − l r̂_ict
```

Residual volatility can be measured at a chosen aggregation level:

* worker
* contract
* worker‑contract

Denote

```
σ_resid_g
```

---

## 5.2 Weather‑Driven Variance

Weather component of earnings:

```
l r̂_weather = W_ict β̂
```

Define average weather earnings

```
l r̂_weather_avg = W̄_g β̂
```

Weather deviation

```
Δ_weather = l r̂_weather − l r̂_weather_avg
```

Weather volatility

```
σ_weather_g
```

---

# 6. Core Exposure Objects

Four objects should initially be constructed.

### Expected Earnings

```
l r̂_ict
```

### Distance to Floor

```
d_ict = l r̂_ict − ln(M_t)
```

### Residual Volatility

```
σ_resid_g
```

### Weather Volatility

```
σ_weather_g
```

These should first be analyzed **separately**.

---

# 7. Optional Probability Transformation

Once distance and variance objects behave sensibly a binding probability can be constructed.

General form:

```
π̂_ict = F((ln(M_t) − l r̂_ict) / σ_g)
```

Possible choices for `σ_g`:

1. residual variance only
2. weather variance only
3. combined variance

Combined variance approximation:

```
σ² = σ_resid² + σ_weather²
```

---

# 8. Empirical Implementation Stages

## Stage 1 — Estimate Pre‑Period Earnings Model

Estimate

```
lr_ict = α_i + γ_c + W_ict β + Z_ict δ + u_ict
```

using **pre‑policy data only**.

Outputs:

* worker effects
* contract effects
* weather coefficients
* residuals

---

## Stage 2 — Construct Expected Earnings

```
l r̂_ict = α̂_i + γ̂_c + W_ict β̂ + Z_ict δ̂
```

---

## Stage 3 — Construct Distance to Floor

```
d_ict = l r̂_ict − ln(M_t)
```

---

## Stage 4 — Construct Variance Measures

Compute

```
σ_resid_g
σ_weather_g
```

---

## Stage 5 — Heterogeneity Analysis

Estimate treatment heterogeneity using

* distance to floor
* residual volatility
* weather volatility

---

## Stage 6 — Construct Binding Probability (Optional)

```
π̂_ict = F(ln(M_t) − l r̂_ict)
```

using either

* empirical residual CDF
* standardized variance version

---

# 9. Key Interpretation Rules

1. Expected earnings are **pre‑determined predictions** under the pre‑policy earnings process.
2. Distance to the floor is the most transparent exposure measure.
3. Variance objects capture distinct forms of earnings risk.
4. Binding probability should only be constructed once the underlying components behave sensibly.

---

# 10. Bottom Line

This framework replaces a reduced‑form top‑up prediction with a structurally motivated exposure measure based on

* expected piece‑rate earnings
* distance to the minimum wage floor
* earnings volatility

This produces a cleaner mapping between the empirical design and the theoretical incentive environment.
