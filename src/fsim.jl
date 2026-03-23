"""
    FireOccurrenceLogistic(; name=:FireOccurrenceLogistic)

Create a logistic regression model for large fire occurrence probability
as a function of the Energy Release Component (ERC).

The model predicts the probability of at least one large fire starting on a given day,
based on the daily ERC(G) value. This is a key component of the FSim large-fire
simulation system.

# Model Description

The probability of fire occurrence is given by the standard logistic function:

```
P(fire | ERC) = 1 / (1 + exp(-(β₀ + β₁ × ERC)))
```

Where:
- `P(fire | ERC)`: Probability of at least one large fire on a day with given ERC value
- `β₀`: Logistic regression intercept (FPU-specific, fitted from historical data)
- `β₁`: Logistic regression slope coefficient (FPU-specific)
- `ERC`: Energy Release Component index value

The logistic regression coefficients are estimated from historical fire records
for each Fire Planning Unit (FPU). Larger fire size thresholds yield lower
probabilities at all ERC levels.

# Reference

Finney, M.A., McHugh, C.W., Grenfell, I.C., Riley, K.L., Short, K.C. (2011).
A simulation of probabilistic wildfire risk components for the continental United States.
*Stochastic Environmental Research and Risk Assessment*, 25:973–1000.
DOI: 10.1007/s00477-011-0462-z (Section 2.2, Fig. 4a)

# Example

```julia
using ModelingToolkit, WildlandFire, NonlinearSolve

sys = FireOccurrenceLogistic()
compiled_sys = mtkcompile(sys)

# Example: probability of fire at ERC = 60 with typical Western US coefficients
prob = NonlinearProblem(compiled_sys, Dict(
    compiled_sys.β₀ => -5.0,
    compiled_sys.β₁ => 0.05,
    compiled_sys.ERC => 60.0
))
sol = solve(prob)
```
"""
@component function FireOccurrenceLogistic(; name = :FireOccurrenceLogistic)
    @constants begin
        one = 1.0, [description = "Dimensionless unity (dimensionless)", unit = u"1"]
    end

    @parameters begin
        β₀, [description = "Logistic regression intercept (dimensionless)", unit = u"1"]
        β₁, [description = "Logistic regression slope coefficient (dimensionless)", unit = u"1"]
        ERC, [description = "Energy Release Component index (dimensionless)", unit = u"1"]
    end

    @variables begin
        P_fire(t), [description = "Probability of at least one large fire (dimensionless)", unit = u"1"]
    end

    eqs = [
        # Logistic regression — Finney et al. (2011), Section 2.2, Fig. 4a
        P_fire ~ one / (one + exp(-(β₀ + β₁ * ERC))),
    ]

    return System(eqs, t; name)
end


"""
    FireContainment(; name=:FireContainment)

Create a statistical fire containment probability model.

This model estimates the probability of fire containment based on fire growth
rate characteristics and fire duration. It represents the influence of modern
fire management policy on large fire patterns.

# Model Description

The containment probability increases with:
1. Lower fire growth rates (low spread intervals)
2. Longer fire duration (more previous intervals)
3. Non-timber fuel types (easier to suppress)

The model uses a logistic-type function of the number of previous intervals (NPI)
and fire spread rate characteristics:

```
P_contain = 1 / (1 + exp(-(α₀ + α₁ × NPI + α₂ × is_high_spread + α₃ × is_timber)))
```

Where:
- `P_contain`: Probability of fire containment at end of current interval
- `NPI`: Number of previous spread intervals (proxy for fire duration)
- `is_high_spread`: Whether current interval is a high-spread interval (0 or 1)
- `is_timber`: Whether timber fuel types are present (0 or 1)
- `α₀, α₁, α₂, α₃`: Model coefficients

# Reference

Finney, M.A., McHugh, C.W., Grenfell, I.C., Riley, K.L., Short, K.C. (2011).
A simulation of probabilistic wildfire risk components for the continental United States.
*Stochastic Environmental Research and Risk Assessment*, 25:973–1000.
DOI: 10.1007/s00477-011-0462-z (Section 2.4, Fig. 6)

Also based on:
Finney, M.A., Grenfell, I.C., McHugh, C.W. (2009). Modeling large fire containment
using generalized linear mixed model analysis. *For Sci* 55(3):249–255.

# Example

```julia
using ModelingToolkit, WildlandFire, NonlinearSolve

sys = FireContainment()
compiled_sys = mtkcompile(sys)

# Example: containment probability after 5 intervals, low spread, no timber
prob = NonlinearProblem(compiled_sys, Dict(
    compiled_sys.α₀ => 0.5,
    compiled_sys.α₁ => 0.3,
    compiled_sys.α₂ => -1.5,
    compiled_sys.α₃ => -0.8,
    compiled_sys.NPI => 5.0,
    compiled_sys.is_high_spread => 0.0,
    compiled_sys.is_timber => 0.0
))
sol = solve(prob)
```
"""
@component function FireContainment(; name = :FireContainment)
    @constants begin
        one = 1.0, [description = "Dimensionless unity (dimensionless)", unit = u"1"]
    end

    @parameters begin
        α₀, [description = "Containment model intercept (dimensionless)", unit = u"1"]
        α₁, [description = "Containment model NPI coefficient (dimensionless)", unit = u"1"]
        α₂, [description = "Containment model high-spread coefficient (dimensionless)", unit = u"1"]
        α₃, [description = "Containment model timber fuel coefficient (dimensionless)", unit = u"1"]
        NPI, [description = "Number of previous spread intervals (dimensionless)", unit = u"1"]
        is_high_spread, [description = "High spread interval indicator, 0 or 1 (dimensionless)", unit = u"1"]
        is_timber, [description = "Timber fuel type indicator, 0 or 1 (dimensionless)", unit = u"1"]
    end

    @variables begin
        P_contain(t), [description = "Probability of fire containment (dimensionless)", unit = u"1"]
    end

    eqs = [
        # Fire containment probability — Finney et al. (2011), Section 2.4, Fig. 6
        # Based on Finney et al. (2009) generalized linear mixed model
        P_contain ~ one / (one + exp(-(α₀ + α₁ * NPI + α₂ * is_high_spread + α₃ * is_timber))),
    ]

    return System(eqs, t; name)
end


"""
    BurnProbability(; name=:BurnProbability)

Create a burn probability calculation model.

Burn probability is the fundamental output of the FSim fire simulation system.
It represents the likelihood that a given location will burn in any given year,
calculated as the ratio of the number of times a cell burns to the total number
of simulated years.

# Model Description

For a grid cell:
```
BP = N_burned / N_years
```

For an FPU:
```
BP_FPU = A_burned / (A_FPU × N_years)
```

The conditional burn probability at a given flame length category is:
```
BP_conditional(FL) = BP × P(FL | burn)
```

Where:
- `BP`: Annual burn probability (dimensionless, 0-1)
- `N_burned`: Number of times the cell burned in the simulation (dimensionless)
- `N_years`: Total number of simulated years (dimensionless)
- `A_burned`: Total area burned across all simulated years (m²)
- `A_FPU`: Total Fire Planning Unit area (m²)
- `P(FL | burn)`: Conditional probability of flame length category given burning

Historical burn probabilities across the U.S. span approximately 10⁻⁵ to 10⁻²,
with western FPUs generally having higher values.

# Reference

Finney, M.A., McHugh, C.W., Grenfell, I.C., Riley, K.L., Short, K.C. (2011).
A simulation of probabilistic wildfire risk components for the continental United States.
*Stochastic Environmental Research and Risk Assessment*, 25:973–1000.
DOI: 10.1007/s00477-011-0462-z (Sections 2.5, 3)

# Example

```julia
using ModelingToolkit, WildlandFire, NonlinearSolve

sys = BurnProbability()
compiled_sys = mtkcompile(sys)

prob = NonlinearProblem(compiled_sys, Dict(
    compiled_sys.N_burned => 50.0,
    compiled_sys.N_years => 10000.0,
    compiled_sys.A_burned => 1.0e9,
    compiled_sys.A_FPU => 1.0e10,
    compiled_sys.P_FL_given_burn => 0.3
))
sol = solve(prob)
```
"""
@component function BurnProbability(; name = :BurnProbability)
    @parameters begin
        N_burned, [description = "Number of times cell burned in simulation (dimensionless)", unit = u"1"]
        N_years, [description = "Total number of simulated years (dimensionless)", unit = u"1"]
        A_burned, [description = "Total area burned across all simulated years", unit = u"m^2"]
        A_FPU, [description = "Total Fire Planning Unit area", unit = u"m^2"]
        P_FL_given_burn, [description = "Conditional probability of flame length category given burning (dimensionless)", unit = u"1"]
    end

    @variables begin
        BP(t), [description = "Cell-level annual burn probability (dimensionless)", unit = u"1"]
        BP_FPU(t), [description = "FPU-level annual burn probability (dimensionless)", unit = u"1"]
        BP_conditional(t), [description = "Conditional burn probability at given flame length (dimensionless)", unit = u"1"]
    end

    eqs = [
        # Cell-level burn probability — Finney et al. (2011), Section 2.5
        BP ~ N_burned / N_years,

        # FPU-level burn probability — Finney et al. (2011), Section 2.5
        BP_FPU ~ A_burned / (A_FPU * N_years),

        # Conditional burn probability by flame length — Finney et al. (2011), Section 3
        BP_conditional ~ BP * P_FL_given_burn,
    ]

    return System(eqs, t; name)
end


"""
    ERCTimeSeries(; name=:ERCTimeSeries)

Create an autoregressive time series model for Energy Release Component (ERC)
generation.

This model generates synthetic daily ERC(G) values using an autoregressive (AR)
process that captures the temporal autocorrelation structure of historical ERC data.
The model is a key component of the FSim weather generation module.

# Model Description

The ERC time series model consists of three components:

1. **Seasonal trend** `f(t)`: A weighted least squares polynomial fit to the daily
   mean ERC(G) over the historical record.

2. **Autoregressive filter**: Captures day-to-day persistence in ERC driven by
   fuel moisture time-lag effects.

3. **White noise**: Random innovations that add daily variability.

The model generates values according to (Eq. 2):
```
ẑ(t) = f(t) + φ₁·a(t-1) + φ₂·a(t-2) + ... + φ_{t*}·a(t-t*) + a(t)
```

The autoregressive coefficients are obtained from the autocorrelation function (Eq. 1):
```
φ = P_{t*}^{-1} ρ_{t*}
```

The variance of the white noise process is (Eq. 4):
```
var(a(t)) = s²(t) × (1 - ρ₁φ₁ - ρ₂φ₂ - ... - ρ_{t*}φ_{t*})
```

This is a simplified single-lag version suitable for demonstrating the concept
within a continuous ModelingToolkit framework. The full model requires discrete
time steps and matrix operations that are better handled procedurally.

# Reference

Finney, M.A., McHugh, C.W., Grenfell, I.C., Riley, K.L., Short, K.C. (2011).
A simulation of probabilistic wildfire risk components for the continental United States.
*Stochastic Environmental Research and Risk Assessment*, 25:973–1000.
DOI: 10.1007/s00477-011-0462-z (Section 2.1, Eqs. 1–4)

Based on methodology from:
Box, G., Jenkins, G. (1976). Time series analysis: forecasting and control.
2nd edition. Holden Day, San Francisco.

# Example

```julia
using ModelingToolkit, WildlandFire, NonlinearSolve

sys = ERCTimeSeries()
compiled_sys = mtkcompile(sys)

# Example: compute variance correction factor
prob = NonlinearProblem(compiled_sys, Dict(
    compiled_sys.f_t => 50.0,
    compiled_sys.a_prev => 3.0,
    compiled_sys.a_current => 1.5,
    compiled_sys.φ₁ => 0.7,
    compiled_sys.ρ₁ => 0.75,
    compiled_sys.s² => 100.0
))
sol = solve(prob)
```
"""
@component function ERCTimeSeries(; name = :ERCTimeSeries)
    @constants begin
        one = 1.0, [description = "Dimensionless unity (dimensionless)", unit = u"1"]
    end

    @parameters begin
        f_t, [description = "Seasonal trend value at current time step (dimensionless)", unit = u"1"]
        a_prev, [description = "White noise value at previous time step (dimensionless)", unit = u"1"]
        a_current, [description = "White noise value at current time step (dimensionless)", unit = u"1"]
        φ₁, [description = "First-order autoregressive coefficient (dimensionless)", unit = u"1"]
        ρ₁, [description = "First-order autocorrelation (dimensionless)", unit = u"1"]
        s², [description = "Observed daily variance of ERC (dimensionless)", unit = u"1"]
    end

    @variables begin
        ERC_hat(t), [description = "Simulated daily ERC value (dimensionless)", unit = u"1"]
        var_a(t), [description = "Variance of white noise process (dimensionless)", unit = u"1"]
        σ²_daily(t), [description = "Daily variance of ERC around seasonal trend (dimensionless)", unit = u"1"]
    end

    eqs = [
        # Simulated ERC value — Eq. 2, Finney et al. (2011)
        # Simplified to first-order autoregressive: ẑ(t) = f(t) + φ₁·a(t-1) + a(t)
        ERC_hat ~ f_t + φ₁ * a_prev + a_current,

        # Variance of white noise — Eq. 4, Finney et al. (2011)
        var_a ~ s² * (one - ρ₁ * φ₁),

        # Daily variance — Eq. 3, Finney et al. (2011)
        σ²_daily ~ var_a / (one - ρ₁ * φ₁),
    ]

    return System(eqs, t; name)
end


"""
    FlameLengthCategory(; name=:FlameLengthCategory)

Create a flame length categorization model.

This model classifies fireline intensity into six flame length categories used
in the FSim fire risk assessment system. Flame length is an empirical
transformation of Byram's fireline intensity and is more interpretable than
raw intensity units for fire management purposes.

The six categories correspond to different fire suppression difficulty levels
and ecological impact thresholds.

# Flame Length Categories

| Category | Flame Length | Description |
|----------|-------------|-------------|
| 1 | ≤ 0.6 m (≤ 2 ft) | Low intensity, easily suppressed |
| 2 | 0.6–1.2 m (2–4 ft) | Moderate, hand line effective |
| 3 | 1.2–1.8 m (4–6 ft) | High, mechanical equipment needed |
| 4 | 1.8–2.4 m (6–8 ft) | Very high, difficult suppression |
| 5 | 2.4–3.7 m (8–12 ft) | Extreme, crown fire potential |
| 6 | > 3.7 m (> 12 ft) | Extreme, full crown fire |

# Reference

Finney, M.A., McHugh, C.W., Grenfell, I.C., Riley, K.L., Short, K.C. (2011).
A simulation of probabilistic wildfire risk components for the continental United States.
*Stochastic Environmental Research and Risk Assessment*, 25:973–1000.
DOI: 10.1007/s00477-011-0462-z (Fig. 10)

Byram, G.M. (1959). Combustion of forest fuels. In: Davis KP (ed) Forest fire:
control and use. McGraw-Hill, New York.

# Example

```julia
using ModelingToolkit, WildlandFire, NonlinearSolve

sys = FlameLengthCategory()
compiled_sys = mtkcompile(sys)

# Classify a flame length of 2.0 m
prob = NonlinearProblem(compiled_sys, Dict(
    compiled_sys.F_L => 2.0
))
sol = solve(prob)
```
"""
@component function FlameLengthCategory(; name = :FlameLengthCategory)
    @constants begin
        one = 1.0, [description = "Dimensionless unity (dimensionless)", unit = u"1"]
        zero = 0.0, [description = "Dimensionless zero (dimensionless)", unit = u"1"]
        # Flame length thresholds — Finney et al. (2011), Fig. 10
        FL_1 = 0.6, [description = "FL category 1 upper bound", unit = u"m"]
        FL_2 = 1.2, [description = "FL category 2 upper bound", unit = u"m"]
        FL_3 = 1.8, [description = "FL category 3 upper bound", unit = u"m"]
        FL_4 = 2.4, [description = "FL category 4 upper bound", unit = u"m"]
        FL_5 = 3.7, [description = "FL category 5 upper bound", unit = u"m"]
    end

    @parameters begin
        F_L, [description = "Flame length (Byram)", unit = u"m"]
    end

    @variables begin
        cat_1(t), [description = "In FL category 1: ≤ 0.6 m (dimensionless)", unit = u"1"]
        cat_2(t), [description = "In FL category 2: 0.6–1.2 m (dimensionless)", unit = u"1"]
        cat_3(t), [description = "In FL category 3: 1.2–1.8 m (dimensionless)", unit = u"1"]
        cat_4(t), [description = "In FL category 4: 1.8–2.4 m (dimensionless)", unit = u"1"]
        cat_5(t), [description = "In FL category 5: 2.4–3.7 m (dimensionless)", unit = u"1"]
        cat_6(t), [description = "In FL category 6: > 3.7 m (dimensionless)", unit = u"1"]
    end

    eqs = [
        # Flame length categorization — Finney et al. (2011), Fig. 10
        # Each category is 1 if the flame length falls in that range, 0 otherwise
        cat_1 ~ ifelse(F_L <= FL_1, one, zero),
        cat_2 ~ ifelse(F_L > FL_1, one, zero) * ifelse(F_L <= FL_2, one, zero),
        cat_3 ~ ifelse(F_L > FL_2, one, zero) * ifelse(F_L <= FL_3, one, zero),
        cat_4 ~ ifelse(F_L > FL_3, one, zero) * ifelse(F_L <= FL_4, one, zero),
        cat_5 ~ ifelse(F_L > FL_4, one, zero) * ifelse(F_L <= FL_5, one, zero),
        cat_6 ~ ifelse(F_L > FL_5, one, zero),
    ]

    return System(eqs, t; name)
end
