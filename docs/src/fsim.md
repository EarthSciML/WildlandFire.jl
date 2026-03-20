# FSim Wildfire Risk Components

## Overview

The FSim (Fire Simulation) system is a large-fire risk assessment framework
designed to estimate probabilistic wildfire risk components including burn
probabilities and fire size distributions. The system was developed for
continental-scale application across the United States, operating on 134 Fire
Planning Units (FPUs).

FSim integrates four main modules:
1. **Weather generation**: Synthetic daily weather using autoregressive time series of the Energy Release Component (ERC)
2. **Fire occurrence**: Logistic regression relating ERC to large fire probability
3. **Fire growth**: Spatial fire spread simulation using the minimum travel time algorithm with Rothermel's spread model
4. **Fire suppression**: Statistical containment model based on fire growth rates and fuel type

The components implemented here represent the analytical/statistical sub-models
from FSim that can be expressed as algebraic equation systems. The spatial fire
growth algorithm (MTT) and Monte Carlo simulation loop are procedural and not
included.

**Reference**: Finney, M.A., McHugh, C.W., Grenfell, I.C., Riley, K.L., Short, K.C. (2011).
A simulation of probabilistic wildfire risk components for the continental United States.
*Stochastic Environmental Research and Risk Assessment*, 25:973--1000.
DOI: [10.1007/s00477-011-0462-z](https://doi.org/10.1007/s00477-011-0462-z)

```@docs
FireOccurrenceLogistic
FireContainment
BurnProbability
ERCTimeSeries
FlameLengthCategory
```

## Implementation

### Fire Occurrence Logistic Model

The fire occurrence model uses logistic regression to predict the probability
of at least one large fire starting on a given day, as a function of the
Energy Release Component (ERC) index (Section 2.2, Fig. 4a).

#### State Variables

```@example fsim
using ModelingToolkit, Symbolics, DataFrames, DynamicQuantities
using WildlandFire

sys = FireOccurrenceLogistic()

vars = unknowns(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [string(dimension(ModelingToolkit.get_unit(v))) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

#### Parameters

```@example fsim
params = parameters(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params],
    :Units => [string(dimension(ModelingToolkit.get_unit(p))) for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

#### Equations

```@example fsim
equations(sys)
```

### Fire Containment Model

The fire containment model estimates the probability of successful fire
suppression based on fire growth rate, duration, and fuel type (Section 2.4, Fig. 6).

#### State Variables

```@example fsim
sys_contain = FireContainment()

vars_c = unknowns(sys_contain)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars_c],
    :Units => [string(dimension(ModelingToolkit.get_unit(v))) for v in vars_c],
    :Description => [ModelingToolkit.getdescription(v) for v in vars_c]
)
```

#### Parameters

```@example fsim
params_c = parameters(sys_contain)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params_c],
    :Units => [string(dimension(ModelingToolkit.get_unit(p))) for p in params_c],
    :Description => [ModelingToolkit.getdescription(p) for p in params_c]
)
```

#### Equations

```@example fsim
equations(sys_contain)
```

### Burn Probability Model

Burn probability is the fundamental output of FSim, representing the annual
likelihood of burning at a given location (Sections 2.5, 3).

#### State Variables

```@example fsim
sys_bp = BurnProbability()

vars_bp = unknowns(sys_bp)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars_bp],
    :Units => [string(dimension(ModelingToolkit.get_unit(v))) for v in vars_bp],
    :Description => [ModelingToolkit.getdescription(v) for v in vars_bp]
)
```

#### Parameters

```@example fsim
params_bp = parameters(sys_bp)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params_bp],
    :Units => [string(dimension(ModelingToolkit.get_unit(p))) for p in params_bp],
    :Description => [ModelingToolkit.getdescription(p) for p in params_bp]
)
```

#### Equations

```@example fsim
equations(sys_bp)
```

### ERC Time Series Model

The autoregressive time series model generates synthetic daily ERC(G) values
that capture seasonal trends and day-to-day persistence (Section 2.1, Eqs. 1--4).

#### State Variables

```@example fsim
sys_erc = ERCTimeSeries()

vars_erc = unknowns(sys_erc)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars_erc],
    :Units => [string(dimension(ModelingToolkit.get_unit(v))) for v in vars_erc],
    :Description => [ModelingToolkit.getdescription(v) for v in vars_erc]
)
```

#### Parameters

```@example fsim
params_erc = parameters(sys_erc)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params_erc],
    :Units => [string(dimension(ModelingToolkit.get_unit(p))) for p in params_erc],
    :Description => [ModelingToolkit.getdescription(p) for p in params_erc]
)
```

#### Equations

```@example fsim
equations(sys_erc)
```

### Flame Length Categorization

The flame length categorization model classifies Byram's fireline intensity
into six operational flame length categories (Fig. 10).

#### State Variables

```@example fsim
sys_fl = FlameLengthCategory()

vars_fl = unknowns(sys_fl)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars_fl],
    :Units => [string(dimension(ModelingToolkit.get_unit(v))) for v in vars_fl],
    :Description => [ModelingToolkit.getdescription(v) for v in vars_fl]
)
```

#### Parameters

```@example fsim
params_fl = parameters(sys_fl)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params_fl],
    :Units => [string(dimension(ModelingToolkit.get_unit(p))) for p in params_fl],
    :Description => [ModelingToolkit.getdescription(p) for p in params_fl]
)
```

#### Equations

```@example fsim
equations(sys_fl)
```

## Analysis

### Fire Occurrence Probability vs. ERC (Fig. 4a)

The logistic regression curves show how fire probability increases with ERC,
with larger fire size thresholds yielding lower probabilities. This reproduces
the pattern shown in Fig. 4a of Finney et al. (2011).

```@example fsim
using NonlinearSolve, Plots

sys = FireOccurrenceLogistic()
cs = mtkcompile(sys)

# Simulate curves for different fire size thresholds
# Larger negative β₀ values represent larger minimum fire sizes
erc_range = 0:1:100
thresholds = [
    (β₀=-3.0, β₁=0.05, label="Small fires (4+ ha)"),
    (β₀=-4.5, β₁=0.05, label="Medium fires (120+ ha)"),
    (β₀=-5.5, β₁=0.05, label="Large fires (250+ ha)"),
    (β₀=-7.0, β₁=0.05, label="Very large fires (1200+ ha)"),
]

p = plot(xlabel="ERC (G)", ylabel="Probability of fire occurrence",
    title="Fire Occurrence Probability vs. ERC\n(cf. Finney et al. 2011, Fig. 4a)",
    legend=:topleft)

for th in thresholds
    P_values = Float64[]
    for erc in erc_range
        prob = NonlinearProblem(cs, Dict(
            cs.β₀ => th.β₀, cs.β₁ => th.β₁, cs.ERC => Float64(erc)))
        sol = solve(prob)
        push!(P_values, sol[cs.P_fire])
    end
    plot!(p, collect(erc_range), P_values, label=th.label, linewidth=2)
end
p
```

### Fire Containment Probability (Fig. 6)

The containment probability varies with the number of previous intervals (fire
duration proxy), spread rate, and fuel type. This reproduces the pattern shown
in Fig. 6 of Finney et al. (2011).

```@example fsim
sys_c = FireContainment()
cs_c = mtkcompile(sys_c)

# Model coefficients calibrated to approximate Fig. 6 patterns
α₀, α₁, α₂, α₃ = 0.5, 0.3, -1.5, -0.8

npi_range = 1:20

# Panel (a): No timber fuel types
p1 = plot(xlabel="Number of Days in Interval",
    ylabel="Probability of Containment",
    title="(a) No Timber Fuel Types Present (cf. Fig. 6a)",
    legend=:bottomright, ylim=(0, 1))

for (spread_val, spread_label) in [(0.0, "Low Spread Intervals"), (1.0, "High Spread Intervals")]
    P_values = Float64[]
    for npi in npi_range
        prob = NonlinearProblem(cs_c, Dict(
            cs_c.α₀ => α₀, cs_c.α₁ => α₁, cs_c.α₂ => α₂, cs_c.α₃ => α₃,
            cs_c.NPI => Float64(npi), cs_c.is_high_spread => spread_val,
            cs_c.is_timber => 0.0))
        sol = solve(prob)
        push!(P_values, sol[cs_c.P_contain])
    end
    plot!(p1, collect(npi_range), P_values, label=spread_label, linewidth=2)
end

# Panel (b): Timber fuel types present
p2 = plot(xlabel="Number of Days in Interval",
    ylabel="Probability of Containment",
    title="(b) Timber Fuel Types Present (cf. Fig. 6b)",
    legend=:bottomright, ylim=(0, 1))

for (spread_val, spread_label) in [(0.0, "Low Spread Intervals"), (1.0, "High Spread Intervals")]
    P_values = Float64[]
    for npi in npi_range
        prob = NonlinearProblem(cs_c, Dict(
            cs_c.α₀ => α₀, cs_c.α₁ => α₁, cs_c.α₂ => α₂, cs_c.α₃ => α₃,
            cs_c.NPI => Float64(npi), cs_c.is_high_spread => spread_val,
            cs_c.is_timber => 1.0))
        sol = solve(prob)
        push!(P_values, sol[cs_c.P_contain])
    end
    plot!(p2, collect(npi_range), P_values, label=spread_label, linewidth=2)
end

plot(p1, p2, layout=(2, 1), size=(700, 600))
```

### ERC Time Series Generation (Fig. 3a)

The autoregressive model generates synthetic ERC values that follow the
seasonal trend with realistic daily variability.

```@example fsim
sys_erc = ERCTimeSeries()
cs_erc = mtkcompile(sys_erc)

# Generate a synthetic year of ERC using the AR(1) model
# Seasonal trend: approximate bell-shaped curve peaking in summer
days = 1:365
f_seasonal = [40.0 * sin(pi * (d - 90) / 180)^2 for d in days]

# AR(1) parameters (typical for a western US station)
phi1 = 0.85
rho1 = 0.85
s2 = 150.0

# Innovation variance
prob_var = NonlinearProblem(cs_erc, Dict(
    cs_erc.f_t => 50.0, cs_erc.a_prev => 0.0, cs_erc.a_current => 0.0,
    cs_erc.φ₁ => phi1, cs_erc.ρ₁ => rho1, cs_erc.s² => s2))
sol_var = solve(prob_var)
sigma_innovation = sqrt(sol_var[cs_erc.var_a])

# Simulate 3 years using the AR(1) process
using Random
Random.seed!(42)
n_years = 3
n_days = 365 * n_years
erc_simulated = zeros(n_days)
a_values = zeros(n_days)

for d in 1:n_days
    day_of_year = mod1(d, 365)
    local f_t = f_seasonal[day_of_year]
    a_t = sigma_innovation * randn()

    a_prev_val = d > 1 ? a_values[d-1] : 0.0

    local prob = NonlinearProblem(cs_erc, Dict(
        cs_erc.f_t => f_t,
        cs_erc.a_prev => a_prev_val,
        cs_erc.a_current => a_t,
        cs_erc.φ₁ => phi1,
        cs_erc.ρ₁ => rho1,
        cs_erc.s² => s2))
    local sol = solve(prob)
    erc_simulated[d] = max(0.0, sol[cs_erc.ERC_hat])
    a_values[d] = a_t
end

# Compute seasonal trend for plotting
f_trend = repeat(f_seasonal, n_years)

p_erc = plot(1:n_days, erc_simulated,
    xlabel="Day", ylabel="ERC (G)",
    title="Simulated ERC Time Series (3 years)\n(cf. Finney et al. 2011, Fig. 3a)",
    label="Simulated ERC", alpha=0.7, color=:blue, legend=:topleft)
plot!(p_erc, 1:n_days, f_trend, label="Seasonal trend", linewidth=2, color=:navy)
p_erc
```

### Flame Length Categories (Fig. 10)

The six flame length categories used in FSim correspond to different levels
of fire suppression difficulty and ecological impact.

```@example fsim
sys_fl = FlameLengthCategory()
cs_fl = mtkcompile(sys_fl)

# Sweep flame lengths and show category assignments
fl_range = 0.1:0.1:5.0
active_cats = Int[]

for fl in fl_range
    local prob = NonlinearProblem(cs_fl, Dict(cs_fl.F_L => fl))
    local sol = solve(prob)
    cats = [sol[cs_fl.cat_1], sol[cs_fl.cat_2], sol[cs_fl.cat_3],
             sol[cs_fl.cat_4], sol[cs_fl.cat_5], sol[cs_fl.cat_6]]
    push!(active_cats, findfirst(cats .> 0.5))
end

p_fl = scatter(collect(fl_range), active_cats,
    xlabel="Flame Length (m)",
    ylabel="Category",
    title="Flame Length Categories\n(cf. Finney et al. 2011, Fig. 10)",
    yticks=(1:6, ["0-0.6m", "0.6-1.2m", "1.2-1.8m", "1.8-2.4m", "2.4-3.7m", ">3.7m"]),
    legend=false, markersize=3, color=active_cats)
p_fl
```

### Burn Probability Range (Fig. 7c)

FSim produces burn probabilities spanning approximately four orders of magnitude
across the U.S. (10^-5 to 10^-2). This example demonstrates the relationship
between burn counts and annual burn probability.

```@example fsim
sys_bp = BurnProbability()
cs_bp = mtkcompile(sys_bp)

# Simulate a range of burn counts
n_years_sim = 10000.0
burn_counts = [1, 5, 10, 50, 100, 500, 1000]

p_bp = plot(xlabel="Number of times cell burned",
    ylabel="Annual Burn Probability",
    title="Cell-Level Burn Probability\n(cf. Finney et al. 2011, Sections 2.5, 3)",
    yscale=:log10, legend=:bottomright)

bp_values = Float64[]
for n_burned in burn_counts
    local prob = NonlinearProblem(cs_bp, Dict(
        cs_bp.N_burned => Float64(n_burned),
        cs_bp.N_years => n_years_sim,
        cs_bp.A_burned => 1.0e9,
        cs_bp.A_FPU => 1.0e10,
        cs_bp.P_FL_given_burn => 1.0))
    local sol = solve(prob)
    push!(bp_values, sol[cs_bp.BP])
end
plot!(p_bp, Float64.(burn_counts), bp_values,
    label="BP = N_burned / N_years", linewidth=2, marker=:circle)
p_bp
```
