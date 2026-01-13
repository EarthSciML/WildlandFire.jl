# Usage Guide

This guide covers common usage patterns and advanced features of WildlandFire.jl.

## Creating Models

### Using the Factory Function

The recommended way to create a Rothermel model is using the `Rothermel()` factory function:

```@example usage1
using WildlandFire

# Create a model with default settings (this is what we'll use later)
sys = Rothermel()

# You can also create a model with custom name
sys_custom = Rothermel(name=:my_fire)

# Or create unsimplified model for further manipulation
sys_unsimplified = Rothermel(simplify=false)
```

### Pre-built Model

For backward compatibility, a pre-built simplified model is available:

```julia
using WildlandFire

# Use the pre-built model
sys_prebuilt = rothermel_simplified
```

## Defining Parameters

Parameters must be provided as a vector of pairs:

```@example usage1
params = [
    h => 8000.0,          # Heat content (Btu/lb)
    S_T => 0.0555,        # Total mineral content
    S_e => 0.010,         # Effective mineral content
    ρ_p => 32.0,          # Particle density (lb/ft³)
    σ => 3500.0,          # Surface-area-to-volume ratio (ft²/ft³)
    w_o => 0.138,         # Fuel load (lb/ft²)
    δ => 1.0,             # Fuel bed depth (ft)
    M_x => 0.12,          # Moisture of extinction
    M_f => 0.05,          # Fuel moisture
    U => 352.0,           # Wind speed (ft/min)
    tan_ϕ => 0.0,         # Slope
]
```

### Required Parameters

All 11 parameters must be specified:

| Parameter | Symbol | Description | Typical Values |
|-----------|--------|-------------|----------------|
| Heat content | `h` | Low heat content of fuel | 7000-9000 Btu/lb |
| Total mineral | `S_T` | Total mineral content | 0.02-0.10 fraction |
| Effective mineral | `S_e` | Effective mineral content | 0.005-0.02 fraction |
| Particle density | `ρ_p` | Oven-dry density | 28-40 lb/ft³ |
| SAV ratio | `σ` | Surface-area-to-volume | 100-4000 ft²/ft³ |
| Fuel load | `w_o` | Oven-dry fuel load | 0.01-1.0 lb/ft² |
| Bed depth | `δ` | Fuel bed depth | 0.1-10 ft |
| Moist. extinct. | `M_x` | Moisture of extinction | 0.12-0.40 fraction |
| Fuel moisture | `M_f` | Current moisture | 0.01-0.30 fraction |
| Wind speed | `U` | Wind at midflame height | 0-2000 ft/min |
| Slope | `tan_ϕ` | Slope steepness (rise/run) | 0-1.0 fraction |

## Solving the Model

### Basic Solution

```@example usage1
using OrdinaryDiffEqDefault

# Create the problem
prob = ODEProblem(sys, [], (0.0, 1.0), params)

# Solve
sol = solve(prob)

# Extract rate of spread
R_value = sol[R][end]  # ft/min
```

### Accessing All Variables

You can access any intermediate variable from the solution:

```@example usage1
# Main output
rate_of_spread = sol[R][end]                    # ft/min

# Reaction characteristics
reaction_intensity = sol[I_R][end]              # Btu/ft²·min
reaction_velocity = sol[Γ][end]                 # min⁻¹

# Fuel bed properties
packing_ratio = sol[β][end]                     # dimensionless
optimum_packing = sol[β_opt][end]               # dimensionless
bulk_density = sol[ρ_b][end]                    # lb/ft³

# Damping coefficients
moisture_damping = sol[η_M][end]                # 0-1
mineral_damping = sol[η_s][end]                 # 0-1

# Wind and slope effects
wind_factor = sol[ϕ_w][end]                     # dimensionless
slope_factor = sol[ϕ_s][end]                    # dimensionless

# Heat transfer
propagating_flux_ratio = sol[ξ][end]            # dimensionless
heat_of_preignition = sol[Q_ig][end]            # Btu/lb
effective_heating = sol[ε][end]                 # dimensionless

# No-wind, no-slope rate
base_rate_of_spread = sol[R_0][end]             # ft/min
```

## Unit Conversions

### Rate of Spread

```@example usage1
# Get rate in ft/min
R_ftmin = sol[R][end]

# Convert to other units
R_mmin = R_ftmin * 0.3048           # meters per minute
R_mhr = R_ftmin * 18.288            # meters per hour
R_kmhr = R_ftmin * 0.018288         # kilometers per hour
R_chains_hr = R_ftmin * 1.829       # chains per hour (wildland fire standard)
R_mph = R_ftmin / 88.0              # miles per hour
```

### Wind Speed

```@example usage7
# Convert from mph to ft/min
wind_mph = 10.0
U = wind_mph * 88.0                 # ft/min

# Convert from m/s to ft/min
wind_ms = 5.0
U = wind_ms * 196.85                # ft/min
```

### Slope

```@example usage8
# Convert from percent to tan(ϕ)
slope_percent = 30.0
tan_ϕ = slope_percent / 100.0       # 0.30

# Convert from degrees to tan(ϕ)
slope_degrees = 20.0
tan_ϕ = tan(deg2rad(slope_degrees)) # 0.364
```

### Fuel Load

```@example usage9
# Convert from kg/m² to lb/ft²
load_kgm2 = 1.0
w_o = load_kgm2 / 4.88              # lb/ft²

# Convert from tons/acre to lb/ft²
load_ton_acre = 2.0
w_o = load_ton_acre * 0.046         # lb/ft²
```

## Parametric Studies

### Varying a Single Parameter

```@example usage10
using WildlandFire
using OrdinaryDiffEqDefault

# Base parameters
base_params = Dict(
    h => 8000.0,
    S_T => 0.0555,
    S_e => 0.010,
    ρ_p => 32.0,
    σ => 3500.0,
    w_o => 0.138,
    δ => 1.0,
    M_x => 0.12,
    M_f => 0.05,
    tan_ϕ => 0.0,
)

# Create model
sys = Rothermel()

# Vary wind speed
wind_speeds = 0:100:1000  # ft/min
results = Float64[]

for U_val in wind_speeds
    params = [base_params..., U => U_val]
    prob = ODEProblem(sys, [], (0.0, 1.0), params)
    sol = solve(prob)
    push!(results, sol[R][end])
end

# Plot results
using Plots
plot(wind_speeds ./ 88.0, results,  # Convert wind to mph
     xlabel="Wind Speed (mph)",
     ylabel="Rate of Spread (ft/min)",
     title="Wind Effect on Fire Spread",
     linewidth=2,
     legend=false)
```

### Varying Multiple Parameters

```@example usage11
using WildlandFire
using OrdinaryDiffEqDefault

# Create parameter grid
moistures = 0.01:0.01:0.11
wind_speeds = [0.0, 352.0, 880.0]  # 0, 4, 10 mph
wind_labels = ["No wind", "4 mph", "10 mph"]

sys = Rothermel()
results = Dict()

for (wind, label) in zip(wind_speeds, wind_labels)
    rates = Float64[]
    for moisture in moistures
        params = [
            h => 8000.0, S_T => 0.0555, S_e => 0.010, ρ_p => 32.0,
            σ => 3500.0, w_o => 0.138, δ => 1.0, M_x => 0.12,
            M_f => moisture, U => wind, tan_ϕ => 0.0,
        ]
        prob = ODEProblem(sys, [], (0.0, 1.0), params)
        sol = solve(prob)
        push!(rates, sol[R][end])
    end
    results[label] = rates
end

# Plot multiple series
using Plots
plot(moistures .* 100, results["No wind"],
     label=wind_labels[1],
     xlabel="Fuel Moisture (%)",
     ylabel="Rate of Spread (ft/min)",
     title="Moisture and Wind Effects",
     linewidth=2)
plot!(moistures .* 100, results["4 mph"], label=wind_labels[2], linewidth=2)
plot!(moistures .* 100, results["10 mph"], label=wind_labels[3], linewidth=2)
```

## Working with Fuel Models

### Standard NFFL Fuel Models

Here are parameters for some standard Northern Forest Fire Laboratory (NFFL) fuel models:

```@example usage12
using WildlandFire

# Fuel Model 1: Short grass
fm1_params = [
    h => 8000.0, S_T => 0.0555, S_e => 0.010, ρ_p => 32.0,
    σ => 3500.0, w_o => 0.138, δ => 1.0, M_x => 0.12,
]

# Fuel Model 2: Timber grass and understory
fm2_params = [
    h => 8000.0, S_T => 0.0555, S_e => 0.010, ρ_p => 32.0,
    σ => 3000.0, w_o => 0.092, δ => 1.0, M_x => 0.15,
]

# Fuel Model 3: Tall grass
fm3_params = [
    h => 8000.0, S_T => 0.0555, S_e => 0.010, ρ_p => 32.0,
    σ => 1500.0, w_o => 0.138, δ => 2.5, M_x => 0.25,
]

# Fuel Model 4: Chaparral
fm4_params = [
    h => 8000.0, S_T => 0.0555, S_e => 0.010, ρ_p => 32.0,
    σ => 2000.0, w_o => 0.230, δ => 6.0, M_x => 0.20,
]

# Fuel Model 5: Brush
fm5_params = [
    h => 8000.0, S_T => 0.0555, S_e => 0.010, ρ_p => 32.0,
    σ => 1750.0, w_o => 0.046, δ => 2.0, M_x => 0.20,
]
```

### Creating Custom Fuel Models

```@example usage13
using WildlandFire

# Define custom fuel parameters
my_fuel_params = [
    # Particle properties (typically constant)
    h => 8000.0,          # Use standard heat content
    S_T => 0.0555,        # Standard mineral content
    S_e => 0.010,         # Standard effective mineral
    ρ_p => 32.0,          # Standard particle density

    # Fuel bed properties (vary by fuel type)
    σ => 2500.0,          # Medium-fine fuel
    w_o => 0.20,          # Moderate fuel load
    δ => 3.0,             # Medium depth fuel bed
    M_x => 0.18,          # Medium moisture of extinction

    # Environmental (vary by conditions)
    M_f => 0.08,          # Current moisture
    U => 440.0,           # Current wind (5 mph)
    tan_ϕ => 0.15,        # Current slope (15%)
]
```

## Batch Processing

### Multiple Scenarios

```@example usage14
using WildlandFire
using OrdinaryDiffEqDefault
using DataFrames

# Define scenarios
scenarios = [
    (name="Low danger", M_f=0.12, U=176.0, tan_ϕ=0.0),   # 2 mph
    (name="Moderate danger", M_f=0.08, U=352.0, tan_ϕ=0.0),  # 4 mph
    (name="High danger", M_f=0.04, U=880.0, tan_ϕ=0.0),      # 10 mph
    (name="Extreme danger", M_f=0.02, U=1320.0, tan_ϕ=0.2),  # 15 mph, 20% slope
]

# Common fuel parameters
base_fuel = Dict(
    h => 8000.0, S_T => 0.0555, S_e => 0.010, ρ_p => 32.0,
    σ => 3500.0, w_o => 0.138, δ => 1.0, M_x => 0.12,
)

sys = Rothermel()

# Run all scenarios
results = DataFrame(
    Scenario = String[],
    Moisture = Float64[],
    WindSpeed_mph = Float64[],
    Slope_pct = Float64[],
    ROS_ftmin = Float64[],
    ROS_mhr = Float64[],
)

for scenario in scenarios
    params = [
        base_fuel...,
        M_f => scenario.M_f,
        U => scenario.U,
        tan_ϕ => scenario.tan_ϕ,
    ]

    prob = ODEProblem(sys, [], (0.0, 1.0), params)
    sol = solve(prob)

    R_ftmin = sol[R][end]

    push!(results, (
        scenario.name,
        scenario.M_f * 100,          # Convert to percent
        scenario.U / 88.0,           # Convert to mph
        scenario.tan_ϕ * 100,        # Convert to percent
        R_ftmin,
        R_ftmin * 18.288,            # Convert to m/hr
    ))
end

println(results)
```

## Error Handling

### Checking for Valid Solutions

```julia
using WildlandFire
using OrdinaryDiffEqDefault
using SciMLBase

params = [...]  # your parameters here

prob = ODEProblem(sys, [], (0.0, 1.0), params)
sol = solve(prob)

# Check if solution succeeded
if SciMLBase.successful_retcode(sol)
    R_value = sol[R][end]
    println("Rate of spread: $R_value ft/min")
else
    println("Warning: Solution did not converge")
    println("Retcode: $(sol.retcode)")
end
```

### Validating Parameters

```@example usage16
function validate_params(params_dict)
    # Check that moisture is less than extinction
    if params_dict[M_f] >= params_dict[M_x]
        @warn "Fuel moisture ($(params_dict[M_f])) exceeds extinction moisture ($(params_dict[M_x]))"
        return false
    end

    # Check for negative values
    for (key, val) in params_dict
        if val < 0
            @warn "Negative value for $key: $val"
            return false
        end
    end

    # Check reasonable ranges
    if params_dict[σ] > 10000
        @warn "Very high SAV ratio: $(params_dict[σ]) ft²/ft³"
    end

    return true
end
```

## Advanced Topics

### Accessing the Equation System

For advanced users who want to inspect the underlying ModelingToolkit system:

```julia
using WildlandFire
using ModelingToolkit

# Get the unsimplified system
sys = Rothermel(simplify=false)

# Access equations
eqs = equations(sys)
println("Number of equations: $(length(eqs))")

# Access variables
vars = unknowns(sys)
println("Number of variables: $(length(vars))")

# Access parameters
params_sys = parameters(sys)
println("Number of parameters: $(length(params_sys))")
```

### Custom Solver Options

```@example usage18
using WildlandFire
using OrdinaryDiffEqDefault

# Create model and define parameters
sys = Rothermel()
params = [
    h => 8000.0, σ => 3500.0, w_o => 0.138, δ => 1.0,
    M_f => 0.05, U => 352.0, tan_ϕ => 0.0,
    S_T => 0.0555, S_e => 0.010, ρ_p => 32.0, M_x => 0.12
]

# Use default solver
prob = ODEProblem(sys, [], (0.0, 1.0), params)
sol = solve(prob)

# With solver options
sol = solve(prob, abstol=1e-8, reltol=1e-6)
```

## Tips and Best Practices

1. **Always check moisture < extinction**: Set M_f < M_x or fire won't spread
2. **Use consistent units**: Wind in ft/min, slope as tan(ϕ), not degrees/percent
3. **Validate results**: Check that rate of spread is reasonable for the fuel type
4. **Start simple**: Begin with no wind/slope, then add complexity
5. **Document assumptions**: Keep notes on parameter sources and assumptions
6. **Compare scenarios**: Use relative comparisons more than absolute values
7. **Validate locally**: If possible, validate against local fire behavior observations
