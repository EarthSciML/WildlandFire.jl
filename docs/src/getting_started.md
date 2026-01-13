# Getting Started

This guide will help you get up and running with WildlandFire.jl.

## Installation

Install the package using Julia's package manager:

```julia
using Pkg
Pkg.add("WildlandFire")
```

Or for the development version:

```julia
Pkg.add(url="https://github.com/EarthSciML/WildlandFire.jl")
```

## Required Dependencies

WildlandFire.jl has the following runtime dependencies:
- [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl) - Symbolic modeling framework
- [OrdinaryDiffEqDefault.jl](https://github.com/SciML/OrdinaryDiffEq.jl) - ODE solvers

For examples and plotting:
- [Plots.jl](https://github.com/JuliaPlots/Plots.jl) - Plotting functionality

## Your First Fire Simulation

Let's simulate a simple grass fire with moderate wind:

```@example getting_started
using WildlandFire
using OrdinaryDiffEqDefault

# Create the model
sys = Rothermel()

# Define fuel parameters for grass
params = [
    # Fuel particle properties
    h => 8000.0,          # Heat content (Btu/lb)
    S_T => 0.0555,        # Total mineral content (5.55%)
    S_e => 0.010,         # Effective mineral content (1%)
    ρ_p => 32.0,          # Particle density (lb/ft³)

    # Fuel bed properties
    σ => 3500.0,          # Surface-area-to-volume ratio (ft²/ft³) - fine grass
    w_o => 0.138,         # Fuel load (lb/ft²)
    δ => 1.0,             # Fuel bed depth (ft)
    M_x => 0.12,          # Moisture of extinction (12%)

    # Environmental conditions
    M_f => 0.05,          # Fuel moisture (5%)
    U => 352.0,           # Wind speed (ft/min) ≈ 4 mph
    tan_ϕ => 0.0,         # Slope (0 = flat)
]

# Solve the model
prob = ODEProblem(sys, [], (0.0, 1.0), params)
sol = solve(prob)

# Extract results
R_ftmin = sol[R][end]                  # Rate of spread in ft/min
R_mhr = R_ftmin * 18.288               # Convert to m/hr
R_chains_hr = R_ftmin * 1.829          # Convert to chains/hr

println("Fire Spread Results:")
println("  Rate of spread: $(round(R_ftmin, digits=2)) ft/min")
println("  Rate of spread: $(round(R_mhr, digits=2)) m/hr")
println("  Rate of spread: $(round(R_chains_hr, digits=2)) chains/hr")
```

## Understanding the Results

The model calculates that this grass fire spreads at about **71 ft/min** (or about 1.3 km/hr). This is a moderate spread rate typical of grass fires under light wind conditions.

### Key Factors Affecting Spread Rate

1. **Fine fuel (σ = 3500)**: Grass has a high surface-area-to-volume ratio, meaning it heats quickly
2. **Low moisture (5%)**: Dry fuel requires less energy to ignite
3. **Moderate wind (4 mph)**: Wind provides oxygen and preheats fuel ahead of the fire
4. **Flat terrain**: No slope to accelerate the fire uphill

## Next Steps

Now that you have a basic simulation running, you can:

1. **Explore different fuel types**: See [Fuel Types Examples](@ref)
2. **Analyze sensitivity**: Learn how parameters affect spread in [Sensitivity Analysis](@ref)
3. **Understand the model**: Read the [Model Description](@ref)
4. **Advanced usage**: Check the [Usage Guide](@ref)

## Common Units and Conversions

### Rate of Spread
- **1 ft/min** = 0.3048 m/min = 18.288 m/hr
- **1 ft/min** = 1.829 chains/hr (wildland fire standard)
- **1 m/hr** = 0.0547 ft/min

### Wind Speed
- **1 mph** ≈ 88 ft/min ≈ 1.47 ft/s
- **1 m/s** ≈ 196.85 ft/min ≈ 2.237 mph
- **Midflame height**: Wind should be measured at roughly half the flame height

### Fuel Load
- **1 lb/ft²** = 4.88 kg/m²
- **1 ton/acre** = 0.2242 kg/m² = 0.046 lb/ft²

### Slope
- **Slope (%)** = tan(ϕ) × 100
- Example: 30% slope → tan_ϕ = 0.30

## Typical Parameter Values

| Parameter | Symbol | Typical Range | Units |
|-----------|--------|---------------|-------|
| Heat content | h | 7000-9000 | Btu/lb |
| Total mineral content | S_T | 0.02-0.10 | fraction |
| Effective mineral content | S_e | 0.005-0.02 | fraction |
| Particle density | ρ_p | 28-40 | lb/ft³ |
| Surface-area-to-volume ratio | σ | 100-4000 | ft²/ft³ |
| Fuel load | w_o | 0.01-1.0 | lb/ft² |
| Fuel bed depth | δ | 0.1-10 | ft |
| Moisture of extinction | M_x | 0.12-0.40 | fraction |
| Fuel moisture | M_f | 0.01-0.30 | fraction |
| Wind speed | U | 0-2000 | ft/min |
| Slope | tan_ϕ | 0-1.0 | fraction |

## Troubleshooting

### Model returns very low rate of spread
- Check that fuel moisture (M_f) is less than moisture of extinction (M_x)
- Ensure fuel load (w_o) and fuel bed depth (δ) are reasonable
- Verify wind speed is in ft/min (not mph)

### Model returns very high rate of spread
- Check wind speed units (should be ft/min, not mph or m/s)
- Verify slope is as fraction (tan_ϕ), not percentage or degrees
- Ensure fuel parameters are realistic

### Solver fails to converge
- Try using the unsimplified system: `Rothermel(simplify=false)`
- Check for invalid parameter values (negative numbers, etc.)
- Ensure moisture < moisture of extinction

## Getting Help

If you encounter issues:

1. Check the [API Reference](@ref) for parameter definitions
2. Review the [Examples](@ref) for similar use cases
3. Search [GitHub Issues](https://github.com/EarthSciML/WildlandFire.jl/issues)
4. Ask a question in [GitHub Discussions](https://github.com/EarthSciML/WildlandFire.jl/discussions)
