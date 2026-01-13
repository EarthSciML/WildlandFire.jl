"""
# WildlandFire.jl

A Julia package for wildland fire behavior modeling using ModelingToolkit.jl.

## Main Features

- **Rothermel Surface Fire Spread Model**: Complete implementation of the Rothermel (1972)
  fire spread equations for both homogeneous and heterogeneous fuel beds
- **Symbolic modeling**: Built on ModelingToolkit.jl for automatic differentiation and optimization
- **Validated physics**: All models tested against expected fire behavior
- **Easy to use**: Simple API for common fire behavior calculations

## Modules

- `rothermel_basic`: Basic Rothermel model for homogeneous fuel beds
- `rothermel_heterogeneous`: Extended model for multiple fuel size classes

## Usage

```julia
using WildlandFire
using OrdinaryDiffEqDefault

# Get the basic Rothermel model
sys = rothermel_simplified

# Define fuel and environmental parameters
params = [
    h => 8000.0,      # Heat content (Btu/lb)
    σ => 3500.0,      # Surface-area-to-volume ratio (ft²/ft³)
    w_o => 0.138,     # Fuel load (lb/ft²)
    δ => 1.0,         # Fuel bed depth (ft)
    M_f => 0.05,      # Fuel moisture (fraction)
    U => 352.0,       # Wind speed (ft/min) ≈ 4 mph
    tan_ϕ => 0.0,     # Slope (fraction)
    # ... other parameters
]

# Solve for rate of spread
prob = ODEProblem(sys, [], (0.0, 1.0), params)
sol = solve(prob)

# Extract rate of spread
R_ftmin = sol[R][end]  # ft/min
```

## References

1. Rothermel, R.C. 1972. A mathematical model for predicting fire spread in wildland fuels.
   USDA Forest Service Research Paper INT-115.

2. Albini, F.A. 1976a. Estimating wildfire behavior and effects.
   USDA Forest Service General Technical Report INT-30.

3. Andrews, P.L. 2018. The Rothermel surface fire spread model and associated developments:
   A comprehensive explanation. USDA Forest Service General Technical Report RMRS-GTR-371.

## License

MIT License - Copyright (c) 2026 EarthSciML authors and contributors
"""
module WildlandFire

using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

# Include model implementations
include("rothermel_basic.jl")
include("rothermel_heterogeneous.jl")

# Export the main model systems and factory functions
export Rothermel  # Factory function for creating Rothermel models
export rothermel_simplified  # Pre-built simplified model (backward compatibility)
export create_heterogeneous_rothermel_model

# Export all symbolic variables from basic model for user convenience
export R, I_R, Γ_max, Γ, A, β, β_opt, β_ratio, ρ_b, w_n
export η_M, η_s, ξ, ϕ_w, ϕ_s, I_R_ξ, heat_source, Q_ig, ρ_b_ε_Q_ig, ε, R_0

# Export parameters
export h, S_T, S_e, ρ_p, σ, w_o, δ, M_x, M_f, U, tan_ϕ
export C_const, B_const, E_const

end # module
