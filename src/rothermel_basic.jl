"""
Rothermel Surface Fire Spread Model Implementation
Based on Rothermel (1972) and Albini (1976a)

This model calculates the rate of spread of a surface fire using the
Rothermel fire spread equation and related parameters.
"""

using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

"""
    Rothermel(; name=:rothermel, simplify=true)

Factory function to create a Rothermel surface fire spread model.

# Arguments
- `name::Symbol=:rothermel`: Name for the ODESystem
- `simplify::Bool=true`: Whether to apply structural_simplify to the system

# Returns
- `ODESystem`: A ModelingToolkit ODESystem representing the Rothermel fire spread model

# Description
This function creates a complete Rothermel surface fire spread model including all
equations from Rothermel (1972) and Albini (1976a). The model calculates the forward
rate of spread based on fuel properties, environmental conditions, and topography.

# Main Equation
The model calculates rate of spread (R) as:
```
R = I_R·ξ(1 + ϕ_w + ϕ_s) / (ρ_b·ε·Q_ig)
```

Where:
- I_R: Reaction intensity (heat release rate)
- ξ: Propagating flux ratio
- ϕ_w: Wind factor
- ϕ_s: Slope factor
- ρ_b: Bulk density
- ε: Effective heating number
- Q_ig: Heat of preignition

# Parameters
The model requires the following parameters:
- `h`: Low heat content (Btu/lb)
- `S_T`: Total mineral content (fraction)
- `S_e`: Effective mineral content (fraction)
- `ρ_p`: Particle density (lb/ft³)
- `σ`: Surface-area-to-volume ratio (ft²/ft³)
- `w_o`: Oven-dry fuel load (lb/ft²)
- `δ`: Fuel bed depth (ft)
- `M_x`: Moisture of extinction (fraction)
- `M_f`: Fuel moisture content (fraction)
- `U`: Wind velocity at midflame height (ft/min)
- `tan_ϕ`: Slope steepness (fraction)

# Example
```julia
using WildlandFire
using OrdinaryDiffEqDefault

# Create the model
sys = Rothermel()

# Define parameters
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
    U => 352.0,           # Wind speed (ft/min) ≈ 4 mph
    tan_ϕ => 0.0,         # Slope (flat)
]

# Solve
prob = ODEProblem(sys, [], (0.0, 1.0), params)
sol = solve(prob)

# Extract rate of spread
R_ftmin = sol[R][end]
println("Rate of spread: \$R_ftmin ft/min")
```

# References
1. Rothermel, R.C. 1972. A mathematical model for predicting fire spread in wildland fuels.
   USDA Forest Service Research Paper INT-115.
2. Albini, F.A. 1976a. Estimating wildfire behavior and effects.
   USDA Forest Service General Technical Report INT-30.
"""
function Rothermel(; name=:rothermel, simplify=true)
    @parameters begin
        # Fuel particle parameters
        h                    # Low heat content (Btu/lb)
        S_T                  # Total mineral content (fraction), generally 0.0555
        S_e                  # Effective mineral content (fraction), generally 0.010
        ρ_p                  # Oven-dry particle density (lb/ft³), generally 32 lb/ft³

        # Fuel array parameters
        σ                    # Surface-area-to-volume ratio (ft²/ft³)
        w_o                  # Oven-dry fuel load (lb/ft²)
        δ                    # Fuel bed depth (ft)
        M_x                  # Dead fuel moisture of extinction (fraction)

        # Environmental parameters
        M_f                  # Moisture content (fraction), dry weight basis
        U                    # Wind velocity at midflame height (ft/min)
        tan_ϕ                # Slope steepness, maximum (fraction) = vertical rise / horizontal distance

        # Constants for wind factor calculation
        C_const = 7.47        # Constant for wind factor
        B_const = 0.02526     # Constant for wind factor
        E_const = 0.715       # Constant for wind factor
    end

    @variables begin
        # Main outputs
        R(t)                 # Rate of spread (ft/min)

        # Intermediate calculations - Reaction intensity
        I_R(t)               # Reaction intensity (Btu/ft²·min)
        Γ_max(t)             # Maximum reaction velocity (min⁻¹)
        Γ(t)                 # Optimum reaction velocity (min⁻¹)
        A(t)                 # Constant for reaction velocity

        # Intermediate calculations - Packing ratios
        β(t)                 # Packing ratio (dimensionless)
        β_opt(t)             # Optimum packing ratio
        β_ratio(t)           # Relative packing ratio β/β_opt

        # Intermediate calculations - Bulk density
        ρ_b(t)               # Bulk density (lb/ft³)

        # Intermediate calculations - Net fuel load
        w_n(t)               # Net fuel load (lb/ft²)

        # Intermediate calculations - Damping coefficients
        η_M(t)               # Moisture damping coefficient
        η_s(t)               # Mineral damping coefficient

        # Intermediate calculations - Propagating flux ratio
        ξ(t)                 # Propagating flux ratio

        # Intermediate calculations - Wind and slope factors
        ϕ_w(t)               # Wind factor
        ϕ_s(t)               # Slope factor

        # Intermediate calculations - Heat source and sink
        I_R_ξ(t)             # No-wind, no-slope propagating flux (Btu/ft²·min)
        heat_source(t)       # Heat source (Btu/ft²/min)
        Q_ig(t)              # Heat of preignition (Btu/lb)
        ρ_b_ε_Q_ig(t)        # Heat sink (Btu/ft³)

        # Effective heating number
        ε(t)                 # Effective heating number (dimensionless)

        # No-wind, no-slope rate of spread
        R_0(t)               # No-wind, no-slope rate of spread (ft/min)
    end

    # Define the system of equations
    eqs = [
        # ============================================================================
        # FUEL BED CALCULATIONS
        # ============================================================================

        # Net fuel load (oven-dry fuel load with mineral reduction)
        # Equation from Albini (1976a): w_n = w_o(1 - S_T)
        w_n ~ w_o * (1 - S_T),

        # Oven-dry bulk density (fuel load divided by fuel bed depth)
        # ρ_b = w_o/δ
        ρ_b ~ w_o / δ,

        # Packing ratio (bulk density divided by particle density)
        # β = ρ_b/ρ_p
        β ~ ρ_b / ρ_p,

        # Optimum packing ratio (function of surface-area-to-volume ratio)
        # β_opt = 3.348σ^(-0.8189)
        β_opt ~ 3.348 * σ^(-0.8189),

        # Relative packing ratio
        β_ratio ~ β / β_opt,

        # ============================================================================
        # REACTION VELOCITY AND INTENSITY
        # ============================================================================

        # Maximum reaction velocity
        # Γ_max = σ^1.5 / (495 + 0.0594σ^1.5)
        Γ_max ~ σ^1.5 / (495 + 0.0594 * σ^1.5),

        # Constant A for optimum reaction velocity
        # A = 133σ^(-0.7913)
        A ~ 133 * σ^(-0.7913),

        # Optimum reaction velocity
        # Γ = Γ_max(β/β_opt)^A * exp[A(1 - β/β_opt)]
        Γ ~ Γ_max * β_ratio^A * exp(A * (1 - β_ratio)),

        # Mineral damping coefficient
        # η_s = 0.174S_e^(-0.19) (max = 1.0)
        η_s ~ min(0.174 * S_e^(-0.19), 1.0),

        # Moisture damping coefficient
        # η_M = 1 - 2.59(M_f/M_x) + 5.11(M_f/M_x)² - 3.52(M_f/M_x)³
        # with constraint η_M = M_f/M_x (max = 1.0)
        η_M ~ ifelse(M_f >= M_x, 0.0,
                     max(1 - 2.59*(M_f/M_x) + 5.11*(M_f/M_x)^2 - 3.52*(M_f/M_x)^3, 0.0)),

        # Reaction intensity (heat release rate per unit area of fire front)
        # I_R = Γ·w_n·h·η_M·η_s
        I_R ~ Γ * w_n * h * η_M * η_s,

        # ============================================================================
        # PROPAGATING FLUX RATIO
        # ============================================================================

        # Propagating flux ratio (proportion of reaction intensity that heats adjacent fuel)
        # ξ = (192 + 0.2595σ)^(-1) * exp[(0.792 + 0.681σ^0.5)(β + 0.1)]
        ξ ~ exp((0.792 + 0.681 * σ^0.5) * (β + 0.1)) / (192 + 0.2595 * σ),

        # ============================================================================
        # WIND AND SLOPE FACTORS
        # ============================================================================

        # Wind factor
        # ϕ_w = C * U^B * (β/β_opt)^(-E)
        # where C = 7.47exp(-0.133σ^0.55)
        #       B = 0.02526σ^0.54
        #       E = 0.715exp(-3.59 × 10^(-4)σ)
        ϕ_w ~ C_const * exp(-0.133 * σ^0.55) * U^(B_const * σ^0.54) *
              β_ratio^(-E_const * exp(-3.59e-4 * σ)),

        # Slope factor
        # ϕ_s = 5.275β^(-0.3)(tan ϕ)²
        ϕ_s ~ 5.275 * β^(-0.3) * tan_ϕ^2,

        # ============================================================================
        # HEAT SOURCE AND SINK
        # ============================================================================

        # Effective heating number (approaches unity for fine fuels, decreases for larger fuels)
        # ε = exp(-138/σ)
        ε ~ exp(-138 / σ),

        # Heat of preignition (energy required to ignite fuel per unit mass)
        # Q_ig = 250 + 1,116M_f
        Q_ig ~ 250 + 1116 * M_f,

        # Heat sink (product of bulk density, heating number, and heat of preignition)
        # ρ_b·ε·Q_ig
        ρ_b_ε_Q_ig ~ ρ_b * ε * Q_ig,

        # No-wind, no-slope propagating flux
        # I_R·ξ
        I_R_ξ ~ I_R * ξ,

        # Heat source (propagating flux with wind and slope effects)
        # I_R·ξ(1 + ϕ_w + ϕ_s)
        heat_source ~ I_R_ξ * (1 + ϕ_w + ϕ_s),

        # ============================================================================
        # RATE OF SPREAD
        # ============================================================================

        # No-wind, no-slope rate of spread
        # R_0 = I_R·ξ / (ρ_b·ε·Q_ig)
        R_0 ~ I_R_ξ / ρ_b_ε_Q_ig,

        # Final rate of spread equation (heat source divided by heat sink)
        # R = I_R·ξ(1 + ϕ_w + ϕ_s) / (ρ_b·ε·Q_ig)
        R ~ heat_source / ρ_b_ε_Q_ig,
    ]

    # Create the system with the specified name
    sys = ODESystem(eqs, t; name=name)

    # Apply structural simplification if requested
    if simplify
        return structural_simplify(sys)
    else
        return sys
    end
end

# Create default instance for backward compatibility
const rothermel_simplified = Rothermel(name=:rothermel, simplify=true)

# Export the symbolic variables for convenience (from the default instance)
# Note: These are defined in the global scope for ease of use
@parameters begin
    h, S_T, S_e, ρ_p, σ, w_o, δ, M_x, M_f, U, tan_ϕ
    C_const = 7.47
    B_const = 0.02526
    E_const = 0.715
end

@variables begin
    R(t), I_R(t), Γ_max(t), Γ(t), A(t)
    β(t), β_opt(t), β_ratio(t), ρ_b(t), w_n(t)
    η_M(t), η_s(t), ξ(t), ϕ_w(t), ϕ_s(t)
    I_R_ξ(t), heat_source(t), Q_ig(t), ρ_b_ε_Q_ig(t), ε(t), R_0(t)
end
