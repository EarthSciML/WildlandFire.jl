"""
    RothermelFireSpread(; name=:RothermelFireSpread)

Create a Rothermel surface fire spread model system.

The Rothermel model is a semi-empirical model for predicting the rate of spread of
surface fires in wildland fuels. It calculates the fire spread rate based on a heat
source/heat sink balance, where the heat source is the reaction intensity modified
by wind and slope factors, and the heat sink is the energy required to raise the
fuel to ignition temperature.

!!! note "Units"
    This model uses US customary units internally as the empirical coefficients were
    calibrated in this unit system. All inputs and outputs are expected in US customary
    units as documented for each parameter and variable. The key unit conversions are:
    - Load: 1 ton/acre = 0.0459137 lb/ft²
    - Wind speed: 1 mi/h = 88 ft/min
    - Energy: 1 Btu = 1055.06 J
    - Length: 1 ft = 0.3048 m
    - Mass: 1 lb = 0.453592 kg

# Arguments
- `name`: System name (default: `:RothermelFireSpread`)

# Model Description

The fundamental equation is:
```
R = IR * ξ * (1 + φw + φs) / (ρb * ε * Qig)
```

Where:
- `R`: Rate of spread (ft/min)
- `IR`: Reaction intensity (Btu/ft²/min)
- `ξ`: Propagating flux ratio (dimensionless)
- `φw`: Wind factor (dimensionless)
- `φs`: Slope factor (dimensionless)
- `ρb`: Bulk density (lb/ft³)
- `ε`: Effective heating number (dimensionless)
- `Qig`: Heat of preignition (Btu/lb)

# Reference

Andrews, Patricia L. 2018. The Rothermel surface fire spread model and associated
developments: A comprehensive explanation. Gen. Tech. Rep. RMRS-GTR-371. Fort Collins,
CO: U.S. Department of Agriculture, Forest Service, Rocky Mountain Research Station. 121 p.

# Example

```julia
using ModelingToolkit, WildlandFire, OrdinaryDiffEqDefault

sys = RothermelFireSpread()
compiled_sys = mtkcompile(sys)

# Set parameters for fuel model 1 (short grass)
# Note: This is an algebraic system, so we use NonlinearProblem
prob = NonlinearProblem(compiled_sys, [], [
    compiled_sys.σ => 3500.0,      # SAV ratio (1/ft)
    compiled_sys.w0 => 0.034,      # Fuel load (lb/ft²)
    compiled_sys.δ => 1.0,         # Fuel bed depth (ft)
    compiled_sys.Mx => 0.12,       # Moisture of extinction (fraction)
    compiled_sys.Mf => 0.05,       # Fuel moisture content (fraction)
    compiled_sys.U => 88.0 * 5,    # Wind speed 5 mi/h = 440 ft/min
    compiled_sys.tanϕ => 0.0       # Slope (flat)
])
sol = solve(prob)
```
"""
@component function RothermelFireSpread(; name=:RothermelFireSpread)
    # Physical constants for fuel particles (Table 3)
    # Note: Using values typical for wood-based fuels
    @constants begin
        h_default = 8000.0, [description = "Low heat content (Btu/lb)"]
        S_T_default = 0.0555, [description = "Total mineral content (dimensionless)"]
        S_e_default = 0.010, [description = "Effective mineral content (dimensionless)"]
        ρ_p_default = 32.0, [description = "Oven-dry particle density (lb/ft³)"]
    end

    # Input parameters - fuel array properties
    @parameters begin
        σ, [description = "Surface-area-to-volume ratio (ft²/ft³ = 1/ft)"]
        w0, [description = "Oven-dry fuel load (lb/ft²)"]
        δ, [description = "Fuel bed depth (ft)"]
        Mx, [description = "Dead fuel moisture of extinction (dimensionless)"]
        h = h_default, [description = "Low heat content (Btu/lb)"]
        S_T = S_T_default, [description = "Total mineral content (dimensionless)"]
        S_e = S_e_default, [description = "Effective mineral content (dimensionless)"]
        ρ_p = ρ_p_default, [description = "Oven-dry particle density (lb/ft³)"]
    end

    # Environmental parameters
    @parameters begin
        Mf, [description = "Fuel moisture content (dimensionless, dry weight basis)"]
        U, [description = "Wind velocity at midflame height (ft/min)"]
        tanϕ, [description = "Slope steepness (dimensionless, rise/run)"]
    end

    # Intermediate variables - fuel bed calculations
    @variables begin
        wn(t), [description = "Net fuel load (lb/ft²)"]
        ρb(t), [description = "Oven-dry bulk density (lb/ft³)"]
        β(t), [description = "Packing ratio (dimensionless)"]
        β_op(t), [description = "Optimum packing ratio (dimensionless)"]
        β_ratio(t), [description = "Relative packing ratio β/β_op (dimensionless)"]
    end

    # Intermediate variables - heat source (numerator)
    @variables begin
        Γ_max(t), [description = "Maximum reaction velocity (1/min)"]
        A_coeff(t), [description = "Coefficient A for reaction velocity (dimensionless)"]
        Γ_prime(t), [description = "Optimum reaction velocity (1/min)"]
        rM(t), [description = "Moisture ratio Mf/Mx (dimensionless)"]
        η_M(t), [description = "Moisture damping coefficient (dimensionless)"]
        η_s(t), [description = "Mineral damping coefficient (dimensionless)"]
        IR(t), [description = "Reaction intensity (Btu/ft²/min)"]
        ξ(t), [description = "Propagating flux ratio (dimensionless)"]
    end

    # Intermediate variables - wind and slope factors
    @variables begin
        C_coeff(t), [description = "Wind coefficient C (dimensionless)"]
        B_coeff(t), [description = "Wind coefficient B (dimensionless)"]
        E_coeff(t), [description = "Wind coefficient E (dimensionless)"]
        φw(t), [description = "Wind factor (dimensionless)"]
        φs(t), [description = "Slope factor (dimensionless)"]
    end

    # Intermediate variables - heat sink (denominator)
    @variables begin
        ε(t), [description = "Effective heating number (dimensionless)"]
        Qig(t), [description = "Heat of preignition (Btu/lb)"]
    end

    # Output variables
    @variables begin
        R0(t), [description = "No-wind no-slope rate of spread (ft/min)"]
        R(t), [description = "Rate of spread (ft/min)"]
        t_r(t), [description = "Flame residence time (min)"]
        HA(t), [description = "Heat per unit area (Btu/ft²)"]
        IB(t), [description = "Fireline intensity (Byram) (Btu/ft/s)"]
        F_L(t), [description = "Flame length (Byram) (ft)"]
    end

    eqs = [
        # Fuel bed calculations - Table 3
        wn ~ w0 * (1 - S_T),                              # Net fuel load
        ρb ~ w0 / δ,                                       # Bulk density
        β ~ ρb / ρ_p,                                      # Packing ratio
        β_op ~ 3.348 * σ^(-0.8189),                       # Optimum packing ratio
        β_ratio ~ β / β_op,                                # Relative packing ratio

        # Reaction intensity components - Table 3
        Γ_max ~ σ^1.5 / (495.0 + 0.0594 * σ^1.5),        # Maximum reaction velocity
        A_coeff ~ 133.0 * σ^(-0.7913),                    # Coefficient A
        Γ_prime ~ Γ_max * β_ratio^A_coeff * exp(A_coeff * (1 - β_ratio)), # Optimum reaction velocity

        # Damping coefficients - Table 3
        rM ~ min(Mf / Mx, 1.0),                           # Moisture ratio (capped at 1.0)
        η_M ~ 1.0 - 2.59*rM + 5.11*rM^2 - 3.52*rM^3,     # Moisture damping coefficient
        η_s ~ min(0.174 * S_e^(-0.19), 1.0),             # Mineral damping coefficient (capped at 1.0)

        # Reaction intensity - Table 3
        IR ~ Γ_prime * wn * h * η_M * η_s,               # Reaction intensity

        # Propagating flux ratio - Table 3
        ξ ~ (192.0 + 0.2595*σ)^(-1) * exp((0.792 + 0.681*sqrt(σ)) * (β + 0.1)),

        # Wind factor coefficients - Table 3
        C_coeff ~ 7.47 * exp(-0.133 * σ^0.55),           # Coefficient C
        B_coeff ~ 0.02526 * σ^0.54,                       # Coefficient B
        E_coeff ~ 0.715 * exp(-3.59e-4 * σ),             # Coefficient E

        # Wind and slope factors - Table 3
        φw ~ C_coeff * U^B_coeff * β_ratio^(-E_coeff),   # Wind factor
        φs ~ 5.275 * β^(-0.3) * tanϕ^2,                  # Slope factor

        # Heat sink - Table 3
        ε ~ exp(-138.0 / σ),                              # Effective heating number
        Qig ~ 250.0 + 1116.0 * Mf,                        # Heat of preignition

        # Rate of spread - Main equation
        R0 ~ (IR * ξ) / (ρb * ε * Qig),                  # No-wind no-slope spread rate
        R ~ R0 * (1 + φw + φs),                           # Rate of spread with wind and slope

        # Related models - Table 7
        t_r ~ 384.0 / σ,                                  # Residence time
        HA ~ IR * t_r,                                    # Heat per unit area
        IB ~ HA * R / 60.0,                               # Fireline intensity (Byram)
        F_L ~ 0.45 * IB^0.46,                             # Flame length (Byram)
    ]

    return System(eqs, t; name)
end


"""
    DynamicFuelLoadTransfer(; name=:DynamicFuelLoadTransfer)

Create a dynamic fuel model load transfer system.

Dynamic fuel models transfer load from live herbaceous fuel to dead herbaceous fuel
based on the live herbaceous fuel moisture content. This simulates the curing of
herbaceous fuels as they dry out during the fire season.

# Transfer Equation

The fraction of live herbaceous load transferred to dead herbaceous is:
```
T = -1.11 * Mf_live_herb + 1.33    (0 ≤ T ≤ 1.0)
```

Where `Mf_live_herb` is the live herbaceous fuel moisture content (fraction).

# Reference

Andrews, Patricia L. 2018. The Rothermel surface fire spread model and associated
developments: A comprehensive explanation. Gen. Tech. Rep. RMRS-GTR-371. Fort Collins,
CO: U.S. Department of Agriculture, Forest Service, Rocky Mountain Research Station. 121 p.
"""
@component function DynamicFuelLoadTransfer(; name=:DynamicFuelLoadTransfer)
    @parameters begin
        w0_live_herb, [description = "Initial live herbaceous fuel load (lb/ft²)"]
        Mf_live_herb, [description = "Live herbaceous fuel moisture content (dimensionless)"]
    end

    @variables begin
        T_fraction(t), [description = "Transfer fraction (dimensionless)"]
        w0_dead_herb(t), [description = "Dead herbaceous fuel load from transfer (lb/ft²)"]
        w0_live_herb_remaining(t), [description = "Remaining live herbaceous fuel load (lb/ft²)"]
    end

    eqs = [
        # Transfer fraction - clamped between 0 and 1
        T_fraction ~ max(0.0, min(1.0, -1.11 * Mf_live_herb + 1.33)),

        # Adjusted fuel loads
        w0_dead_herb ~ T_fraction * w0_live_herb,
        w0_live_herb_remaining ~ w0_live_herb - w0_dead_herb,
    ]

    return System(eqs, t; name)
end


"""
    LiveFuelMoistureExtinction(; name=:LiveFuelMoistureExtinction)

Create a system for calculating live fuel moisture of extinction.

The live fuel moisture of extinction is calculated based on the ratio of dead to live
fuel loading weighted by effective heating numbers, and the dead fuel moisture content.

# Reference

Andrews, Patricia L. 2018. The Rothermel surface fire spread model and associated
developments: A comprehensive explanation. Gen. Tech. Rep. RMRS-GTR-371. Fort Collins,
CO: U.S. Department of Agriculture, Forest Service, Rocky Mountain Research Station. 121 p.
"""
@component function LiveFuelMoistureExtinction(; name=:LiveFuelMoistureExtinction)
    @parameters begin
        Mx_dead, [description = "Dead fuel moisture of extinction (dimensionless)"]
        W_ratio, [description = "Dead-to-live effective fuel load ratio (dimensionless)"]
        Mf_dead, [description = "Fine dead fuel moisture content (dimensionless)"]
    end

    @variables begin
        Mx_live(t), [description = "Live fuel moisture of extinction (dimensionless)"]
    end

    eqs = [
        # Live fuel moisture of extinction - Eq. from Table 6b
        # Minimum value is the dead fuel moisture of extinction
        Mx_live ~ max(Mx_dead, 2.9 * W_ratio * (1 - Mf_dead / Mx_dead) - 0.226),
    ]

    return System(eqs, t; name)
end


"""
    EffectiveMidflameWindSpeed(; name=:EffectiveMidflameWindSpeed)

Create a system for calculating effective midflame wind speed.

The effective midflame wind speed combines the effects of actual wind and slope
into an equivalent wind speed for fire behavior calculations.

# Reference

Andrews, Patricia L. 2018. The Rothermel surface fire spread model and associated
developments: A comprehensive explanation. Gen. Tech. Rep. RMRS-GTR-371. Fort Collins,
CO: U.S. Department of Agriculture, Forest Service, Rocky Mountain Research Station. 121 p.
"""
@component function EffectiveMidflameWindSpeed(; name=:EffectiveMidflameWindSpeed)
    @parameters begin
        C_coeff, [description = "Wind coefficient C (dimensionless)"]
        B_coeff, [description = "Wind coefficient B (dimensionless)"]
        E_coeff, [description = "Wind coefficient E (dimensionless)"]
        β_ratio, [description = "Relative packing ratio (dimensionless)"]
        φw, [description = "Wind factor (dimensionless)"]
        φs, [description = "Slope factor (dimensionless)"]
    end

    @variables begin
        φE(t), [description = "Combined wind and slope factor (dimensionless)"]
        UE(t), [description = "Effective midflame wind speed (ft/min)"]
    end

    eqs = [
        # Combined wind and slope factor
        φE ~ φw + φs,

        # Effective midflame wind speed - Table 7
        # UE = [(φE * (β/βop)^E) / C]^(1/B)
        UE ~ (φE * β_ratio^E_coeff / C_coeff)^(1/B_coeff),
    ]

    return System(eqs, t; name)
end


"""
    WindLimit(; name=:WindLimit, use_corrected=true)

Create a system for calculating the maximum effective wind speed limit.

The wind limit represents the maximum effective wind speed that can influence
fire spread rate. Beyond this limit, the wind factor is capped.

# Arguments
- `use_corrected`: If true (default), use the corrected equation from Andrews et al. 2013.
                   If false, use the original equation.

# Reference

Andrews, Patricia L. 2018. The Rothermel surface fire spread model and associated
developments: A comprehensive explanation. Gen. Tech. Rep. RMRS-GTR-371. Fort Collins,
CO: U.S. Department of Agriculture, Forest Service, Rocky Mountain Research Station. 121 p.
"""
@component function WindLimit(; name=:WindLimit, use_corrected=true)
    @parameters begin
        IR, [description = "Reaction intensity (Btu/ft²/min)"]
    end

    @variables begin
        U_limit(t), [description = "Maximum effective wind speed (ft/min)"]
    end

    if use_corrected
        # Corrected equation (Andrews et al. 2013)
        eqs = [
            U_limit ~ 96.8 * IR^(1/3),
        ]
    else
        # Original equation
        eqs = [
            U_limit ~ 0.9 * IR,
        ]
    end

    return System(eqs, t; name)
end


export RothermelFireSpread, DynamicFuelLoadTransfer, LiveFuelMoistureExtinction
export EffectiveMidflameWindSpeed, WindLimit
