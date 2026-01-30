"""
    RothermelFireSpread(; name=:RothermelFireSpread)

Create a Rothermel surface fire spread model system.

The Rothermel model is a semi-empirical model for predicting the rate of spread of
surface fires in wildland fuels. It calculates the fire spread rate based on a heat
source/heat sink balance, where the heat source is the reaction intensity modified
by wind and slope factors, and the heat sink is the energy required to raise the
fuel to ignition temperature.

!!! note "Units"
    This implementation uses SI units. The original Rothermel equations were calibrated
    in US customary units, so all empirical coefficients have been converted to SI.
    Key unit conversions from the original formulation:
    - Length: 1 ft = 0.3048 m
    - Mass: 1 lb = 0.453592 kg
    - Energy: 1 Btu = 1055.06 J
    - Load: 1 lb/ft² = 4.88243 kg/m²
    - Density: 1 lb/ft³ = 16.0185 kg/m³

# Arguments
- `name`: System name (default: `:RothermelFireSpread`)

# Model Description

The fundamental equation is:
```
R = IR * ξ * (1 + φw + φs) / (ρb * ε * Qig)
```

Where:
- `R`: Rate of spread (m/s)
- `IR`: Reaction intensity (W/m²)
- `ξ`: Propagating flux ratio (dimensionless)
- `φw`: Wind factor (dimensionless)
- `φs`: Slope factor (dimensionless)
- `ρb`: Bulk density (kg/m³)
- `ε`: Effective heating number (dimensionless)
- `Qig`: Heat of preignition (J/kg)

# Reference

Andrews, Patricia L. 2018. The Rothermel surface fire spread model and associated
developments: A comprehensive explanation. Gen. Tech. Rep. RMRS-GTR-371. Fort Collins,
CO: U.S. Department of Agriculture, Forest Service, Rocky Mountain Research Station. 121 p.

# Example

```julia
using ModelingToolkit, WildlandFire, OrdinaryDiffEqDefault

sys = RothermelFireSpread()
compiled_sys = mtkcompile(sys)

# Set parameters for fuel model 1 (short grass) in SI units
# Note: This is an algebraic system, so we use NonlinearProblem
prob = NonlinearProblem(compiled_sys, [], [
    compiled_sys.σ => 11483.0,     # SAV ratio (1/m), converted from 3500 1/ft
    compiled_sys.w0 => 0.166,      # Fuel load (kg/m²), converted from 0.034 lb/ft²
    compiled_sys.δ => 0.3048,      # Fuel bed depth (m), converted from 1.0 ft
    compiled_sys.Mx => 0.12,       # Moisture of extinction (fraction)
    compiled_sys.Mf => 0.05,       # Fuel moisture content (fraction)
    compiled_sys.U => 2.235,       # Wind speed (m/s), converted from 5 mi/h
    compiled_sys.tanϕ => 0.0       # Slope (flat)
])
sol = solve(prob)
```
"""
@component function RothermelFireSpread(; name=:RothermelFireSpread)
    # Physical constants for fuel particles (Table 3) - converted to SI
    # Note: The Rothermel model is a semi-empirical model with coefficients calibrated
    # in US customary units. All coefficients have been converted to SI.
    @constants begin
        h_default = 18608000.0, [description = "Low heat content (8000 Btu/lb converted to J/kg)"]
        S_T_default = 0.0555, [description = "Total mineral content (dimensionless)"]
        S_e_default = 0.010, [description = "Effective mineral content (dimensionless)"]
        ρ_p_default = 512.6, [description = "Oven-dry particle density (32 lb/ft³ converted to kg/m³)"]
    end

    # Empirical coefficients converted to SI units
    @constants begin
        # β_op = 9.189 * σ^(-0.8189) (SI, original: 3.348 * σ^(-0.8189) in 1/ft)
        c_beta_op = 9.189, [description = "Optimum packing ratio coefficient (SI)"]

        # Γ_max = σ^1.5 / (176502 + 3.564*σ^1.5) (SI, units: 1/s)
        c_Gamma_denom1 = 176502.0, [description = "Γ_max denominator constant 1 (SI)"]
        c_Gamma_denom2 = 3.564, [description = "Γ_max denominator constant 2 (SI)"]

        # A_coeff = 343.2 * σ^(-0.7913) (SI, original: 133 * σ^(-0.7913) in 1/ft)
        c_A = 343.2, [description = "Coefficient A constant (SI)"]

        # Moisture damping polynomial coefficients (dimensionless)
        c_etaM_1 = 2.59, [description = "Moisture damping coefficient 1"]
        c_etaM_2 = 5.11, [description = "Moisture damping coefficient 2"]
        c_etaM_3 = 3.52, [description = "Moisture damping coefficient 3"]

        # Mineral damping: η_s = 0.174 * S_e^(-0.19)
        c_etas = 0.174, [description = "Mineral damping coefficient"]

        # ξ coefficients (SI)
        c_xi_1 = 192.0, [description = "ξ coefficient 1"]
        c_xi_2 = 0.07908, [description = "ξ coefficient 2 (SI)"]
        c_xi_3 = 0.792, [description = "ξ coefficient 3"]
        c_xi_4 = 0.376, [description = "ξ coefficient 4 (SI)"]

        # Wind coefficients (SI)
        c_C_1 = 7.47, [description = "Wind coefficient C constant 1"]
        c_C_2 = 0.0717, [description = "Wind coefficient C constant 2 (SI)"]
        c_B = 0.01317, [description = "Wind coefficient B constant (SI)"]
        c_E_1 = 0.715, [description = "Wind coefficient E constant 1"]
        c_E_2 = 1.094e-4, [description = "Wind coefficient E constant 2 (SI)"]

        # Slope factor coefficient
        c_phis = 5.275, [description = "Slope factor coefficient"]

        # ε = exp(-452.7/σ) (SI, original: exp(-138/σ) in 1/ft)
        c_eps = 452.7, [description = "Effective heating number coefficient (SI)"]

        # Qig = 581500 + 2595816*Mf (SI J/kg, original: 250 + 1116*Mf in Btu/lb)
        c_Qig_1 = 581500.0, [description = "Heat of preignition constant 1 (SI J/kg)"]
        c_Qig_2 = 2595816.0, [description = "Heat of preignition constant 2 (SI J/kg)"]

        # t_r = 7023.8/σ (SI s, original: 384/σ in min with σ in 1/ft)
        c_tr = 7023.8, [description = "Residence time coefficient (SI)"]

        # F_L = 0.003145 * IB^0.46 (SI m, original: 0.45 * IB^0.46 with IB in Btu/ft/s)
        c_FL = 0.003145, [description = "Flame length coefficient (SI)"]
    end

    # Input parameters - fuel array properties (SI units)
    @parameters begin
        σ, [description = "Surface-area-to-volume ratio (1/m)"]
        w0, [description = "Oven-dry fuel load (kg/m²)"]
        δ, [description = "Fuel bed depth (m)"]
        Mx, [description = "Dead fuel moisture of extinction (dimensionless)"]
        h = h_default, [description = "Low heat content (J/kg)"]
        S_T = S_T_default, [description = "Total mineral content (dimensionless)"]
        S_e = S_e_default, [description = "Effective mineral content (dimensionless)"]
        ρ_p = ρ_p_default, [description = "Oven-dry particle density (kg/m³)"]
    end

    # Environmental parameters (SI units)
    @parameters begin
        Mf, [description = "Fuel moisture content (dimensionless, dry weight basis)"]
        U, [description = "Wind velocity at midflame height (m/s)"]
        tanϕ, [description = "Slope steepness (dimensionless, rise/run)"]
    end

    # Intermediate variables - fuel bed calculations
    @variables begin
        wn(t), [description = "Net fuel load (kg/m²)"]
        ρb(t), [description = "Oven-dry bulk density (kg/m³)"]
        β(t), [description = "Packing ratio (dimensionless)"]
        β_op(t), [description = "Optimum packing ratio (dimensionless)"]
        β_ratio(t), [description = "Relative packing ratio β/β_op (dimensionless)"]
    end

    # Intermediate variables - heat source (numerator)
    @variables begin
        Γ_max(t), [description = "Maximum reaction velocity (1/s)"]
        A_coeff(t), [description = "Coefficient A for reaction velocity (dimensionless)"]
        Γ_prime(t), [description = "Optimum reaction velocity (1/s)"]
        rM(t), [description = "Moisture ratio Mf/Mx (dimensionless)"]
        η_M(t), [description = "Moisture damping coefficient (dimensionless)"]
        η_s(t), [description = "Mineral damping coefficient (dimensionless)"]
        IR(t), [description = "Reaction intensity (W/m²)"]
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
        Qig(t), [description = "Heat of preignition (J/kg)"]
    end

    # Output variables
    @variables begin
        R0(t), [description = "No-wind no-slope rate of spread (m/s)"]
        R(t), [description = "Rate of spread (m/s)"]
        t_r(t), [description = "Flame residence time (s)"]
        HA(t), [description = "Heat per unit area (J/m²)"]
        IB(t), [description = "Fireline intensity (Byram) (W/m)"]
        F_L(t), [description = "Flame length (Byram) (m)"]
    end

    eqs = [
        # Fuel bed calculations - Table 3
        wn ~ w0 * (1.0 - S_T),                              # Net fuel load
        ρb ~ w0 / δ,                                        # Bulk density
        β ~ ρb / ρ_p,                                       # Packing ratio
        β_op ~ c_beta_op * σ^(-0.8189),                     # Optimum packing ratio (SI)
        β_ratio ~ β / β_op,                                 # Relative packing ratio

        # Reaction intensity components - Table 3 (SI coefficients)
        Γ_max ~ σ^1.5 / (c_Gamma_denom1 + c_Gamma_denom2 * σ^1.5),  # Maximum reaction velocity (1/s)
        A_coeff ~ c_A * σ^(-0.7913),                        # Coefficient A (SI)
        Γ_prime ~ Γ_max * β_ratio^A_coeff * exp(A_coeff * (1.0 - β_ratio)), # Optimum reaction velocity

        # Damping coefficients - Table 3
        rM ~ min(Mf / Mx, 1.0),                             # Moisture ratio (capped at 1.0)
        η_M ~ 1.0 - c_etaM_1*rM + c_etaM_2*rM^2 - c_etaM_3*rM^3, # Moisture damping coefficient
        η_s ~ min(c_etas * S_e^(-0.19), 1.0),               # Mineral damping coefficient (capped at 1.0)

        # Reaction intensity - Table 3
        IR ~ Γ_prime * wn * h * η_M * η_s,                 # Reaction intensity (W/m²)

        # Propagating flux ratio - Table 3 (SI coefficients)
        ξ ~ (c_xi_1 + c_xi_2*σ)^(-1) * exp((c_xi_3 + c_xi_4*sqrt(σ)) * (β + 0.1)),

        # Wind factor coefficients - Table 3 (SI coefficients)
        C_coeff ~ c_C_1 * exp(-c_C_2 * σ^0.55),            # Coefficient C (SI)
        B_coeff ~ c_B * σ^0.54,                             # Coefficient B (SI)
        E_coeff ~ c_E_1 * exp(-c_E_2 * σ),                 # Coefficient E (SI)

        # Wind and slope factors - Table 3
        φw ~ C_coeff * U^B_coeff * β_ratio^(-E_coeff),     # Wind factor
        φs ~ c_phis * β^(-0.3) * tanϕ^2,                   # Slope factor

        # Heat sink - Table 3 (SI coefficients)
        ε ~ exp(-c_eps / σ),                               # Effective heating number
        Qig ~ c_Qig_1 + c_Qig_2 * Mf,                      # Heat of preignition (J/kg)

        # Rate of spread - Main equation
        R0 ~ (IR * ξ) / (ρb * ε * Qig),                    # No-wind no-slope spread rate (m/s)
        R ~ R0 * (1.0 + φw + φs),                          # Rate of spread with wind and slope (m/s)

        # Related models - Table 7 (SI coefficients)
        t_r ~ c_tr / σ,                                    # Residence time (s)
        HA ~ IR * t_r,                                     # Heat per unit area (J/m²)
        IB ~ HA * R,                                       # Fireline intensity (W/m)
        F_L ~ c_FL * IB^0.46,                              # Flame length (m)
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
    @constants begin
        # Transfer fraction coefficients (dimensionless)
        c_T_1 = -1.11, [description = "Transfer fraction coefficient 1"]
        c_T_2 = 1.33, [description = "Transfer fraction coefficient 2"]
    end

    @parameters begin
        w0_live_herb, [description = "Initial live herbaceous fuel load (kg/m²)"]
        Mf_live_herb, [description = "Live herbaceous fuel moisture content (dimensionless)"]
    end

    @variables begin
        T_fraction(t), [description = "Transfer fraction (dimensionless)"]
        w0_dead_herb(t), [description = "Dead herbaceous fuel load from transfer (kg/m²)"]
        w0_live_herb_remaining(t), [description = "Remaining live herbaceous fuel load (kg/m²)"]
    end

    eqs = [
        # Transfer fraction - clamped between 0 and 1
        T_fraction ~ max(0.0, min(1.0, c_T_1 * Mf_live_herb + c_T_2)),

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
    @constants begin
        # Coefficients for Mx_live calculation (dimensionless)
        c_Mx_1 = 2.9, [description = "Live Mx coefficient 1"]
        c_Mx_2 = 0.226, [description = "Live Mx coefficient 2"]
    end

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
        Mx_live ~ max(Mx_dead, c_Mx_1 * W_ratio * (1.0 - Mf_dead / Mx_dead) - c_Mx_2),
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
        UE(t), [description = "Effective midflame wind speed (m/s)"]
    end

    eqs = [
        # Combined wind and slope factor
        φE ~ φw + φs,

        # Effective midflame wind speed - Table 7
        # UE = [(φE * (β/βop)^E) / C]^(1/B)
        UE ~ (φE * β_ratio^E_coeff / C_coeff)^(1.0/B_coeff),
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
    # Wind limit coefficients converted to SI
    # Original: U_limit (ft/min) = 96.8 * IR^(1/3) where IR is in Btu/ft²/min
    # or U_limit = 0.9 * IR
    #
    # For corrected SI:
    # IR_us = IR_si / 189.27, U_limit_us = U_limit_si / 0.00508
    # U_limit_si / 0.00508 = 96.8 * (IR_si / 189.27)^(1/3)
    # U_limit_si = 0.00508 * 96.8 / 189.27^(1/3) = 0.00508 * 96.8 / 5.741 = 0.0857 * IR_si^(1/3)
    #
    # For original SI:
    # U_limit_si / 0.00508 = 0.9 * IR_si / 189.27
    # U_limit_si = 0.00508 * 0.9 / 189.27 * IR_si = 2.42e-5 * IR_si
    @constants begin
        c_Ulim_corr = 0.0857, [description = "Corrected wind limit coefficient (SI)"]
        c_Ulim_orig = 2.42e-5, [description = "Original wind limit coefficient (SI)"]
    end

    @parameters begin
        IR, [description = "Reaction intensity (W/m²)"]
    end

    @variables begin
        U_limit(t), [description = "Maximum effective wind speed (m/s)"]
    end

    if use_corrected
        # Corrected equation (Andrews et al. 2013) - SI
        eqs = [
            U_limit ~ c_Ulim_corr * IR^(1.0/3.0),
        ]
    else
        # Original equation - SI
        eqs = [
            U_limit ~ c_Ulim_orig * IR,
        ]
    end

    return System(eqs, t; name)
end


export RothermelFireSpread, DynamicFuelLoadTransfer, LiveFuelMoistureExtinction
export EffectiveMidflameWindSpeed, WindLimit
