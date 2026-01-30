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
prob = NonlinearProblem(compiled_sys, Dict(
    compiled_sys.σ => 11483.0,     # SAV ratio (1/m), converted from 3500 1/ft
    compiled_sys.w0 => 0.166,      # Fuel load (kg/m²), converted from 0.034 lb/ft²
    compiled_sys.δ => 0.3048,      # Fuel bed depth (m), converted from 1.0 ft
    compiled_sys.Mx => 0.12,       # Moisture of extinction (fraction)
    compiled_sys.Mf => 0.05,       # Fuel moisture content (fraction)
    compiled_sys.U => 2.235,       # Wind speed (m/s), converted from 5 mi/h
    compiled_sys.tanϕ => 0.0       # Slope (flat)
))
sol = solve(prob)
```
"""
@component function RothermelFireSpread(; name=:RothermelFireSpread)
    # Physical constants for fuel particles (Table 3, Andrews 2018) - converted to SI
    # Note: The Rothermel model is a semi-empirical model with coefficients calibrated
    # in US customary units. All coefficients have been converted to SI.
    @constants begin
        h_default = 18608000.0, [description = "Low heat content (8000 Btu/lb converted to J/kg)", unit = u"J/kg"]
        S_T_default = 0.0555, [description = "Total mineral content (dimensionless)", unit = u"1"]
        S_e_default = 0.010, [description = "Effective mineral content (dimensionless)", unit = u"1"]
        ρ_p_default = 512.6, [description = "Oven-dry particle density (32 lb/ft³ converted to kg/m³)", unit = u"kg/m^3"]
    end

    # Empirical coefficients converted to SI units (Table 3, Andrews 2018)
    # Note: Many of these coefficients arise from empirical calibration and have
    # complex dimensional relationships. We use dimensionless unit annotations
    # where appropriate since these are empirical calibration constants.
    @constants begin
        # β_op = 3.348 * σ^(-0.8189) (US), converted to SI: 8.858 * σ^(-0.8189)
        # Conversion: 3.348 * ft_to_m^(-0.8189) = 3.348 * 0.3048^(-0.8189) = 8.858
        c_beta_op = 8.858, [description = "Optimum packing ratio coefficient (SI, empirical)", unit = u"1"]

        # Γ_max = σ^1.5 / (495 + 0.0594*σ^1.5) (US), converted to SI: σ^1.5 / (176502 + 3.564*σ^1.5)
        # Note: These coefficients are empirical and unit-absorbing
        c_Gamma_denom1 = 176502.0, [description = "Γ_max denominator constant 1 (SI, empirical)", unit = u"1"]
        c_Gamma_denom2 = 3.564, [description = "Γ_max denominator constant 2 (SI, empirical)", unit = u"1"]

        # A = 133 * σ^(-0.7913) (US), converted to SI: 340.5 * σ^(-0.7913)
        # Conversion: 133 * ft_to_m^(-0.7913) = 133 * 0.3048^(-0.7913) = 340.5
        c_A = 340.5, [description = "Coefficient A constant (SI, empirical)", unit = u"1"]

        # Moisture damping polynomial coefficients (dimensionless) - Table 3
        c_etaM_1 = 2.59, [description = "Moisture damping coefficient 1 (dimensionless)", unit = u"1"]
        c_etaM_2 = 5.11, [description = "Moisture damping coefficient 2 (dimensionless)", unit = u"1"]
        c_etaM_3 = 3.52, [description = "Moisture damping coefficient 3 (dimensionless)", unit = u"1"]

        # Mineral damping: η_s = 0.174 * S_e^(-0.19) - Table 3
        c_etas = 0.174, [description = "Mineral damping coefficient (empirical)", unit = u"1"]

        # ξ coefficients (SI) - Table 3
        # ξ = (192 + 0.2595σ)^(-1) * exp[(0.792 + 0.681σ^0.5)(β + 0.1)]
        # SI: (192 + 0.07908σ)^(-1) * exp[(0.792 + 0.376σ^0.5)(β + 0.1)]
        c_xi_1 = 192.0, [description = "ξ coefficient 1 (empirical)", unit = u"1"]
        c_xi_2 = 0.07908, [description = "ξ coefficient 2 (SI, empirical)", unit = u"1"]
        c_xi_3 = 0.792, [description = "ξ coefficient 3 (empirical)", unit = u"1"]
        c_xi_4 = 0.376, [description = "ξ coefficient 4 (SI, empirical)", unit = u"1"]

        # Wind coefficients (SI) - Table 3
        # C = 7.47 * exp(-0.133 * σ^0.55) (US), SI: 7.47 * exp(-0.0692 * σ^0.55)
        # Conversion: 0.133 * ft_to_m^0.55 = 0.133 * 0.3048^0.55 = 0.0692
        c_C_1 = 7.47, [description = "Wind coefficient C constant 1 (empirical)", unit = u"1"]
        c_C_2 = 0.0692, [description = "Wind coefficient C constant 2 (SI, empirical)", unit = u"1"]
        # B = 0.02526 * σ^0.54 (US), SI: 0.01317 * σ^0.54
        c_B = 0.01317, [description = "Wind coefficient B constant (SI, empirical)", unit = u"1"]
        # E = 0.715 * exp(-3.59e-4 * σ) (US), SI: 0.715 * exp(-1.094e-4 * σ)
        c_E_1 = 0.715, [description = "Wind coefficient E constant 1 (empirical)", unit = u"1"]
        c_E_2 = 1.094e-4, [description = "Wind coefficient E constant 2 (SI, empirical)", unit = u"1"]

        # Slope factor coefficient - Table 3
        c_phis = 5.275, [description = "Slope factor coefficient (empirical)", unit = u"1"]

        # ε = exp(-138/σ) (US), SI: exp(-452.7/σ) - Table 3
        c_eps = 452.7, [description = "Effective heating number coefficient (SI, empirical)", unit = u"1"]

        # Qig = 250 + 1116*Mf (US Btu/lb), SI: 581500 + 2595816*Mf (J/kg) - Table 3
        c_Qig_1 = 581500.0, [description = "Heat of preignition constant 1 (SI)", unit = u"J/kg"]
        c_Qig_2 = 2595816.0, [description = "Heat of preignition constant 2 (SI)", unit = u"J/kg"]

        # t_r = 384/σ (US min with σ in 1/ft), SI: 75590.6/σ (s) - Table 7
        # Conversion: 384 * 60 / ft_to_m = 384 * 60 / 0.3048 = 75590.6
        c_tr = 75590.6, [description = "Residence time coefficient (SI, empirical)", unit = u"1"]

        # F_L = 0.45 * IB^0.46 (US ft with IB in Btu/ft/s), SI: 0.00323 * IB^0.46 - Table 7
        # Conversion: 0.45 * ft_to_m / (Btu_to_J / ft_to_m)^0.46 = 0.45 * 0.3048 / 3461.5^0.46 = 0.00323
        c_FL = 0.00323, [description = "Flame length coefficient (SI, empirical)", unit = u"1"]
    end

    # Input parameters - fuel array properties (SI units)
    @parameters begin
        σ, [description = "Surface-area-to-volume ratio", unit = u"1/m"]
        w0, [description = "Oven-dry fuel load", unit = u"kg/m^2"]
        δ, [description = "Fuel bed depth", unit = u"m"]
        Mx, [description = "Dead fuel moisture of extinction (dimensionless)", unit = u"1"]
        h = h_default, [description = "Low heat content", unit = u"J/kg"]
        S_T = S_T_default, [description = "Total mineral content (dimensionless)", unit = u"1"]
        S_e = S_e_default, [description = "Effective mineral content (dimensionless)", unit = u"1"]
        ρ_p = ρ_p_default, [description = "Oven-dry particle density", unit = u"kg/m^3"]
    end

    # Environmental parameters (SI units)
    @parameters begin
        Mf, [description = "Fuel moisture content (dimensionless, dry weight basis)", unit = u"1"]
        U, [description = "Wind velocity at midflame height", unit = u"m/s"]
        tanϕ, [description = "Slope steepness (dimensionless, rise/run)", unit = u"1"]
    end

    # Intermediate variables - fuel bed calculations
    @variables begin
        wn(t), [description = "Net fuel load", unit = u"kg/m^2"]
        ρb(t), [description = "Oven-dry bulk density", unit = u"kg/m^3"]
        β(t), [description = "Packing ratio (dimensionless)", unit = u"1"]
        β_op(t), [description = "Optimum packing ratio (dimensionless)", unit = u"1"]
        β_ratio(t), [description = "Relative packing ratio β/β_op (dimensionless)", unit = u"1"]
    end

    # Intermediate variables - heat source (numerator)
    @variables begin
        Γ_max(t), [description = "Maximum reaction velocity", unit = u"1/s"]
        A_coeff(t), [description = "Coefficient A for reaction velocity (dimensionless)", unit = u"1"]
        Γ_prime(t), [description = "Optimum reaction velocity", unit = u"1/s"]
        rM(t), [description = "Moisture ratio Mf/Mx (dimensionless)", unit = u"1"]
        η_M(t), [description = "Moisture damping coefficient (dimensionless)", unit = u"1"]
        η_s(t), [description = "Mineral damping coefficient (dimensionless)", unit = u"1"]
        IR(t), [description = "Reaction intensity", unit = u"W/m^2"]
        ξ(t), [description = "Propagating flux ratio (dimensionless)", unit = u"1"]
    end

    # Intermediate variables - wind and slope factors
    @variables begin
        C_coeff(t), [description = "Wind coefficient C (dimensionless)", unit = u"1"]
        B_coeff(t), [description = "Wind coefficient B (dimensionless)", unit = u"1"]
        E_coeff(t), [description = "Wind coefficient E (dimensionless)", unit = u"1"]
        φw(t), [description = "Wind factor (dimensionless)", unit = u"1"]
        φs(t), [description = "Slope factor (dimensionless)", unit = u"1"]
    end

    # Intermediate variables - heat sink (denominator)
    @variables begin
        ε(t), [description = "Effective heating number (dimensionless)", unit = u"1"]
        Qig(t), [description = "Heat of preignition", unit = u"J/kg"]
    end

    # Output variables
    @variables begin
        R0(t), [description = "No-wind no-slope rate of spread", unit = u"m/s"]
        R(t), [description = "Rate of spread", unit = u"m/s"]
        t_r(t), [description = "Flame residence time", unit = u"s"]
        HA(t), [description = "Heat per unit area", unit = u"J/m^2"]
        IB(t), [description = "Fireline intensity (Byram)", unit = u"W/m"]
        F_L(t), [description = "Flame length (Byram)", unit = u"m"]
    end

    eqs = [
        # Fuel bed calculations - Table 3, Andrews (2018)
        wn ~ w0 * (1.0 - S_T),                              # Eq. wn, Net fuel load
        ρb ~ w0 / δ,                                        # Eq. ρb, Bulk density
        β ~ ρb / ρ_p,                                       # Eq. β, Packing ratio
        β_op ~ c_beta_op * σ^(-0.8189),                     # Eq. βop, Optimum packing ratio (SI)
        β_ratio ~ β / β_op,                                 # Eq. β/βop, Relative packing ratio

        # Reaction intensity components - Table 3, Andrews (2018) (SI coefficients)
        Γ_max ~ σ^1.5 / (c_Gamma_denom1 + c_Gamma_denom2 * σ^1.5),  # Eq. Γmax, Maximum reaction velocity
        A_coeff ~ c_A * σ^(-0.7913),                        # Eq. A, Coefficient A (SI)
        Γ_prime ~ Γ_max * β_ratio^A_coeff * exp(A_coeff * (1.0 - β_ratio)), # Eq. Γ', Optimum reaction velocity

        # Damping coefficients - Table 3, Andrews (2018)
        rM ~ min(Mf / Mx, 1.0),                             # Eq. rM, Moisture ratio (capped at 1.0)
        η_M ~ 1.0 - c_etaM_1*rM + c_etaM_2*rM^2 - c_etaM_3*rM^3, # Eq. ηM, Moisture damping coefficient
        η_s ~ min(c_etas * S_e^(-0.19), 1.0),               # Eq. ηs, Mineral damping coefficient (capped at 1.0)

        # Reaction intensity - Table 3, Andrews (2018)
        IR ~ Γ_prime * wn * h * η_M * η_s,                 # Eq. IR, Reaction intensity

        # Propagating flux ratio - Table 3, Andrews (2018) (SI coefficients)
        ξ ~ (c_xi_1 + c_xi_2*σ)^(-1) * exp((c_xi_3 + c_xi_4*sqrt(σ)) * (β + 0.1)), # Eq. ξ

        # Wind factor coefficients - Table 3, Andrews (2018) (SI coefficients)
        C_coeff ~ c_C_1 * exp(-c_C_2 * σ^0.55),            # Eq. C, Wind coefficient C (SI)
        B_coeff ~ c_B * σ^0.54,                             # Eq. B, Wind coefficient B (SI)
        E_coeff ~ c_E_1 * exp(-c_E_2 * σ),                 # Eq. E, Wind coefficient E (SI)

        # Wind and slope factors - Table 3, Andrews (2018)
        φw ~ C_coeff * U^B_coeff * β_ratio^(-E_coeff),     # Eq. φw, Wind factor
        φs ~ c_phis * β^(-0.3) * tanϕ^2,                   # Eq. φs, Slope factor

        # Heat sink - Table 3, Andrews (2018) (SI coefficients)
        ε ~ exp(-c_eps / σ),                               # Eq. ε, Effective heating number
        Qig ~ c_Qig_1 + c_Qig_2 * Mf,                      # Eq. Qig, Heat of preignition

        # Rate of spread - Eq. 1, Table 3, Andrews (2018)
        R0 ~ (IR * ξ) / (ρb * ε * Qig),                    # Eq. R (no-wind no-slope)
        R ~ R0 * (1.0 + φw + φs),                          # Eq. R (with wind and slope)

        # Related models - Table 7, Andrews (2018) (SI coefficients)
        t_r ~ c_tr / σ,                                    # Eq. tr, Residence time
        HA ~ IR * t_r,                                     # Eq. HA, Heat per unit area
        IB ~ HA * R,                                       # Eq. IB, Fireline intensity (Byram)
        F_L ~ c_FL * IB^0.46,                              # Eq. FL, Flame length (Byram)
    ]

    # Note: Unit checking is disabled for this system because the Rothermel model
    # is semi-empirical with power-law relationships that don't have consistent
    # dimensional analysis (e.g., σ^(-0.8189), σ^0.54, etc.)
    return System(eqs, t; name, checks=false)
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
        # Transfer fraction coefficients (dimensionless) - Table 6c, Andrews (2018)
        c_T_1 = -1.11, [description = "Transfer fraction coefficient 1 (dimensionless)", unit = u"1"]
        c_T_2 = 1.33, [description = "Transfer fraction coefficient 2 (dimensionless)", unit = u"1"]
    end

    @parameters begin
        w0_live_herb, [description = "Initial live herbaceous fuel load", unit = u"kg/m^2"]
        Mf_live_herb, [description = "Live herbaceous fuel moisture content (dimensionless)", unit = u"1"]
    end

    @variables begin
        T_fraction(t), [description = "Transfer fraction (dimensionless)", unit = u"1"]
        w0_dead_herb(t), [description = "Dead herbaceous fuel load from transfer", unit = u"kg/m^2"]
        w0_live_herb_remaining(t), [description = "Remaining live herbaceous fuel load", unit = u"kg/m^2"]
    end

    eqs = [
        # Transfer fraction - clamped between 0 and 1 - Table 6c, Andrews (2018)
        T_fraction ~ max(0.0, min(1.0, c_T_1 * Mf_live_herb + c_T_2)),  # Eq. T

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
        # Coefficients for Mx_live calculation (dimensionless) - Table 6b, Andrews (2018)
        c_Mx_1 = 2.9, [description = "Live Mx coefficient 1 (dimensionless)", unit = u"1"]
        c_Mx_2 = 0.226, [description = "Live Mx coefficient 2 (dimensionless)", unit = u"1"]
    end

    @parameters begin
        Mx_dead, [description = "Dead fuel moisture of extinction (dimensionless)", unit = u"1"]
        W_ratio, [description = "Dead-to-live effective fuel load ratio (dimensionless)", unit = u"1"]
        Mf_dead, [description = "Fine dead fuel moisture content (dimensionless)", unit = u"1"]
    end

    @variables begin
        Mx_live(t), [description = "Live fuel moisture of extinction (dimensionless)", unit = u"1"]
    end

    eqs = [
        # Live fuel moisture of extinction - Eq. Mxlive, Table 6b, Andrews (2018)
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
        C_coeff, [description = "Wind coefficient C (dimensionless)", unit = u"1"]
        B_coeff, [description = "Wind coefficient B (dimensionless)", unit = u"1"]
        E_coeff, [description = "Wind coefficient E (dimensionless)", unit = u"1"]
        β_ratio, [description = "Relative packing ratio (dimensionless)", unit = u"1"]
        φw, [description = "Wind factor (dimensionless)", unit = u"1"]
        φs, [description = "Slope factor (dimensionless)", unit = u"1"]
    end

    @variables begin
        φE(t), [description = "Combined wind and slope factor (dimensionless)", unit = u"1"]
        UE(t), [description = "Effective midflame wind speed", unit = u"m/s"]
    end

    eqs = [
        # Combined wind and slope factor
        φE ~ φw + φs,

        # Effective midflame wind speed - Eq. UE, Table 7, Andrews (2018)
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
    # Wind limit coefficients converted to SI - Table 7, Andrews (2018)
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
        c_Ulim_corr = 0.0857, [description = "Corrected wind limit coefficient (SI)", unit = u"m/s/(W/m^2)^(1/3)"]
        c_Ulim_orig = 2.42e-5, [description = "Original wind limit coefficient (SI)", unit = u"m/s/(W/m^2)"]
    end

    @parameters begin
        IR, [description = "Reaction intensity", unit = u"W/m^2"]
    end

    @variables begin
        U_limit(t), [description = "Maximum effective wind speed", unit = u"m/s"]
    end

    if use_corrected
        # Corrected equation (Andrews et al. 2013) - Eq. Ulimit, Table 7, Andrews (2018)
        eqs = [
            U_limit ~ c_Ulim_corr * IR^(1.0/3.0),
        ]
    else
        # Original equation - Eq. Ulimit (original), Table 7, Andrews (2018)
        eqs = [
            U_limit ~ c_Ulim_orig * IR,
        ]
    end

    return System(eqs, t; name)
end


export RothermelFireSpread, DynamicFuelLoadTransfer, LiveFuelMoistureExtinction
export EffectiveMidflameWindSpeed, WindLimit
