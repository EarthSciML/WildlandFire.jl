"""
    Clark1996FireSpread(; name=:Clark1996FireSpread)

Create a coupled atmosphere-fire model component implementing the fire spread rate
and fuel consumption equations from Clark et al. (1996).

This implements a simple dry eucalyptus forest fire model designed to study the coupling
between fire dynamics and atmospheric circulations. The model uses the McArthur empirical
fire spread formula (Noble et al., 1980) to compute fire spread rates as a function of
wind speed, and tracks fuel consumption for four fuel types (litter, trash, scrub, and
canopy) with corresponding heat flux outputs.

# Model Description

The fire spread rate follows the McArthur formula (Eq. 9):

```math
S_f = S_a \\cdot \\exp(0.08424 \\cdot |\\vec{V}_A|)
```

where ``S_a = 0.18`` m/s is the backing fire speed and ``|\\vec{V}_A|`` is the wind speed
at 15 m above ground level.

The burn rates are modulated by a ratio ``B_{ratio}`` (Eq. 8) that accounts for oxygen
limitation at low wind speeds:

```math
B_{ratio} = \\sqrt{\\frac{|\\vec{V}_A| + 1}{|\\vec{V}_A| + 4}}
```

Each fuel type (litter, trash, scrub, canopy) is consumed at its nominal burn rate
multiplied by ``B_{ratio}``. The sensible heat flux is computed from the combustion
coefficient ``H_c = 1.7 \\times 10^7`` J/kg applied to the total fuel burn rate, with 3%
of the sensible heat going to evaporate moisture in the ground fuel. The latent heat flux
accounts for both the evaporated fuel moisture and the water vapor produced by cellulose
combustion (56% of dry fuel mass).

# Canopy Ignition

The canopy automatically ignites when the cumulative ground heat flux exceeds a threshold
of 170 kJ/m² (as stated on p. 180 of Clark et al., 1996). Before ignition, only ground
fuels (litter, trash, scrub) burn. The model tracks the cumulative heat energy released
and triggers canopy burning once the threshold is exceeded.

# Arguments
- `name`: System name (default `:Clark1996FireSpread`)

# Reference

Clark, T.L., Jenkins, M.A., Coen, J.L., and Packham, D.R. (1996). A Coupled
Atmosphere-Fire Model: Role of the Convective Froude Number and Dynamic Fingering at the
Fireline. *Int. J. Wildland Fire*, 6(4), 177-190.

Noble, I.R., Bary, G.A.V., and Gill, A.M. (1980). McArthur's fire-danger meters expressed
as equations. *Aust. J. Ecol.*, 5, 201-203.

Walker, J. (1981). Fuel dynamics in Australian vegetation. In *Fire in the Australian
Biota*, ed. A.M. Gill, R.H. Groves, and I.R. Noble. Australian Academy of Science, 582 pp.

# Example

```julia
using WildlandFire, ModelingToolkit, DynamicQuantities, OrdinaryDiffEqDefault
using ModelingToolkit: t, D

sys = Clark1996FireSpread()
compiled = mtkcompile(sys)

# Set initial fuel masses and run at wind speed of 3 m/s
prob = ODEProblem(compiled,
    [compiled.M_litter => 2.0, compiled.M_trash => 0.5,
     compiled.M_scrub => 0.2, compiled.M_canopy => 1.2,
     compiled.Q_cumulative => 0.0],
    (0.0, 120.0),
    [compiled.V_A => 3.0])
sol = solve(prob)

# Canopy ignition occurs automatically when cumulative heat flux reaches 170 kJ/m²
```
"""
@component function Clark1996FireSpread(; name = :Clark1996FireSpread)
    @constants begin
        # Fire spread rate constants — Eq. 9, Clark et al. (1996)
        S_a = 0.18, [description = "Backing/base fire spread rate (McArthur)", unit = u"m/s"]
        k_spread = 0.08424, [description = "McArthur spread rate coefficient", unit = u"s/m"]

        # Fuel parameters — Walker (1981), p. 180 of Clark et al. (1996)
        R_litter = 0.04, [description = "Litter nominal burn rate", unit = u"kg/(m^2*s)"]
        R_trash = 0.005, [description = "Trash nominal burn rate", unit = u"kg/(m^2*s)"]
        R_scrub = 0.004, [description = "Scrub nominal burn rate", unit = u"kg/(m^2*s)"]
        R_canopy = 0.02, [description = "Dry canopy fuel nominal burn rate", unit = u"kg/(m^2*s)"]

        # Combustion constants — p. 180 of Clark et al. (1996)
        H_c = 1.7e7, [description = "Combustion coefficient", unit = u"J/kg"]
        f_evap = 0.03, [description = "Fraction of sensible heat for moisture evaporation (dimensionless)", unit = u"1"]
        f_water = 0.56, [description = "Fraction of dry fuel mass converted to water vapor (dimensionless)", unit = u"1"]

        # Latent heat of vaporization of water
        L_v = 2.5e6, [description = "Latent heat of vaporization of water", unit = u"J/kg"]

        # Reference values for non-dimensionalization
        one_vel = 1.0, [description = "Reference velocity for non-dimensionalization", unit = u"m/s"]
        four_vel = 4.0, [description = "Reference velocity for O2 limitation", unit = u"m/s"]
        one_dimless = 1.0, [description = "Dimensionless one", unit = u"1"]
        zero_dimless = 0.0, [description = "Dimensionless zero", unit = u"1"]
        zero_fuel = 0.0, [description = "Zero fuel mass", unit = u"kg/m^2"]

        # Heat flux extinction depth — Eq. 10, Clark et al. (1996)
        alpha_ext = 50.0, [description = "Heat flux e-folding extinction depth", unit = u"m"]

        # Canopy ignition threshold — p. 180, Clark et al. (1996)
        Q_canopy_threshold = 170000.0, [description = "Cumulative heat flux threshold for canopy ignition", unit = u"J/m^2"]
    end

    @parameters begin
        V_A, [description = "Tracer-advecting wind speed at 15 m AGL", unit = u"m/s"]
    end

    @variables begin
        M_litter(t), [description = "Litter fuel mass remaining", unit = u"kg/m^2"]
        M_trash(t), [description = "Trash fuel mass remaining", unit = u"kg/m^2"]
        M_scrub(t), [description = "Scrub fuel mass remaining", unit = u"kg/m^2"]
        M_canopy(t), [description = "Dry canopy fuel mass remaining", unit = u"kg/m^2"]
        Q_cumulative(t), [description = "Cumulative ground heat flux for canopy ignition", unit = u"J/m^2"]
        canopy_burning(t), [description = "Whether canopy is burning (1=yes, 0=no) (dimensionless)", unit = u"1"]
        B_ratio(t), [description = "Burn rate modulation factor (dimensionless)", unit = u"1"]
        S_f(t), [description = "Forward fire spread rate (McArthur)", unit = u"m/s"]
        ground_burn_rate(t), [description = "Total ground fuel burn rate", unit = u"kg/(m^2*s)"]
        total_burn_rate(t), [description = "Total fuel burn rate (ground + canopy)", unit = u"kg/(m^2*s)"]
        F_s(t), [description = "Sensible heat flux at surface from fire", unit = u"W/m^2"]
        F_l(t), [description = "Latent heat flux at surface from fire", unit = u"W/m^2"]
    end

    eqs = [
        # Burn rate modulation factor — Eq. 8, Clark et al. (1996)
        # B_ratio = sqrt((|V_A| + 1) / (|V_A| + 4))
        # Note: random factor r_f omitted (r_f ≈ 1 ± 0.05) for deterministic implementation
        B_ratio ~ sqrt((V_A + one_vel) / (V_A + four_vel)),

        # Forward fire spread rate — Eq. 9, Clark et al. (1996)
        # S_f = S_a * exp(k_spread * |V_A|)
        S_f ~ S_a * exp(k_spread * V_A),

        # Ground fuel consumption ODEs — derived from burn rates on p. 180
        # Each fuel type consumed at nominal rate * B_ratio while fuel remains
        D(M_litter) ~ -R_litter * B_ratio * ifelse(M_litter > zero_fuel, one_dimless, zero_dimless),
        D(M_trash) ~ -R_trash * B_ratio * ifelse(M_trash > zero_fuel, one_dimless, zero_dimless),
        D(M_scrub) ~ -R_scrub * B_ratio * ifelse(M_scrub > zero_fuel, one_dimless, zero_dimless),

        # Cumulative heat flux from ground fuels — p. 180, Clark et al. (1996)
        # Track cumulative heat energy for canopy ignition trigger
        D(Q_cumulative) ~ (one_dimless - f_evap) * H_c * ground_burn_rate,

        # Canopy ignition logic — p. 180, Clark et al. (1996)
        # Canopy ignites when cumulative ground heat flux exceeds 170 kJ/m²
        canopy_burning ~ ifelse(Q_cumulative > Q_canopy_threshold, one_dimless, zero_dimless),

        # Canopy fuel consumption ODE
        D(M_canopy) ~ -R_canopy * B_ratio * canopy_burning * ifelse(M_canopy > zero_fuel, one_dimless, zero_dimless),

        # Total burn rates
        ground_burn_rate ~ R_litter * B_ratio * ifelse(M_litter > zero_fuel, one_dimless, zero_dimless) +
            R_trash * B_ratio * ifelse(M_trash > zero_fuel, one_dimless, zero_dimless) +
            R_scrub * B_ratio * ifelse(M_scrub > zero_fuel, one_dimless, zero_dimless),
        total_burn_rate ~ ground_burn_rate +
            R_canopy * B_ratio * canopy_burning * ifelse(M_canopy > zero_fuel, one_dimless, zero_dimless),

        # Sensible heat flux at surface — derived from p. 180
        # 97% of combustion heat goes to sensible heat (3% to moisture evaporation)
        F_s ~ (one_dimless - f_evap) * H_c * total_burn_rate,

        # Latent heat flux at surface — derived from p. 180
        # Water vapor from combustion of cellulose (56% of dry fuel mass) plus evaporated moisture
        F_l ~ (f_water * total_burn_rate + f_evap * H_c * total_burn_rate / L_v) * L_v,
    ]

    return System(eqs, t; name)
end


"""
    Clark1996HeatFluxProfile(; name=:Clark1996HeatFluxProfile)

Compute the vertical profile of heat fluxes from the fire using a simple
extinction depth formulation (Eq. 10, Clark et al., 1996).

The sensible and latent heat fluxes decay exponentially with height:

```math
F_s(\\vec{x}, t) = F_s(x, y, 0, t) \\cdot \\exp(-z/\\alpha)
```

```math
F_l(\\vec{x}, t) = F_l(x, y, 0, t) \\cdot \\exp(-z/\\alpha)
```

where ``\\alpha = 50`` m is the e-folding extinction depth.

# Reference

Clark, T.L., Jenkins, M.A., Coen, J.L., and Packham, D.R. (1996). A Coupled
Atmosphere-Fire Model: Role of the Convective Froude Number and Dynamic Fingering at the
Fireline. *Int. J. Wildland Fire*, 6(4), 177-190.
"""
@component function Clark1996HeatFluxProfile(; name = :Clark1996HeatFluxProfile)
    @constants begin
        alpha_ext = 50.0, [description = "Heat flux e-folding extinction depth", unit = u"m"]
        one_m = 1.0, [description = "Reference length", unit = u"m"]
    end

    @parameters begin
        z, [description = "Height above ground level", unit = u"m"]
        F_s_sfc, [description = "Sensible heat flux at surface", unit = u"W/m^2"]
        F_l_sfc, [description = "Latent heat flux at surface", unit = u"W/m^2"]
    end

    @variables begin
        F_s(t), [description = "Sensible heat flux at height z", unit = u"W/m^2"]
        F_l(t), [description = "Latent heat flux at height z", unit = u"W/m^2"]
    end

    eqs = [
        # Eq. 10 — Clark et al. (1996)
        F_s ~ F_s_sfc * exp(-z / alpha_ext),
        F_l ~ F_l_sfc * exp(-z / alpha_ext),
    ]

    return System(eqs, t; name)
end


"""
    Clark1996ConvectiveFroudeNumber(; name=:Clark1996ConvectiveFroudeNumber)

Compute the square of the convective Froude number (Eq. 1, Clark et al., 1996),
a diagnostic parameter that characterizes the coupling strength between fire and
atmospheric dynamics.

```math
F_c^2 = \\frac{(U - S_f)^2}{g \\cdot \\frac{\\Delta\\theta}{\\bar{\\theta}} \\cdot W_f}
```

Large ``F_c^2`` indicates weak coupling (wind dominates); small ``F_c^2`` indicates
strong coupling (fire-induced buoyancy dominates).

# Reference

Clark, T.L., Jenkins, M.A., Coen, J.L., and Packham, D.R. (1996). A Coupled
Atmosphere-Fire Model: Role of the Convective Froude Number and Dynamic Fingering at the
Fireline. *Int. J. Wildland Fire*, 6(4), 177-190.
"""
@component function Clark1996ConvectiveFroudeNumber(; name = :Clark1996ConvectiveFroudeNumber)
    @constants begin
        g = 9.81, [description = "Gravitational acceleration", unit = u"m/s^2"]
        one_dimless = 1.0, [description = "Dimensionless one", unit = u"1"]
    end

    @parameters begin
        U, [description = "Ambient wind speed", unit = u"m/s"]
        S_f, [description = "Fire spread rate", unit = u"m/s"]
        delta_theta_over_theta, [description = "Convective buoyancy Δθ/θ̄ (dimensionless)", unit = u"1"]
        W_f, [description = "Fire width", unit = u"m"]
    end

    @variables begin
        F_c_sq(t), [description = "Square of convective Froude number (dimensionless)", unit = u"1"]
    end

    eqs = [
        # Eq. 1 — Clark et al. (1996)
        F_c_sq ~ (U - S_f)^2 / (g * delta_theta_over_theta * W_f),
    ]

    return System(eqs, t; name)
end


"""
    Clark1996WindProfile(; name=:Clark1996WindProfile)

Compute the hyperbolic tangent wind profile used in the FR7CS1 experiment
(Eq. 11, Clark et al., 1996).

```math
U(z) = 3 - 3 \\left[1 + \\tanh\\left(\\frac{z - 500}{100}\\right)\\right]
```

This profile gives ``U \\approx +3`` m/s near the surface, reverses sign at ``z = 500`` m,
and is asymptotic to ``-3`` m/s aloft. It allows fire-induced convective motions to
propagate upstream and re-enter the boundary layer flow that ventilates the fire.

# Reference

Clark, T.L., Jenkins, M.A., Coen, J.L., and Packham, D.R. (1996). A Coupled
Atmosphere-Fire Model: Role of the Convective Froude Number and Dynamic Fingering at the
Fireline. *Int. J. Wildland Fire*, 6(4), 177-190.
"""
@component function Clark1996WindProfile(; name = :Clark1996WindProfile)
    @constants begin
        U_base = 3.0, [description = "Base wind speed", unit = u"m/s"]
        z_reversal = 500.0, [description = "Height of wind reversal", unit = u"m"]
        z_scale = 100.0, [description = "Vertical scale of wind transition", unit = u"m"]
        one_dimless = 1.0, [description = "Dimensionless one", unit = u"1"]
    end

    @parameters begin
        z, [description = "Height above ground level", unit = u"m"]
    end

    @variables begin
        U_z(t), [description = "Wind speed at height z", unit = u"m/s"]
    end

    eqs = [
        # Eq. 11 — Clark et al. (1996)
        U_z ~ U_base - U_base * (one_dimless + tanh((z - z_reversal) / z_scale)),
    ]

    return System(eqs, t; name)
end
