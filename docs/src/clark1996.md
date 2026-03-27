# Clark et al. (1996) Coupled Atmosphere-Fire Model

## Overview

This module implements the fire model component of the coupled atmosphere-fire model
described by Clark et al. (1996). The model simulates fuel consumption and heat flux
generation for a simple dry eucalyptus forest fire, using the McArthur empirical fire
spread formula (Noble et al., 1980).

The model is designed to study the coupling between fire dynamics and atmospheric
circulations, characterized by the square of the convective Froude number, ``F_c^2``.
Small ``F_c^2`` indicates strong coupling between fire and atmosphere; large ``F_c^2``
indicates weak coupling where wind dominates.

The fire model tracks four fuel types — litter, trash, scrub, and canopy — with
corresponding burn rates and heat flux outputs. It includes an oxygen limitation
factor (``B_{ratio}``) that modulates burn rates at low wind speeds, and computes both
sensible and latent heat fluxes from combustion.

**Reference**: Clark, T.L., Jenkins, M.A., Coen, J.L., and Packham, D.R. (1996). A Coupled
Atmosphere-Fire Model: Role of the Convective Froude Number and Dynamic Fingering at the
Fireline. *Int. J. Wildland Fire*, 6(4), 177-190.

```@docs
Clark1996FireSpread
Clark1996HeatFluxProfile
Clark1996ConvectiveFroudeNumber
Clark1996WindProfile
```

## Implementation

The fire model consists of four components:

1. **`Clark1996FireSpread`**: Main fire model with fuel consumption ODEs and heat flux computation
2. **`Clark1996HeatFluxProfile`**: Vertical extinction profile for heat fluxes (Eq. 10)
3. **`Clark1996ConvectiveFroudeNumber`**: Diagnostic coupling parameter (Eq. 1)
4. **`Clark1996WindProfile`**: Hyperbolic tangent wind profile for FR7CS1 experiment (Eq. 11)

### State Variables

```@example clark1996_impl
using DataFrames, ModelingToolkit, Symbolics, DynamicQuantities
using WildlandFire

sys = Clark1996FireSpread()
vars = unknowns(sys)
units = [ModelingToolkit.get_unit(v) for v in vars]
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [isnothing(u) ? "dimensionless" : string(u) for u in units],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

### Parameters

```@example clark1996_impl
params = parameters(sys)
punits = [ModelingToolkit.get_unit(p) for p in params]
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params],
    :Units => [isnothing(u) ? "dimensionless" : string(u) for u in punits],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

### Equations

```@example clark1996_impl
equations(sys)
```

## Analysis

### Fire Spread Rate vs Wind Speed (Figure 3)

The McArthur fire spread rate (Eq. 9) and convective Froude number squared (Eq. 1)
are plotted as functions of the ambient wind speed ``U_0``. The spread rate increases
exponentially with wind speed according to the McArthur formula. The convective Froude
number squared increases more rapidly, indicating weakening fire-atmosphere coupling at
higher wind speeds.

The markers correspond to the experimental values from Table 1 of Clark et al. (1996).

```@example clark1996_impl
using Plots

# McArthur spread rate — Eq. 9
S_a = 0.18  # m/s
k_spread = 0.08424  # s/m
U_0_range = range(0, 30, length=200)
S_f_calc = [S_a * exp(k_spread * U) for U in U_0_range]

# Convective Froude number squared — Eq. 1
# Back-calculate g*Δθ/θ̄*W_f from Table 1: at U_0=3, S_f=0.23, F_c²=0.25
g_delta_W = (3.0 - 0.23)^2 / 0.25
F_c_sq_calc = [(U - S_a * exp(k_spread * U))^2 / g_delta_W for U in U_0_range]

# Table 1 data points
U_0_table = [1.0, 2.0, 3.0, 4.0, 5.0, 10.0, 15.0, 20.0]
S_f_table = [0.20, 0.21, 0.23, 0.25, 0.27, 0.42, 0.64, 0.97]
F_c_sq_table = [0.025, 0.11, 0.25, 0.43, 0.62, 1.7, 2.5, 2.9]
S_fa_table = [0.40, 0.39, 0.43, 0.46, 0.47, 0.56, 0.72, 1.04]
F_ca_sq_table = [0.0069, 0.051, 0.12, 0.21, 0.33, 1.2, 2.2, 2.6]

p = plot(collect(U_0_range), F_c_sq_calc, label="F²c", lw=2, color=:blue,
    xlabel="U₀ (m s⁻¹)", ylabel="F²c, Sf (m s⁻¹)",
    title="Fire Spread Rate and Convective Froude Number vs Wind Speed")
plot!(p, collect(U_0_range), S_f_calc, label="Sf", lw=2, color=:red, ls=:dash)
scatter!(p, U_0_table, F_c_sq_table, label="F²c (Table 1)", color=:blue, marker=:circle)
scatter!(p, U_0_table, S_f_table, label="Sf (Table 1)", color=:red, marker=:diamond)
scatter!(p, U_0_table, F_ca_sq_table, label="F²c_a (actual)", color=:cyan, marker=:utriangle)
scatter!(p, U_0_table, S_fa_table, label="Sf_a (actual)", color=:orange, marker=:utriangle)
p
```

### Fuel Consumption Time Series

Simulation of fuel consumption for the standard experiment conditions (FIR7CR: ``U_0 = 3`` m/s,
single 420 m fire line). Ground fuels (litter, trash, scrub) are consumed at different rates,
with litter being the fastest. Canopy fuel consumption is shown both with and without
canopy ignition.

```@example clark1996_impl
using OrdinaryDiffEqDefault

sys = Clark1996FireSpread()
compiled = mtkcompile(sys)

# Ground-only burning
prob_ground = ODEProblem(compiled,
    [compiled.M_litter => 2.0,
     compiled.M_trash => 0.5,
     compiled.M_scrub => 0.2,
     compiled.M_canopy => 1.2],
    (0.0, 200.0),
    [compiled.V_A => 3.0,
     compiled.canopy_burning => 0.0])
sol_ground = solve(prob_ground)

# Ground + canopy burning
prob_full = ODEProblem(compiled,
    [compiled.M_litter => 2.0,
     compiled.M_trash => 0.5,
     compiled.M_scrub => 0.2,
     compiled.M_canopy => 1.2],
    (0.0, 200.0),
    [compiled.V_A => 3.0,
     compiled.canopy_burning => 1.0])
sol_full = solve(prob_full)

p1 = plot(sol_ground.t, sol_ground[compiled.M_litter], label="Litter", lw=2, color=:brown,
    xlabel="Time (s)", ylabel="Fuel Mass (kg m⁻²)",
    title="Ground Fuel Only (U₀ = 3 m/s)")
plot!(p1, sol_ground.t, sol_ground[compiled.M_trash], label="Trash", lw=2, color=:orange)
plot!(p1, sol_ground.t, sol_ground[compiled.M_scrub], label="Scrub", lw=2, color=:green)
plot!(p1, sol_ground.t, sol_ground[compiled.M_canopy], label="Canopy (not burning)", lw=2, color=:darkgreen, ls=:dash)

p2 = plot(sol_full.t, sol_full[compiled.M_litter], label="Litter", lw=2, color=:brown,
    xlabel="Time (s)", ylabel="Fuel Mass (kg m⁻²)",
    title="Ground + Canopy Burning (U₀ = 3 m/s)")
plot!(p2, sol_full.t, sol_full[compiled.M_trash], label="Trash", lw=2, color=:orange)
plot!(p2, sol_full.t, sol_full[compiled.M_scrub], label="Scrub", lw=2, color=:green)
plot!(p2, sol_full.t, sol_full[compiled.M_canopy], label="Canopy", lw=2, color=:darkgreen)

plot(p1, p2, layout=(1,2), size=(900, 400))
```

### Sensible Heat Flux vs Wind Speed

The sensible heat flux at the surface varies with wind speed through the burn rate
modulation factor ``B_{ratio}`` (Eq. 8). Higher wind speeds increase ``B_{ratio}``
toward 1.0, increasing the effective burn rate and heat release.

```@example clark1996_impl
V_A_range = range(0, 20, length=100)
H_c = 1.7e7  # J/kg
f_evap = 0.03

R_ground = 0.040 + 0.005 + 0.004
R_canopy_rate = 0.020

B_ratio_vals = [sqrt((v + 1.0) / (v + 4.0)) for v in V_A_range]
F_s_ground = [(1.0 - f_evap) * H_c * R_ground * b / 1000 for b in B_ratio_vals]
F_s_total = [(1.0 - f_evap) * H_c * (R_ground + R_canopy_rate) * b / 1000 for b in B_ratio_vals]

p = plot(collect(V_A_range), F_s_ground, label="Ground fuel only", lw=2, color=:orange,
    xlabel="Wind Speed V_A (m s⁻¹)", ylabel="Sensible Heat Flux (kW m⁻²)",
    title="Initial Sensible Heat Flux vs Wind Speed")
plot!(p, collect(V_A_range), F_s_total, label="Ground + canopy", lw=2, color=:red)
p
```

### Burn Rate Modulation Factor

The burn rate ratio ``B_{ratio}`` (Eq. 8) models the oxygen limitation effect at low wind
speeds. It approaches 0.5 at zero wind and asymptotically approaches 1.0 at high wind speeds.

```@example clark1996_impl
V_A_range = range(0, 30, length=200)
B_vals = [sqrt((v + 1.0) / (v + 4.0)) for v in V_A_range]

p = plot(collect(V_A_range), B_vals, lw=2, color=:black, legend=false,
    xlabel="Wind Speed V_A (m s⁻¹)", ylabel="B_ratio (dimensionless)",
    title="Burn Rate Modulation Factor (Eq. 8)")
hline!(p, [0.5], ls=:dash, color=:gray, label="V_A = 0 limit")
hline!(p, [1.0], ls=:dot, color=:gray, label="V_A → ∞ limit")
p
```

### Heat Flux Vertical Profile

The heat flux decays exponentially with height above the fire, with an e-folding depth
of ``\alpha = 50`` m (Eq. 10). This simple extinction depth formulation parameterizes
the effects of radiation and small-scale turbulence.

```@example clark1996_impl
z_range = range(0, 200, length=100)
alpha = 50.0
decay = [exp(-z / alpha) for z in z_range]

p = plot(collect(z_range), decay, lw=2, color=:red, legend=false,
    xlabel="Height z (m)", ylabel="F(z) / F(0)",
    title="Heat Flux Extinction Profile (Eq. 10, α = 50 m)")
vline!(p, [50.0], ls=:dash, color=:gray)
p
```

### Wind Profile (Eq. 11, FR7CS1 Experiment)

The hyperbolic tangent wind profile used in the FR7CS1 experiment, which reverses
direction with height. This allows fire-induced convective motions to propagate
upstream and re-enter the fire ventilation zone.

```@example clark1996_impl
z_range = range(0, 1500, length=200)
U_z_vals = [3.0 - 3.0 * (1.0 + tanh((z - 500.0) / 100.0)) for z in z_range]

p = plot(U_z_vals, collect(z_range), lw=2, color=:blue,
    xlabel="U(z) (m s⁻¹)", ylabel="Height z (m)",
    title="Wind Profile (Eq. 11, FR7CS1)", legend=false)
vline!(p, [0.0], ls=:dash, color=:gray)
hline!(p, [500.0], ls=:dot, color=:gray)
p
```
