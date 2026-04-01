# Level-Set Fire Spread Model

## Overview

The level-set method for fire front propagation tracks the fire perimeter implicitly
as the zero contour of a level-set function ``\psi(x, y, t)``. The burning region is
defined by ``\psi \leq 0`` and unburned fuel by ``\psi > 0``. The fire front evolves
according to the anisotropic Hamilton-Jacobi equation (based on Eq. 9, Mandel et al. 2011):

```math
\frac{\partial \psi}{\partial t} + S(\hat{n}) \|\nabla \psi\| = 0
```

where the direction-dependent spread rate is given by the elliptical fire shape model
(Andrews 2018, Table 26):

```math
S(\gamma) = R_H \frac{1 - e}{1 - e \cos(\gamma)}
```

Here ``\hat{n} = \nabla\psi / \|\nabla\psi\|`` is the outward fire front normal,
``\gamma`` is the angle between ``\hat{n}`` and the head fire direction ``\alpha``,
``R_H`` is the head fire rate of spread, and ``e = \sqrt{Z^2 - 1}/Z`` is the
eccentricity derived from the fire length-to-width ratio ``Z``. When ``Z = 1``
(no wind or slope asymmetry), ``e = 0`` and the equation reduces to the isotropic
case ``S = R_H``.

This approach enables natural handling of merging fire fronts, complex terrain,
and realistic wind-driven elliptical fire shapes.

The implementation also includes fuel consumption and fire heat flux models that
couple the fire front to atmospheric and fuel state.

## Implementation Notes and Limitations

This implementation provides the core level-set fire spread algorithm using ModelingToolkit.jl
and MethodOfLines.jl for spatial discretization. While based on Muñoz-Esparza et al. (2018),
this version differs from the full WRF-Fire implementation in one key aspect:

**Current Implementation:**
- Fifth-order WENO (Weighted Essentially Non-Oscillatory) spatial discretization via MethodOfLines.jl
- Adaptive temporal integration via OrdinaryDiffEq.jl
- Anisotropic (elliptical) fire spread coupled through `FireSpreadDirection`

**Future work (full Muñoz-Esparza et al. 2018 algorithm):**
- Level-set reinitialization equation to maintain signed distance property

This implementation provides accurate solutions for both circular and elliptical fire
spread using the WENO scheme recommended by Muñoz-Esparza et al. (2018).

## Accuracy Considerations

Muñoz-Esparza et al. (2018) demonstrate that common level-set implementations with
standard finite differences yield rate-of-spread errors in the range 10–35% for typical
grid sizes (Δ = 12.5–100 m) and considerably underestimate fire area. The fifth-order
WENO scheme used in this implementation significantly reduces these errors.

**References**:

Muñoz-Esparza, D., Kosović, B., Jiménez, P.A., and Coen, J.L. (2018). An accurate
fire-spread algorithm in the Weather Research and Forecasting model using the level-set
method. *J. Adv. Model. Earth Syst.*, 10, 908–926. doi:10.1002/2017MS001108

Mandel, J., Beezley, J.D., and Kochanski, A.K. (2011). Coupled atmosphere-wildland fire
modeling with WRF 3.3 and SFIRE 2011. *Geosci. Model Dev.*, 4, 591–610.
doi:10.5194/gmd-4-591-2011

```@docs
LevelSetFireSpread
FuelConsumption
FireHeatFlux
anderson_fuel_coefficients
```

## Implementation

### Level-Set PDE System

The `LevelSetFireSpread` function returns a `PDESystem` representing the anisotropic
level-set equation on a 2D spatial domain with time. This implementation uses:

- **Spatial Discretization**: MethodOfLines.jl with fifth-order WENO scheme
- **Temporal Integration**: OrdinaryDiffEq.jl with adaptive time stepping
- **Boundary Conditions**: Neumann (zero gradient) by default, following Mandel et al. (2011) Sect. 3.4
- **Anisotropic Spread**: Elliptical fire shape model from Andrews (2018) Table 26

The anisotropic Hamilton-Jacobi equation is discretized as:
```math
\frac{\partial \psi}{\partial t} = -S(\gamma) \sqrt{\left(\frac{\partial \psi}{\partial x}\right)^2 + \left(\frac{\partial \psi}{\partial y}\right)^2}
```

where the direction-dependent speed ``S(\gamma) = R_H (1-e)/(1-e\cos\gamma)`` uses the
angle ``\gamma = \text{atan2}(\psi_y, \psi_x) - \alpha`` between the gradient and head fire direction.

```@example levelset
using DataFrames, ModelingToolkit, Symbolics, DynamicQuantities
using ModelingToolkit: t
using WildlandFire, EarthSciMLBase, DomainSets

@parameters x [unit = u"m"]
@parameters y [unit = u"m"]
domain = DomainInfo(
    constIC(0.0, t ∈ Interval(0.0, 10.0)),
    constBC(0.0, x ∈ Interval(0.0, 100.0), y ∈ Interval(0.0, 100.0)),
)
sys = LevelSetFireSpread(domain;
    initial_condition = (x, y) -> sqrt((x - 50.0)^2 + (y - 50.0)^2) - 10.0,
)

# Independent variables
DataFrame(
    :Name => [string(v) for v in sys.ivs],
)
```

```@example levelset
# Dependent variables
DataFrame(
    :Name => [string(v) for v in sys.dvs],
)
```

```@example levelset
# Parameters
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in sys.ps],
    :Description => [ModelingToolkit.getdescription(p) for p in sys.ps]
)
```

### Equations

```@example levelset
equations(sys)
```

### Boundary Conditions

```@example levelset
sys.bcs
```

### Fuel Consumption

The `FuelConsumption` component models exponential fuel decay after ignition
(Eq. 3, Mandel et al. 2011):

```math
F(t) = \exp\left(-\frac{t - t_i}{T_f}\right), \quad t > t_i
```

where ``t_i`` is the ignition time and ``T_f = w / 0.8514`` is the fuel burn time
constant derived from the fuel weight parameter ``w``. The burn time constant is
provided by `FuelModelLookup` via coupling, so it varies spatially by fuel type.

The effective fuel load ``w_{0,\text{eff}} = F \cdot w_{0,\text{initial}}`` is used
by `FireHeatFlux` to compute sensible and latent heat fluxes released to the
atmosphere. Per Mandel et al. (2011) Section 3.2, the fire spread rate depends on
the original fuel model properties, not the remaining fuel fraction.

```@example levelset
fc = FuelConsumption()
vars = unknowns(fc)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [string(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

```@example levelset
equations(fc)
```

### Fire Heat Flux

The `FireHeatFlux` component computes sensible and latent heat fluxes from burning
fuel (Eqs. 4–5, Mandel et al. 2011).

Sensible heat flux (Eq. 4):
```math
\phi_h = \frac{-dF/dt \cdot w_\ell \cdot h}{1 + M_f} \quad (W\,m^{-2})
```

Latent heat flux (Eq. 5):
```math
\phi_q = \frac{-dF/dt \cdot (M_f + 0.56) \cdot L \cdot w_\ell}{1 + M_f} \quad (W\,m^{-2})
```

where 0.56 is the estimated mass ratio of water output from combustion to dry fuel,
and ``L = 2.5 \times 10^6`` J/kg is the latent heat of condensation of water at 0°C.

```@example levelset
fhf = FireHeatFlux()
vars = unknowns(fhf)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [string(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

```@example levelset
equations(fhf)
```

### Fire Spread Rate (Eq. 2)

The fire spread rate is computed using the modified Rothermel formula
(Eq. 2, Mandel et al. 2011):

```math
S = \max\{S_0,\; R_0 + c \cdot \min\{e,\; \max\{0, U\}\}^b + d \cdot \max\{0,\; \tan\phi\}^2\}
```

where ``U`` is the wind speed normal to the fire front (m/s), ``\tan\phi`` is the
terrain slope in the normal direction, and ``S_0, R_0, b, c, d, e`` are
fuel-dependent coefficients computed from the Rothermel (1972) model following
Table 2 of Mandel et al. (2011).

### Anderson Fuel Model Coefficients

The `anderson_fuel_coefficients` function precomputes these Rothermel spread rate
coefficients for the 13 Anderson (1982) fuel models:

```@example levelset
coeffs = anderson_fuel_coefficients(1)
DataFrame(
    :Field => [string(k) for k in keys(coeffs)],
    :Value => [coeffs[k] for k in keys(coeffs)],
    :Units => ["m/s", "m/s", "—", "(m/s)^(1-b)", "m/s", "m/s", "s", "kg/m²", "J/kg", "m"],
)
```

### Fuel Properties (Table 1, Mandel et al. 2011)

The fuel is characterized by the quantities listed in Table 1 of
Mandel et al. (2011). The 13 Anderson (1982) fuel categories provide
preset vectors of these properties. Note that some quantities are stored
in English units per the Rothermel (1972) convention (surface-area-to-volume
ratio ``\sigma`` in 1/ft, fuel particle density ``\rho_p`` in lb/ft³, and
heat content ``h`` in BTU/lb).

```@example levelset
fuel_names = [
    "Short grass", "Timber grass/understory", "Tall grass",
    "Chaparral", "Brush", "Dormant brush",
    "Southern rough", "Compact timber litter", "Hardwood litter",
    "Timber (understory)", "Light logging slash", "Medium logging slash",
    "Heavy logging slash"
]
DataFrame(
    :Model => 1:13,
    :Name => fuel_names,
    Symbol("w_ℓ (kg/m²)") => [ANDERSON_FUEL_DATA[i].fgi for i in 1:13],
    Symbol("δ (m)") => [ANDERSON_FUEL_DATA[i].depth for i in 1:13],
    Symbol("σ (1/ft)") => [ANDERSON_FUEL_DATA[i].savr for i in 1:13],
    Symbol("M_x") => [ANDERSON_FUEL_DATA[i].mce for i in 1:13],
    Symbol("ρ_p (lb/ft³)") => [ANDERSON_FUEL_DATA[i].dens for i in 1:13],
    Symbol("S_T") => [ANDERSON_FUEL_DATA[i].st for i in 1:13],
    Symbol("S_E") => [ANDERSON_FUEL_DATA[i].se for i in 1:13],
    Symbol("h (BTU/lb)") => [ANDERSON_FUEL_DATA[i].h for i in 1:13],
    Symbol("w (s)") => [ANDERSON_FUEL_DATA[i].weight for i in 1:13],
)
```

### Derived Quantities

The computed spread rate coefficients and burn parameters for each fuel model:

```@example levelset
all_coeffs = [anderson_fuel_coefficients(i) for i in 1:13]
DataFrame(
    :Model => 1:13,
    :Name => fuel_names,
    Symbol("R₀ (m/s)") => [round(c.R_0, sigdigits=4) for c in all_coeffs],
    Symbol("T_f (s)") => [round(c.T_f, sigdigits=4) for c in all_coeffs],
    Symbol("w_l (kg/m²)") => [c.w_l for c in all_coeffs],
    Symbol("δ (m)") => [c.delta for c in all_coeffs],
)
```

### Rothermel Spread Rate Computation (Table 2, Mandel et al. 2011)

The computation of the fire spread rate factors from fuel properties
follows the equations in Table 2 of Mandel et al. (2011), originally
from Rothermel (1972). Here we show the full coefficients for fuel
model 1 (short grass) at 8% moisture:

```@example levelset
coeffs1 = anderson_fuel_coefficients(1; M_f = 0.08)
DataFrame(
    :Coefficient => ["S₀", "R₀", "b", "c", "d", "e"],
    :Description => [
        "Minimum spread rate",
        "No-wind no-slope spread rate",
        "Wind speed exponent",
        "Wind factor coefficient (used as c·U^b)",
        "Slope factor coefficient",
        "Maximum wind speed",
    ],
    :Value => [coeffs1.S_0, coeffs1.R_0, coeffs1.b, coeffs1.c, coeffs1.d, coeffs1.e],
    :Units => ["m/s", "m/s", "—", "(m/s)^(1-b)", "m/s", "m/s"],
)
```

## Analysis

### Circular Fire Spread (Isotropic)

With a constant spread rate ``S`` and no wind or slope, the level-set equation
produces circular fire spread. Starting from a circular ignition of radius
``r_0``, the fire radius at time ``t`` should be ``r_0 + S \cdot t``.
This is the fundamental analytical solution for the level-set equation
(Mandel et al. 2011, Sect. 3.4).

```@example levelset
using MethodOfLines, OrdinaryDiffEqDefault, Plots

r0 = 10.0       # initial radius (m)
S_val = 1.0      # spread rate (m/s)
domain_size = 100.0
center = domain_size / 2.0
t_end = 10.0

domain_circ = DomainInfo(
    constIC(0.0, t ∈ Interval(0.0, t_end)),
    constBC(0.0, x ∈ Interval(0.0, domain_size), y ∈ Interval(0.0, domain_size)),
)
sys = LevelSetFireSpread(domain_circ;
    initial_condition = (x, y) -> sqrt((x - center)^2 + (y - center)^2) - r0,
    spread_rate = S_val,
)

dx = 5.0
discretization = MOLFiniteDifference([sys.ivs[2] => dx, sys.ivs[3] => dx], sys.ivs[1];
    advection_scheme = WENOScheme())
prob = MethodOfLines.discretize(sys, discretization; checks = false)
sol = solve(prob; saveat = [0.0, t_end / 2, t_end])

psi = sol[sys.dvs[1]]
x_grid = sol[sys.ivs[2]]
y_grid = sol[sys.ivs[3]]

p1 = contour(x_grid, y_grid, psi[1, :, :]', levels = [0.0],
    title = "t = 0 s", xlabel = "x (m)", ylabel = "y (m)",
    aspect_ratio = :equal, linewidth = 2, label = "Fire front")
p2 = contour(x_grid, y_grid, psi[2, :, :]', levels = [0.0],
    title = "t = $(t_end/2) s", xlabel = "x (m)", ylabel = "y (m)",
    aspect_ratio = :equal, linewidth = 2, label = "Fire front")
p3 = contour(x_grid, y_grid, psi[end, :, :]', levels = [0.0],
    title = "t = $t_end s", xlabel = "x (m)", ylabel = "y (m)",
    aspect_ratio = :equal, linewidth = 2, label = "Fire front")

plot(p1, p2, p3, layout = (1, 3), size = (900, 300))
savefig("circular_spread.png"); nothing # hide
```

![Circular fire spread: the level-set method propagates the ψ=0 contour outward at rate S, producing approximately circular expansion from a circular initial condition (Eq. 9, Mandel et al. 2011).](circular_spread.png)

### Elliptical Fire Spread (Anisotropic)

When wind or slope creates asymmetric fire behavior (``Z > 1``), the fire shape
becomes elliptical. The direction-dependent spread rate ``S(\gamma) = R_H (1-e)/(1-e\cos\gamma)``
(Andrews 2018, Table 26) produces faster spread in the head fire direction (``\alpha``)
and slower spread in the backing direction. Here we demonstrate this with ``Z = 2``
and head fire direction along the +x axis (``\alpha = 0``):

```@example levelset
r0 = 10.0       # initial radius (m)
R_H_val = 1.0    # head fire rate of spread (m/s)
Z_val = 2.0      # length-to-width ratio
domain_size = 150.0
center = domain_size / 2.0
t_end = 10.0

domain_ellip = DomainInfo(
    constIC(0.0, t ∈ Interval(0.0, t_end)),
    constBC(0.0, x ∈ Interval(0.0, domain_size), y ∈ Interval(0.0, domain_size)),
)
sys_ellip = LevelSetFireSpread(domain_ellip;
    initial_condition = (x, y) -> sqrt((x - center)^2 + (y - center)^2) - r0,
    spread_rate = R_H_val,
)

# Set Z > 1 for elliptical spread, α = 0 for head fire along +x
for p in sys_ellip.ps
    sym = Symbolics.tosymbol(p, escape = false)
    if sym == :Z
        sys_ellip.initial_conditions[Symbolics.unwrap(p)] = Z_val
    elseif sym == :α
        sys_ellip.initial_conditions[Symbolics.unwrap(p)] = 0.0
    end
end

dx = 5.0
disc_ellip = MOLFiniteDifference(
    [sys_ellip.ivs[2] => dx, sys_ellip.ivs[3] => dx], sys_ellip.ivs[1];
    advection_scheme = WENOScheme())
prob_ellip = MethodOfLines.discretize(sys_ellip, disc_ellip; checks = false)
sol_ellip = solve(prob_ellip; saveat = [0.0, t_end / 2, t_end])

psi_e = sol_ellip[sys_ellip.dvs[1]]
x_grid_e = sol_ellip[sys_ellip.ivs[2]]
y_grid_e = sol_ellip[sys_ellip.ivs[3]]

p1 = contour(x_grid_e, y_grid_e, psi_e[1, :, :]', levels = [0.0],
    title = "t = 0 s", xlabel = "x (m)", ylabel = "y (m)",
    aspect_ratio = :equal, linewidth = 2, label = "Fire front")
p2 = contour(x_grid_e, y_grid_e, psi_e[2, :, :]', levels = [0.0],
    title = "t = $(t_end/2) s", xlabel = "x (m)", ylabel = "y (m)",
    aspect_ratio = :equal, linewidth = 2, label = "Fire front")
p3 = contour(x_grid_e, y_grid_e, psi_e[end, :, :]', levels = [0.0],
    title = "t = $t_end s", xlabel = "x (m)", ylabel = "y (m)",
    aspect_ratio = :equal, linewidth = 2, label = "Fire front")

plot(p1, p2, p3, layout = (1, 3), size = (900, 300),
    plot_title = "Elliptical Fire Spread (Z = $Z_val, α = 0)")
savefig("elliptical_spread.png"); nothing # hide
```

![Elliptical fire spread: with Z = 2 and head fire along +x (α = 0), the fire elongates in the downwind direction. The backing fire (−x) spreads more slowly, producing the characteristic elliptical shape (Andrews 2018, Table 26).](elliptical_spread.png)

### Fuel Consumption Dynamics (Eq. 3)

The fuel fraction ``F(t)`` decays exponentially with time constant ``T_f``
after ignition, following Eq. 3 of Mandel et al. (2011). The fuel weight
parameter ``w`` determines the burn time through ``T_f = w / 0.8514``.
For ``w = 1000`` s, fuel burns to 60% of its original quantity in 600 s
(Table 1, Mandel et al. 2011).

```@example levelset
fc = FuelConsumption()
compiled_fc = mtkcompile(fc)

T_f_val = 10.0
prob_fc = ODEProblem(
    compiled_fc,
    Dict(compiled_fc.F => 1.0),
    (0.0, 40.0),
    Dict(compiled_fc.T_f => T_f_val, compiled_fc.is_burning => 1.0, compiled_fc.w0_initial => 1.0),
)
sol_fc = solve(prob_fc)

t_plot = 0:0.1:40
F_numerical = [sol_fc(ti, idxs = compiled_fc.F) for ti in t_plot]
F_analytical = exp.(-t_plot ./ T_f_val)

plot(t_plot, F_numerical, label = "Numerical", linewidth = 2)
plot!(t_plot, F_analytical, label = "Analytical: exp(-t/T_f)", linestyle = :dash, linewidth = 2)
xlabel!("Time (s)")
ylabel!("Fuel fraction F")
title!("Fuel Consumption (Eq. 3, Mandel et al. 2011, T_f = $(T_f_val) s)")
savefig("fuel_consumption.png"); nothing # hide
```

![Fuel consumption: exponential decay of fuel fraction after ignition. At t = T_f, F ≈ 1/e ≈ 0.368 (Eq. 3, Mandel et al. 2011).](fuel_consumption.png)

### Anderson Fuel Model Comparison

Comparison of no-wind no-slope spread rates ``R_0`` across the 13 Anderson
fuel models at 8% fuel moisture, computed using the Rothermel coefficients
from Table 2 of Mandel et al. (2011):

```@example levelset
R0_vals = [anderson_fuel_coefficients(i).R_0 for i in 1:13]

bar(1:13, R0_vals .* 60,  # Convert m/s to m/min
    xticks = (1:13, string.(1:13)),
    xlabel = "Anderson Fuel Model",
    ylabel = "R₀ (m/min)",
    title = "No-Wind No-Slope Spread Rate (8% moisture)\n(Table 2, Mandel et al. 2011)",
    legend = false,
    bar_width = 0.7)
savefig("anderson_fuel_models.png"); nothing # hide
```

![Anderson fuel model comparison: no-wind no-slope spread rates for the 13 Anderson (1982) fuel categories, computed following the Rothermel model as described in Table 2 of Mandel et al. (2011).](anderson_fuel_models.png)

### Moisture Sensitivity (Table 2, Mandel et al. 2011)

Higher fuel moisture reduces spread rate by absorbing heat that would
otherwise preheat adjacent fuel. The moisture damping coefficient
``\eta_M`` (Eq. 29/30, Rothermel 1972) decreases the reaction intensity
as moisture approaches the moisture of extinction ``M_x``:

```@example levelset
moisture_range = 0.02:0.01:0.11
R0_fm1 = [anderson_fuel_coefficients(1; M_f = m).R_0 for m in moisture_range]
R0_fm3 = [anderson_fuel_coefficients(3; M_f = m).R_0 for m in moisture_range]

plot(moisture_range .* 100, R0_fm1 .* 60,
    label = "FM 1 (Short grass)", linewidth = 2, marker = :circle)
plot!(moisture_range .* 100, R0_fm3 .* 60,
    label = "FM 3 (Tall grass)", linewidth = 2, marker = :square)
xlabel!("Fuel Moisture (%)")
ylabel!("R₀ (m/min)")
title!("Spread Rate vs Fuel Moisture\n(Table 2, Mandel et al. 2011)")
savefig("moisture_sensitivity.png"); nothing # hide
```

![Moisture sensitivity: higher fuel moisture reduces the fire spread rate through the moisture damping coefficient η_M (Table 2, Mandel et al. 2011).](moisture_sensitivity.png)
