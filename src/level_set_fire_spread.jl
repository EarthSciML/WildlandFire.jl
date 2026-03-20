"""
    LevelSetFireSpread(; name=:LevelSetFireSpread, x_domain, y_domain,
                         t_domain, initial_condition, boundary_conditions=nothing,
                         spread_rate=1.0)

Create a level-set fire front propagation PDE system.

This implements the level-set method for tracking fire front evolution on a 2D domain,
as described in Muñoz-Esparza et al. (2018) and Mandel et al. (2011). The fire front
is represented implicitly as the zero contour of a level-set function ψ(x, y, t), where
ψ ≤ 0 denotes the burning region and ψ > 0 denotes unburned fuel.

The level-set function evolves according to the Hamilton-Jacobi equation (Eq. 9,
Mandel et al. 2011):

```math
\\frac{\\partial \\psi}{\\partial t} + S \\|\\nabla \\psi\\| = 0
```

where S is the fire spread rate (m/s). Spatial discretization is handled by
MethodOfLines.jl and time stepping by OrdinaryDiffEq.jl.

# Arguments
- `name`: System name (default `:LevelSetFireSpread`)
- `x_domain`: Tuple (x_min, x_max) defining the x-extent of the domain (m)
- `y_domain`: Tuple (y_min, y_max) defining the y-extent of the domain (m)
- `t_domain`: Tuple (t_min, t_max) defining the time interval (s)
- `initial_condition`: Function `(x, y) -> ψ₀` giving the initial level-set field.
  Negative values indicate initial fire region, positive values indicate unburned fuel.
- `boundary_conditions`: Optional vector of boundary condition equations. If not provided,
  Neumann (zero-gradient) boundary conditions are used.
- `spread_rate`: Default value for the fire spread rate parameter S (m/s). Default is 1.0.

# Returns
A `PDESystem` representing the level-set fire spread equation.

# References

Muñoz-Esparza, D., Kosović, B., Jiménez, P.A., and Coen, J.L. (2018). An accurate
fire-spread algorithm in the Weather Research and Forecasting model using the level-set
method. *J. Adv. Model. Earth Syst.*, 10, 908–926. doi:10.1002/2017MS001108

Mandel, J., Beezley, J.D., and Kochanski, A.K. (2011). Coupled atmosphere-wildland fire
modeling with WRF 3.3 and SFIRE 2011. *Geosci. Model Dev.*, 4, 591–610.
doi:10.5194/gmd-4-591-2011

# Example

```julia
using WildlandFire, MethodOfLines, DomainSets, OrdinaryDiffEqDefault

# Domain: 500m x 500m, 60 seconds
x_domain = (0.0, 500.0)
y_domain = (0.0, 500.0)
t_domain = (0.0, 60.0)

# Circular ignition at center (radius 10m)
initial_condition(x, y) = sqrt((x - 250.0)^2 + (y - 250.0)^2) - 10.0

# Create PDE system
sys = LevelSetFireSpread(;
    x_domain, y_domain, t_domain, initial_condition,
)

# Discretize and solve
dx = 5.0
discretization = MOLFiniteDifference([sys.ivs[2] => dx, sys.ivs[3] => dx], sys.ivs[1])
prob = MethodOfLines.discretize(sys, discretization; checks=false)
sol = solve(prob)
```
"""
function LevelSetFireSpread(;
        name = :LevelSetFireSpread,
        x_domain::Tuple{Float64, Float64},
        y_domain::Tuple{Float64, Float64},
        t_domain::Tuple{Float64, Float64},
        initial_condition,
        boundary_conditions = nothing,
        spread_rate = 1.0
    )

    @parameters x [description = "Horizontal coordinate x", unit = u"m"]
    @parameters y [description = "Horizontal coordinate y", unit = u"m"]

    # Fire spread rate with default value
    @parameters S = spread_rate [description = "Fire spread rate", unit = u"m/s"]

    # Level-set function
    @variables ψ(..) [description = "Level-set function (fire front at ψ=0)", unit = u"m"]

    # Spatial derivative operators
    Dx = Differential(x)
    Dy = Differential(y)

    # Gradient components — Eq. 6, Mandel et al. (2011)
    # ψ has units m, x and y have units m, so Dx(ψ) and Dy(ψ) are dimensionless
    ψ_x = Dx(ψ(t, x, y))
    ψ_y = Dy(ψ(t, x, y))

    # Level-set equation — Eq. 9, Mandel et al. (2011)
    # ∂ψ/∂t + S‖∇ψ‖ = 0
    # ‖∇ψ‖ = sqrt(ψ_x² + ψ_y²) is dimensionless
    # D(ψ) has units m/s, S has units m/s, so S * ‖∇ψ‖ has units m/s
    eq = [
        D(ψ(t, x, y)) ~ -S * sqrt(ψ_x^2 + ψ_y^2),  # Eq. 9
    ]

    # Domain definitions
    @constants begin
        ψ_ref = 1.0, [description = "Reference length for initial condition", unit = u"m"]
    end

    domains = [
        t ∈ Interval(t_domain[1], t_domain[2]),
        x ∈ Interval(x_domain[1], x_domain[2]),
        y ∈ Interval(y_domain[1], y_domain[2]),
    ]

    # Initial condition
    ic = ψ(t_domain[1], x, y) ~ initial_condition(x, y) * ψ_ref

    # Boundary conditions: default to Neumann (zero gradient) — Mandel et al. (2011) Sect. 3.4
    if boundary_conditions === nothing
        bcs = [
            ic,
            Dx(ψ(t, x_domain[1], y)) ~ 0.0,
            Dx(ψ(t, x_domain[2], y)) ~ 0.0,
            Dy(ψ(t, x, y_domain[1])) ~ 0.0,
            Dy(ψ(t, x, y_domain[2])) ~ 0.0,
        ]
    else
        bcs = [ic; boundary_conditions]
    end

    return PDESystem(
        eq, bcs, domains, [t, x, y], [ψ(t, x, y)],
        [S, ψ_ref]; name = name
    )
end


"""
    FuelConsumption(; name=:FuelConsumption)

Create a fuel consumption model system.

The fuel fraction remaining at a point decreases exponentially after ignition,
following Eq. 3 of Mandel et al. (2011):

```math
F(t) = \\exp\\left(-\\frac{t - t_i}{T_f}\\right), \\quad t > t_i
```

where ``t_i`` is the ignition time and ``T_f`` is the fuel burn time constant.

The fuel weight parameter ``w`` relates to ``T_f`` by:
``T_f \\approx w / 0.8514``

# Reference

Mandel, J., Beezley, J.D., and Kochanski, A.K. (2011). Coupled atmosphere-wildland fire
modeling with WRF 3.3 and SFIRE 2011. *Geosci. Model Dev.*, 4, 591–610.
doi:10.5194/gmd-4-591-2011
"""
@component function FuelConsumption(; name = :FuelConsumption)
    @parameters begin
        T_f, [description = "Fuel burn time constant (1/e decay time)", unit = u"s"]
    end

    @variables begin
        F(t), [description = "Fuel fraction remaining (1=full, 0=consumed)", unit = u"1"]
        is_burning(t), [description = "Whether the fuel is burning (1=yes, 0=no) (dimensionless)", unit = u"1"]
    end

    # Fuel fraction decreases exponentially — Eq. 3, Mandel et al. (2011)
    # dF/dt = -F/T_f when burning, 0 otherwise
    eqs = [
        D(F) ~ -is_burning * F / T_f,  # Eq. 3 (differential form)
    ]

    return System(eqs, t; name)
end


"""
    FireHeatFlux(; name=:FireHeatFlux)

Create a fire heat flux model system.

Computes the sensible and latent heat fluxes released by burning fuel.
The sensible heat flux follows Eq. 4 and the latent heat flux follows
Eq. 5 of Mandel et al. (2011).

Sensible heat flux:
```math
\\phi_h = \\frac{-dF/dt \\cdot w_\\ell \\cdot h}{1 + M_f}
```

Latent heat flux:
```math
\\phi_q = \\frac{-dF/dt \\cdot (M_f + 0.56) \\cdot L \\cdot w_\\ell}{1 + M_f}
```

where 0.56 is the estimated mass ratio of water output from combustion to dry fuel,
and L = 2.5 times 10^6 J/kg is the latent heat of condensation of water at 0 degrees C.

# Reference

Mandel, J., Beezley, J.D., and Kochanski, A.K. (2011). Coupled atmosphere-wildland fire
modeling with WRF 3.3 and SFIRE 2011. *Geosci. Model Dev.*, 4, 591–610.
doi:10.5194/gmd-4-591-2011
"""
@component function FireHeatFlux(; name = :FireHeatFlux)
    @constants begin
        one = 1.0, [description = "Dimensionless one", unit = u"1"]
        # Estimated mass ratio of water output from combustion to dry fuel
        water_combustion_ratio = 0.56, [description = "Water-to-fuel mass ratio from combustion (dimensionless)", unit = u"1"]
        # Latent heat of condensation of water at 0 degrees C
        L_water = 2.5e6, [description = "Latent heat of condensation of water", unit = u"J/kg"]
    end

    @parameters begin
        w_l, [description = "Total fuel load", unit = u"kg/m^2"]
        h_fuel, [description = "Heat content of dry fuel", unit = u"J/kg"]
        M_f, [description = "Fuel moisture content (dimensionless)", unit = u"1"]
    end

    @variables begin
        fuel_burn_rate(t), [description = "Rate of fuel consumption (-dF/dt)", unit = u"1/s"]
        phi_h(t), [description = "Sensible heat flux density", unit = u"W/m^2"]
        phi_q(t), [description = "Latent heat (moisture) flux density", unit = u"W/m^2"]
    end

    eqs = [
        # Sensible heat flux — Eq. 4, Mandel et al. (2011)
        phi_h ~ fuel_burn_rate * w_l * h_fuel / (one + M_f),

        # Latent heat flux — Eq. 5, Mandel et al. (2011)
        phi_q ~ fuel_burn_rate * (M_f + water_combustion_ratio) / (one + M_f) * L_water * w_l,
    ]

    return System(eqs, t; name)
end


# Anderson 13 fuel model data from WRF-SFIRE (Table 1, Mandel et al. 2011)
# fgi: total fuel load (kg/m²), depth: fuel bed depth (m),
# savr: surface-area-to-volume ratio (1/ft), mce: moisture of extinction (1),
# dens: ovendry fuel particle density (lb/ft³), st: total mineral content (1),
# se: effective mineral content (1), weight: fuel weight parameter (s),
# h: heat content (BTU/lb)
const ANDERSON_FUEL_DATA = Dict(
    1 => (fgi = 0.166, depth = 0.305, savr = 3500, mce = 0.12, dens = 32.0, st = 0.0555, se = 0.01, weight = 7, h = 8000.0),
    2 => (fgi = 0.896, depth = 0.305, savr = 3000, mce = 0.15, dens = 32.0, st = 0.0555, se = 0.01, weight = 7, h = 8000.0),
    3 => (fgi = 0.674, depth = 0.762, savr = 1500, mce = 0.25, dens = 32.0, st = 0.0555, se = 0.01, weight = 7, h = 8000.0),
    4 => (fgi = 3.591, depth = 1.829, savr = 1739, mce = 0.2, dens = 32.0, st = 0.0555, se = 0.01, weight = 180, h = 8000.0),
    5 => (fgi = 0.784, depth = 0.61, savr = 1683, mce = 0.2, dens = 32.0, st = 0.0555, se = 0.01, weight = 25, h = 8000.0),
    6 => (fgi = 1.344, depth = 0.762, savr = 1564, mce = 0.25, dens = 32.0, st = 0.0555, se = 0.01, weight = 25, h = 8000.0),
    7 => (fgi = 1.091, depth = 0.762, savr = 1562, mce = 0.4, dens = 32.0, st = 0.0555, se = 0.01, weight = 25, h = 8000.0),
    8 => (fgi = 1.12, depth = 0.061, savr = 1889, mce = 0.3, dens = 32.0, st = 0.0555, se = 0.01, weight = 900, h = 8000.0),
    9 => (fgi = 0.78, depth = 0.061, savr = 2484, mce = 0.25, dens = 32.0, st = 0.0555, se = 0.01, weight = 900, h = 8000.0),
    10 => (fgi = 2.694, depth = 0.305, savr = 1764, mce = 0.25, dens = 32.0, st = 0.0555, se = 0.01, weight = 900, h = 8000.0),
    11 => (fgi = 2.582, depth = 0.305, savr = 1182, mce = 0.15, dens = 32.0, st = 0.0555, se = 0.01, weight = 900, h = 8000.0),
    12 => (fgi = 7.749, depth = 0.701, savr = 1145, mce = 0.2, dens = 32.0, st = 0.0555, se = 0.01, weight = 900, h = 8000.0),
    13 => (fgi = 13.024, depth = 0.914, savr = 1159, mce = 0.25, dens = 32.0, st = 0.0555, se = 0.01, weight = 900, h = 8000.0),
)

"""
    anderson_fuel_coefficients(fuel_model::Int; M_f=0.08)

Return precomputed Rothermel spread rate coefficients for an Anderson (1982) fuel model,
along with fuel consumption parameters, following the form used in WRF-SFIRE.

The fire spread rate in WRF-SFIRE is computed as (Eq. 2, Mandel et al. 2011):

    S = max{S_0, R_0 + c * min{e, max{0, U}}^b + d * max{0, tan(phi)}^2}

where U is the normal wind component (m/s) and tan(phi) is the terrain slope in the
fire front normal direction.

# Arguments
- `fuel_model::Int`: Anderson fuel model number (1-13)
- `M_f`: Fuel moisture content (default 0.08, i.e. 8%)

# Returns
A `NamedTuple` with fields:
- `S_0`: Minimum spread rate (m/s)
- `R_0`: No-wind no-slope spread rate (m/s)
- `b`: Wind exponent (dimensionless)
- `c`: Wind factor coefficient, used as c·U^b in Eq. 2 (units: (m/s)^(1-b))
- `d`: Slope factor coefficient (m/s)
- `e`: Maximum wind speed for spread rate (m/s)
- `T_f`: Fuel burn time constant (s)
- `w_l`: Total fuel load (kg/m^2)
- `h`: Heat content (J/kg)
- `delta`: Fuel bed depth (m)

# Reference

Anderson, H. E. (1982). Aids to determining fuel models for estimating fire behavior.
USDA Forest Service, Intermountain Forest and Range Experiment Station. Gen. Tech. Rep. INT-122.

Mandel, J., Beezley, J.D., and Kochanski, A.K. (2011). Coupled atmosphere-wildland fire
modeling with WRF 3.3 and SFIRE 2011. *Geosci. Model Dev.*, 4, 591-610.
"""
function anderson_fuel_coefficients(fuel_model::Int; M_f = 0.08)
    if !haskey(ANDERSON_FUEL_DATA, fuel_model)
        error("Unknown fuel model: $fuel_model. Valid models are 1-13.")
    end

    fd = ANDERSON_FUEL_DATA[fuel_model]

    # Compute Rothermel spread rate coefficients in English units (BTU-lb-ft-min)
    # then convert to SI, following Table 2 of Mandel et al. (2011).
    # The Rothermel (1972) model was calibrated in English units, so we perform
    # the computation in English units and convert the final results to SI.
    #
    # Unit conversions used:
    #   1 ft  = 0.3048 m
    #   1 lb  = 0.453592 kg
    #   1 BTU = 1055.06 J
    #   1 lb/ft^2 = 4.88243 kg/m^2  (so 1 kg/m^2 = 0.204816 lb/ft^2)
    #   1 ft/min  = 0.00508 m/s     (= 0.3048/60)
    #   1 m/s     = 196.85 ft/min   (= 60/0.3048)

    # Convert inputs from SI to English units
    fgi_lbft2 = fd.fgi * 0.204816     # kg/m^2 → lb/ft^2
    depth_ft = fd.depth / 0.3048       # m → ft
    sigma = fd.savr                    # 1/ft (already in English units)
    rho_p = fd.dens                    # lb/ft^3 (already in English units)
    h_us = fd.h                          # BTU/lb (heat content from fuel model data)

    # Fuel bed properties — Table 2, Mandel et al. (2011)
    w_0 = fgi_lbft2 / (1.0 + M_f)     # Total fuel load net of moisture (Eq. T2-7)
    w_n = fgi_lbft2 / (1.0 + fd.st)   # Fuel loading net of minerals (Eq. T2-6)
    rho_b = w_0 / depth_ft            # Oven dry bulk density (Eq. T2-11)
    beta = rho_b / rho_p              # Packing ratio (Eq. T2-10)
    beta_op = 3.348 * sigma^(-0.8189) # Optimum packing ratio (Eq. T2-21)
    beta_ratio = beta / beta_op

    # Reaction intensity components — Table 2, Mandel et al. (2011)
    coeff_A = 1.0 / (4.774 * sigma^0.1 - 7.27)  # Eq. T2-12
    Gamma_max = sigma^1.5 / (495.0 + 0.0594 * sigma^1.5)  # Eq. T2-9
    Gamma = Gamma_max * beta_ratio^coeff_A * exp(coeff_A * (1.0 - beta_ratio))  # Eq. T2-8

    rM = min(M_f / fd.mce, 1.0)
    eta_M = 1.0 - 2.59 * rM + 5.11 * rM^2 - 3.52 * rM^3  # Eq. T2-5
    eta_s = min(0.174 * fd.se^(-0.19), 1.0)  # Eq. T2-4
    I_R = Gamma * w_n * h_us * eta_M * eta_s  # Eq. T2-3

    # Propagating flux ratio — Eq. T2-2
    xi = exp((0.792 + 0.681 * sqrt(sigma)) * (beta + 0.1)) / (192.0 + 0.2595 * sigma)

    # Heat sink — Table 2, Mandel et al. (2011)
    eps = exp(-138.0 / sigma)           # Eq. T2-13
    Q_ig = 250.0 + 1116.0 * M_f        # Eq. T2-14

    # Base spread rate (ft/min) — Eq. T2-1
    R0_us = (I_R * xi) / (rho_b * eps * Q_ig)

    # Wind factor coefficients — Table 2, Mandel et al. (2011)
    C_wind = 7.47 * exp(-0.133 * sigma^0.55)  # Eq. T2-16
    B_wind = 0.02526 * sigma^0.54  # Eq. T2-17
    E_wind = 0.715 * exp(-3.59e-4 * sigma)  # Eq. T2-19

    # Coefficients for Eq. 2 form: S = max{S_0, R_0 + c*min{e, max{0,U}}^b + d*max{0,tan(phi)}^2}
    # In English units: c_us * U_ftmin^b gives ft/min, d_us * tan²φ gives ft/min
    c_wind = R0_us * C_wind * beta_ratio^(-E_wind)  # ft/min (Eq. 2 wind coefficient)
    d_slope_coeff = R0_us * 5.275 * beta^(-0.3)     # ft/min (Eq. 2 slope coefficient)

    # Convert from English (ft/min) to SI (m/s)
    ft_min_to_m_s = 0.3048 / 60.0  # = 0.00508 m/s per ft/min
    m_s_to_ft_min = 60.0 / 0.3048  # = 196.85 ft/min per m/s
    R0_si = R0_us * ft_min_to_m_s
    # For c: c_si * U_ms^b = c_us * (U_ms * m_s_to_ft_min)^b * ft_min_to_m_s
    c_wind_si = c_wind * m_s_to_ft_min^B_wind * ft_min_to_m_s
    # For d: slope tan²φ is dimensionless, so only the rate needs conversion
    d_slope_si = d_slope_coeff * ft_min_to_m_s
    e_max_si = 30.0  # m/s (maximum wind speed, following CAWFE convention)

    # Fuel burn time constant — Mandel et al. (2011)
    T_f = fd.weight / 0.8514

    # Heat content in SI: convert BTU/lb to J/kg
    # 1 BTU = 1055.06 J, 1 lb = 0.453592 kg
    h_si = fd.h * 1055.06 / 0.453592  # J/kg

    return (
        S_0 = 0.0,
        R_0 = R0_si,
        b = B_wind,
        c = c_wind_si,
        d = d_slope_si,
        e = e_max_si,
        T_f = T_f,
        w_l = fd.fgi,
        h = h_si,
        delta = fd.depth,
    )
end
