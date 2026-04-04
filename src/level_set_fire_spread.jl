"""
    LevelSetFireSpread(domain::DomainInfo; name=:LevelSetFireSpread,
                       initial_condition, boundary_conditions=nothing,
                       spread_rate=1.0)

Create a level-set fire front propagation PDE system using the Mandel et al. (2011) /
Muñoz-Esparza et al. (2018) normal-projection approach.

The fire front is tracked implicitly as the zero contour of a level-set function
ψ(x, y, t), where ψ ≤ 0 is the burning region and ψ > 0 is unburned fuel. The
level-set evolves according to the Hamilton-Jacobi equation (Mandel 2011, Eq. 9):

```math
\\frac{\\partial \\psi}{\\partial t} + S(\\hat{n}) \\|\\nabla \\psi\\| = 0
```

The direction-dependent spread rate is evaluated by projecting wind and terrain slope
onto the fire front normal ``\\hat{n} = \\nabla\\psi / \\|\\nabla\\psi\\|``
(Mandel 2011, Eq. 2; Muñoz-Esparza 2018, Eq. 10):

```math
S(\\hat{n}) = R_0 \\left(1 + C \\left(\\frac{\\max(0,\\, \\mathbf{u}_f \\cdot \\hat{n})}{U_{\\text{ref}}}\\right)^B \\beta_{\\text{ratio}}^{-E}
  + \\phi_{s,\\text{coeff}} \\, \\max(0,\\, \\nabla z \\cdot \\hat{n})^2 \\right)
```

When wind and slope are zero the equation reduces to the isotropic case ``S = R_0``.

## Implementation Details

Uses MethodOfLines.jl with WENO spatial discretization following Muñoz-Esparza et al.
(2018). For production use, consider adding level-set reinitialization.

# Arguments
- `domain`: A `DomainInfo` with at least 2 spatial dimensions.
- `name`: System name (default `:LevelSetFireSpread`)
- `initial_condition`: Function `(x, y) -> ψ₀`. Negative = burning, positive = unburned.
- `boundary_conditions`: Optional BC equations. Default: Neumann (zero-gradient).
- `spread_rate`: Default no-wind no-slope rate of spread R_0 (m/s). Default 1.0.

# References

Mandel, J., Beezley, J.D., and Kochanski, A.K. (2011). Coupled atmosphere-wildland fire
modeling with WRF 3.3 and SFIRE 2011. *Geosci. Model Dev.*, 4, 591–610.
doi:10.5194/gmd-4-591-2011

Muñoz-Esparza, D., Kosović, B., Jiménez, P.A., and Coen, J.L. (2018). An accurate
fire-spread algorithm in the Weather Research and Forecasting model using the level-set
method. *J. Adv. Model. Earth Syst.*, 10, 908–926. doi:10.1002/2017MS001108

# Example

```julia
using WildlandFire, EarthSciMLBase, ModelingToolkit, DynamicQuantities
using ModelingToolkit: t
using MethodOfLines, DomainSets, OrdinaryDiffEqDefault

@parameters x [unit = u"m"]
@parameters y [unit = u"m"]
domain = DomainInfo(
    constIC(0.0, t ∈ Interval(0.0, 60.0)),
    constBC(0.0, x ∈ Interval(0.0, 500.0), y ∈ Interval(0.0, 500.0)),
)
initial_condition(x, y) = sqrt((x - 250.0)^2 + (y - 250.0)^2) - 10.0
sys = LevelSetFireSpread(domain; initial_condition, spread_rate=0.5)

dx = 5.0
discretization = MOLFiniteDifference([sys.ivs[2] => dx, sys.ivs[3] => dx], sys.ivs[1];
    advection_scheme = WENOScheme())
prob = MethodOfLines.discretize(sys, discretization; checks=false)
sol = solve(prob)
```
"""
function LevelSetFireSpread(
        domain::DomainInfo;
        name = :LevelSetFireSpread,
        initial_condition,
        boundary_conditions = nothing,
        spread_rate = 1.0
    )

    # Extract domain information
    t_domain = get_tspan(domain)
    spatial_eps = EarthSciMLBase.endpoints(domain)
    @assert length(spatial_eps) >= 2 "LevelSetFireSpread requires at least 2 spatial dimensions, got $(length(spatial_eps))"
    x_domain = spatial_eps[1]
    y_domain = spatial_eps[2]

    # Get spatial variables from domain (use only the first 2 for the 2D fire PDE)
    spatial_vars = EarthSciMLBase.pvars(domain)
    x = spatial_vars[1]
    y = spatial_vars[2]

    # Rothermel fire spread coefficients — Mandel (2011) Eq. 2
    @parameters R_0 = spread_rate [description = "No-wind no-slope rate of spread", unit = u"m/s"]
    @parameters C_wind = 1.0 [description = "Wind coefficient C (dimensionless)", unit = u"1"]
    @parameters B_wind = 0.5 [description = "Wind exponent B (dimensionless)", unit = u"1"]
    @parameters E_wind = 0.5 [description = "Wind coefficient E (dimensionless)", unit = u"1"]
    @parameters β_ratio = 1.0 [description = "Relative packing ratio β/β_op (dimensionless)", unit = u"1"]
    @parameters φs_coeff = 0.0 [description = "Slope factor coefficient 5.275*β^(-0.3) (dimensionless)", unit = u"1"]

    # Wind vector at midflame height
    @parameters u_x = 0.0 [description = "Midflame wind x-component", unit = u"m/s"]
    @parameters u_y = 0.0 [description = "Midflame wind y-component", unit = u"m/s"]

    # Terrain gradient
    @parameters dzdx = 0.0 [description = "Terrain gradient in x direction (dimensionless)", unit = u"1"]
    @parameters dzdy = 0.0 [description = "Terrain gradient in y direction (dimensionless)", unit = u"1"]

    @constants begin
        one = 1.0, [description = "Dimensionless one for unit balancing", unit = u"1"]
        ψ_ref = 1.0, [description = "Reference length for initial condition", unit = u"m"]
        grad_eps = 1.0e-10, [description = "Small value for numerical stability (dimensionless)", unit = u"1"]
        # Reference wind speed: 1 ft/min in m/s — matches Rothermel calibration units
        U_ref = 0.3048 / 60, [description = "Reference wind speed (1 ft/min in m/s)", unit = u"m/s"]
        zero_ms = 0.0, [description = "Zero wind speed", unit = u"m/s"]
        zero_1 = 0.0, [description = "Zero dimensionless", unit = u"1"]
        β_ratio_floor = 1.0e-10, [description = "Minimum β_ratio to avoid singularity (dimensionless)", unit = u"1"]
    end

    # Level-set function
    @variables ψ(..) [description = "Level-set function (fire front at ψ=0)", unit = u"m"]

    # Spatial derivative operators with coordinate transforms
    δs = EarthSciMLBase.partialderivatives(domain)
    transforms = EarthSciMLBase.partialderivative_transforms(domain)
    Dx = Differential(x)
    Dy = Differential(y)

    # Collect symbolic constants/parameters introduced by coordinate transforms
    ivs_set = Set(Symbolics.unwrap.([t, x, y]))
    transform_params_set = Set{Any}()
    transform_params = Num[]
    for tf in transforms
        tf isa Integer && continue
        for v in Symbolics.get_variables(tf)
            uv = Symbolics.unwrap(v)
            if uv ∉ ivs_set && uv ∉ transform_params_set
                push!(transform_params_set, uv)
                push!(transform_params, v)
            end
        end
    end

    # Gradient components — Eq. 6, Mandel et al. (2011)
    ψ_x = δs[1](ψ(t, x, y))
    ψ_y = δs[2](ψ(t, x, y))

    # Gradient magnitude with regularization to avoid division by zero
    grad_ψ = sqrt(ψ_x^2 + ψ_y^2)
    grad_ψ_safe = grad_ψ + grad_eps

    # Normal-projection fire spread rate — Mandel (2011) Eq. 2, Esparza (2018) Eq. 10
    #
    # Fire front normal: n = ∇ψ / |∇ψ|
    # Normal wind component: U_n = u_f · n
    # Normal slope component: tanφ_n = ∇z · n
    #
    # Wind factor: φ_W(n) = C * (max(0, U_n) / U_ref)^B * β_ratio^(-E)
    # Slope factor: φ_S(n) = φs_coeff * max(0, tanφ_n)²
    # Spread rate: S(n) = R_0 * (1 + φ_W(n) + φ_S(n))

    # Project wind and slope onto fire front normal
    U_n = (u_x * ψ_x + u_y * ψ_y) / grad_ψ_safe + zero_ms
    tanφ_n = (dzdx * ψ_x + dzdy * ψ_y) / grad_ψ_safe + zero_1

    # Wind factor in normal direction (Rothermel form)
    φ_W_n = C_wind * (max(zero_ms, U_n) / U_ref)^B_wind * max(β_ratio, β_ratio_floor)^(-E_wind)

    # Slope factor in normal direction
    φ_S_n = φs_coeff * max(zero_1, tanφ_n)^2

    # Total spread rate — Mandel (2011) Eq. 1
    S_n = R_0 * (one + φ_W_n + φ_S_n)

    # Level-set equation: ∂ψ/∂t + S(n)|∇ψ| = 0 — Mandel (2011) Eq. 9
    eq = [
        D(ψ(t, x, y)) ~ -S_n * grad_ψ,
    ]

    pde_domains = EarthSciMLBase.domains(domain)

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
        eq, bcs, pde_domains, [t, x, y], [ψ(t, x, y)],
        [
            R_0, C_wind, B_wind, E_wind, β_ratio, φs_coeff,
            u_x, u_y, dzdx, dzdy,
            ψ_ref, one, grad_eps, U_ref, zero_ms, zero_1, β_ratio_floor,
            transform_params...,
        ]; name = name,
        metadata = Dict(EarthSciMLBase.CoupleType => LevelSetCoupler)
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

The effective fuel load is computed as ``w_{0,\\text{eff}} = F \\cdot w_{0,\\text{initial}}``,
which can be coupled to `RothermelFireSpread` to reduce the fire spread rate as fuel
is consumed.

When coupled to `LevelSetFireSpread`, the `is_burning` parameter is driven by the
level-set function ψ via a smooth Heaviside approximation:
``\\text{is\\_burning} = 0.5 (1 - \\tanh(\\psi / \\varepsilon))``

# Reference

Mandel, J., Beezley, J.D., and Kochanski, A.K. (2011). Coupled atmosphere-wildland fire
modeling with WRF 3.3 and SFIRE 2011. *Geosci. Model Dev.*, 4, 591–610.
doi:10.5194/gmd-4-591-2011
"""
@component function FuelConsumption(; name = :FuelConsumption)
    @parameters begin
        T_f, [description = "Fuel burn time constant (1/e decay time)", unit = u"s"]
        w0_initial, [description = "Initial oven-dry fuel load", unit = u"kg/m^2"]
        is_burning = 0.0, [description = "Whether the fuel is burning (1=yes, 0=no) (dimensionless)", unit = u"1"]
    end

    @variables begin
        F(t), [description = "Fuel fraction remaining (1=full, 0=consumed)", unit = u"1"]
        w0_effective(t), [description = "Effective fuel load after consumption", unit = u"kg/m^2"]
    end

    eqs = [
        D(F) ~ -is_burning * F / T_f,  # Eq. 3, Mandel et al. (2011), differential form
        w0_effective ~ F * w0_initial,   # Fuel load feedback
    ]

    return System(
        eqs, t; name,
        metadata = Dict(CoupleType => FuelConsumptionCoupler)
    )
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

    # Unit conversion constants - numerical for computational efficiency
    ft_to_m = 0.3048        # feet to meters conversion
    lb_to_kg = 0.453592     # pounds to kilograms conversion
    btu_to_j = 1055.06      # BTU to Joules conversion
    kgm2_to_lbft2 = 0.204816  # kg/m² to lb/ft² conversion
    min_to_s = 60.0         # minutes to seconds conversion
    one = 1.0               # dimensionless one

    # Convert inputs from SI to English units
    fgi_lbft2 = fd.fgi * kgm2_to_lbft2     # kg/m^2 → lb/ft^2
    depth_ft = fd.depth / ft_to_m          # m → ft
    sigma = fd.savr                        # 1/ft (already in English units)
    rho_p = fd.dens                        # lb/ft^3 (already in English units)
    h_us = fd.h                            # BTU/lb (heat content from fuel model data)

    # Empirical constants from Rothermel (1972) formulas — Table 2, Mandel et al. (2011)
    # Eq. T2-21 coefficients for optimum packing ratio
    beta_op_coeff = 3.348
    beta_op_exp = -0.8189

    # Eq. T2-12 coefficients
    coeff_A_c1 = 4.774
    coeff_A_exp = 0.1
    coeff_A_c2 = 7.27

    # Eq. T2-9 coefficients for maximum reaction velocity
    gamma_c1 = 495.0
    gamma_c2 = 0.0594
    gamma_exp = 1.5

    # Eq. T2-5 coefficients for moisture damping
    eta_M_c1 = 2.59
    eta_M_c2 = 5.11
    eta_M_c3 = 3.52

    # Eq. T2-4 coefficients for mineral damping
    eta_s_c1 = 0.174
    eta_s_exp = -0.19

    # Eq. T2-2 coefficients for propagating flux ratio
    xi_c1 = 0.792
    xi_c2 = 0.681
    xi_c3 = 0.1
    xi_c4 = 192.0
    xi_c5 = 0.2595

    # Fuel bed properties — Table 2, Mandel et al. (2011)
    w_0 = fgi_lbft2 / (one + M_f)     # Total fuel load net of moisture (Eq. T2-7)
    w_n = fgi_lbft2 / (one + fd.st)   # Fuel loading net of minerals (Eq. T2-6)
    rho_b = w_0 / depth_ft            # Oven dry bulk density (Eq. T2-11)
    beta = rho_b / rho_p              # Packing ratio (Eq. T2-10)
    beta_op = beta_op_coeff * sigma^beta_op_exp # Optimum packing ratio (Eq. T2-21)
    beta_ratio = beta / beta_op

    # Reaction intensity components — Table 2, Mandel et al. (2011)
    coeff_A = one / (coeff_A_c1 * sigma^coeff_A_exp - coeff_A_c2)  # Eq. T2-12
    Gamma_max = sigma^gamma_exp / (gamma_c1 + gamma_c2 * sigma^gamma_exp)  # Eq. T2-9
    Gamma = Gamma_max * beta_ratio^coeff_A * exp(coeff_A * (one - beta_ratio))  # Eq. T2-8

    rM = min(M_f / fd.mce, one)
    eta_M = one - eta_M_c1 * rM + eta_M_c2 * rM^2 - eta_M_c3 * rM^3  # Eq. T2-5
    eta_s = min(eta_s_c1 * fd.se^eta_s_exp, one)  # Eq. T2-4
    I_R = Gamma * w_n * h_us * eta_M * eta_s  # Eq. T2-3

    # Propagating flux ratio — Eq. T2-2
    xi = exp((xi_c1 + xi_c2 * sqrt(sigma)) * (beta + xi_c3)) / (xi_c4 + xi_c5 * sigma)

    # More empirical constants from Rothermel (1972) formulas
    # Eq. T2-13 heat sink coefficient
    eps_coeff = -138.0

    # Eq. T2-14 ignition energy coefficients
    Q_ig_c1 = 250.0
    Q_ig_c2 = 1116.0

    # Eq. T2-16 wind factor coefficient
    C_wind_coeff = 7.47
    C_wind_exp_coeff = -0.133
    C_wind_exp = 0.55

    # Eq. T2-17 wind speed exponent
    B_wind_coeff = 0.02526
    B_wind_exp = 0.54

    # Eq. T2-19 wind speed coefficient
    E_wind_coeff = 0.715
    E_wind_exp_coeff = -3.59e-4

    # Slope coefficient
    slope_coeff = 5.275
    slope_exp = -0.3

    # Unit conversion rate constants
    ft_min_to_m_s = 0.00508  # ft/min to m/s conversion
    m_s_to_ft_min = 196.85   # m/s to ft/min conversion

    # Heat sink — Table 2, Mandel et al. (2011)
    eps = exp(eps_coeff / sigma)           # Eq. T2-13
    Q_ig = Q_ig_c1 + Q_ig_c2 * M_f        # Eq. T2-14

    # Base spread rate (ft/min) — Eq. T2-1
    R0_us = (I_R * xi) / (rho_b * eps * Q_ig)

    # Wind factor coefficients — Table 2, Mandel et al. (2011)
    C_wind = C_wind_coeff * exp(C_wind_exp_coeff * sigma^C_wind_exp)  # Eq. T2-16
    B_wind = B_wind_coeff * sigma^B_wind_exp  # Eq. T2-17
    E_wind = E_wind_coeff * exp(E_wind_exp_coeff * sigma)  # Eq. T2-19

    # Coefficients for Eq. 2 form: S = max{S_0, R_0 + c*min{e, max{0,U}}^b + d*max{0,tan(phi)}^2}
    # In English units: c_us * U_ftmin^b gives ft/min, d_us * tan²φ gives ft/min
    c_wind = R0_us * C_wind * beta_ratio^(-E_wind)  # ft/min (Eq. 2 wind coefficient)
    d_slope_coeff = R0_us * slope_coeff * beta^slope_exp     # ft/min (Eq. 2 slope coefficient)

    # Convert from English (ft/min) to SI (m/s)
    R0_si = R0_us * ft_min_to_m_s
    # More unit conversion constants
    # Maximum wind speed (from CAWFE convention)
    e_max = 30.0  # m/s

    # Fuel burn time conversion factor
    fuel_burn_factor = 0.8514

    # Heat content conversion to SI: 1 BTU/lb = 2326.11 J/kg
    btu_lb_to_j_kg = 2326.11  # J/kg per BTU/lb

    # For c: c_si * U_ms^b = c_us * (U_ms * m_s_to_ft_min)^b * ft_min_to_m_s
    c_wind_si = c_wind * m_s_to_ft_min^B_wind * ft_min_to_m_s
    # For d: slope tan²φ is dimensionless, so only the rate needs conversion
    d_slope_si = d_slope_coeff * ft_min_to_m_s
    e_max_si = e_max  # m/s (maximum wind speed, following CAWFE convention)

    # Fuel burn time constant — Mandel et al. (2011)
    T_f = fd.weight / fuel_burn_factor

    # Heat content in SI: convert BTU/lb to J/kg
    h_si = fd.h * btu_lb_to_j_kg  # J/kg

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
