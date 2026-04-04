export RothermelCoupler, TerrainSlope, TerrainSlopeCoupler, MidflameWind, MidflameWindCoupler
export FuelModelLookup, FuelModelLookupCoupler
export EMCCoupler, OneHourFuelMoistureCoupler
export LevelSetCoupler, FuelConsumptionCoupler
export FireSpreadDirectionCoupler

using EarthSciMLBase
using EarthSciMLBase: CoupleType, ConnectorSystem, param_to_var
import EarthSciMLBase: couple2

# ---- Coupler structs ---------------------------------------------------------

struct RothermelCoupler
    sys::Any
end

struct TerrainSlopeCoupler
    sys::Any
end

struct MidflameWindCoupler
    sys::Any
end

struct FuelModelLookupCoupler
    sys::Any
end

struct EMCCoupler
    sys::Any
end

struct OneHourFuelMoistureCoupler
    sys::Any
end

struct LevelSetCoupler
    sys::Any
end

struct FuelConsumptionCoupler
    sys::Any
end

struct FireSpreadDirectionCoupler
    sys::Any
end

# ---- Fuel model lookup functions (registered for symbolic use) ---------------

# LANDFIRE non-burnable fuel model codes (urban, snow/ice, agriculture, water, barren).
const _NONBURNABLE_CODES = Set([91, 92, 93, 98, 99])

# Build lookup tables from ANDERSON_FUEL_DATA once at module load time.
const _FUEL_SAVR = let
    d = Dict{Int, Float64}()
    for (k, v) in ANDERSON_FUEL_DATA
        # ANDERSON_FUEL_DATA.savr is in 1/ft; convert to 1/m
        d[k] = v.savr / 0.3048
    end
    d
end

const _FUEL_LOAD = let
    d = Dict{Int, Float64}()
    for (k, v) in ANDERSON_FUEL_DATA
        d[k] = v.fgi  # already kg/m²
    end
    d
end

const _FUEL_DEPTH = let
    d = Dict{Int, Float64}()
    for (k, v) in ANDERSON_FUEL_DATA
        d[k] = v.depth  # already m
    end
    d
end

const _FUEL_MCE = let
    d = Dict{Int, Float64}()
    for (k, v) in ANDERSON_FUEL_DATA
        d[k] = v.mce  # dimensionless
    end
    d
end

const _FUEL_HEAT = let
    d = Dict{Int, Float64}()
    for (k, v) in ANDERSON_FUEL_DATA
        # ANDERSON_FUEL_DATA.h is in BTU/lb; convert to J/kg
        d[k] = v.h * 2326.0
    end
    for code in _NONBURNABLE_CODES
        d[code] = 0.0
    end
    d
end

function _lookup_fuel(table::Dict{Int, Float64}, code::Real, default::Float64)
    c = round(Int, code)
    return get(table, c, default)
end

# Wrapper functions for each fuel property.
fuel_savr(code) = _lookup_fuel(_FUEL_SAVR, code, 3500.0 / 0.3048)
fuel_load(code) = _lookup_fuel(_FUEL_LOAD, code, 0.166)
fuel_depth(code) = _lookup_fuel(_FUEL_DEPTH, code, 0.305)
fuel_mce(code) = _lookup_fuel(_FUEL_MCE, code, 0.12)
fuel_heat(code) = _lookup_fuel(_FUEL_HEAT, code, 0.0)

# Register for use in symbolic equations.
@register_symbolic fuel_savr(code)
@register_symbolic fuel_load(code)
@register_symbolic fuel_depth(code)
@register_symbolic fuel_mce(code)
@register_symbolic fuel_heat(code)

# Unit-validation pass-throughs for ModelingToolkit.
fuel_savr(::DynamicQuantities.AbstractQuantity) = 1.0
fuel_load(::DynamicQuantities.AbstractQuantity) = 1.0
fuel_depth(::DynamicQuantities.AbstractQuantity) = 1.0
fuel_mce(::DynamicQuantities.AbstractQuantity) = 1.0
fuel_heat(::DynamicQuantities.AbstractQuantity) = 1.0

# ---- FuelModelLookup component -----------------------------------------------

"""
    FuelModelLookup(; name = :FuelModelLookup)

A component that maps an integer fuel model code to continuous Rothermel fuel
parameters using the Anderson 13 fuel model data.

Takes `fuel_model` as input (dimensionless integer code from LANDFIRE) and
outputs the five core Rothermel parameters: `σ`, `w0`, `δ`, `Mx`, `h`.

Non-burnable LANDFIRE codes (91=urban, 92=snow/ice, 93=agriculture, 98=water,
99=barren) and unrecognized codes return zero heat content (`h=0`), which
guarantees the Rothermel model computes zero reaction intensity (`IR=0`) and
therefore zero rate of spread (`R=0`).
"""
@component function FuelModelLookup(; name = :FuelModelLookup)
    @constants begin
        σ_unit = 1.0, [description = "SAV ratio unit", unit = u"1/m"]
        w0_unit = 1.0, [description = "Fuel load unit", unit = u"kg/m^2"]
        δ_unit = 1.0, [description = "Fuel depth unit", unit = u"m"]
        h_unit = 1.0, [description = "Heat content unit", unit = u"J/kg"]
        one = 1.0, [description = "Dimensionless one", unit = u"1"]
    end

    @parameters begin
        fuel_model, [description = "Anderson 13 fuel model code (dimensionless)", unit = u"1"]
    end

    @variables begin
        σ(t), [description = "Surface-area-to-volume ratio", unit = u"1/m"]
        w0(t), [description = "Oven-dry fuel load", unit = u"kg/m^2"]
        δ(t), [description = "Fuel bed depth", unit = u"m"]
        Mx(t), [description = "Dead fuel moisture of extinction (dimensionless)", unit = u"1"]
        h(t), [description = "Low heat content", unit = u"J/kg"]
    end

    eqs = [
        σ ~ fuel_savr(fuel_model) * σ_unit,
        w0 ~ fuel_load(fuel_model) * w0_unit,
        δ ~ fuel_depth(fuel_model) * δ_unit,
        Mx ~ fuel_mce(fuel_model) * one,
        h ~ fuel_heat(fuel_model) * h_unit,
    ]

    return System(
        eqs, t; name,
        metadata = Dict(CoupleType => FuelModelLookupCoupler)
    )
end

# ---- TerrainSlope component --------------------------------------------------

"""
    TerrainSlope(; name = :TerrainSlope)

A component that computes terrain slope magnitude and aspect from elevation
gradients.

Outputs:
- `tanϕ`: slope steepness (rise/run, dimensionless)
- `slope_aspect`: angle of steepest descent from east (rad)
"""
@component function TerrainSlope(; name = :TerrainSlope)
    @constants begin
        zero_rad = 0.0, [description = "Zero angle", unit = u"rad"]
        zero_1 = 0.0, [description = "Zero dimensionless", unit = u"1"]
    end

    @parameters begin
        dzdx = 0.0, [description = "Elevation gradient in x (east) direction (dimensionless)", unit = u"1"]
        dzdy = 0.0, [description = "Elevation gradient in y (north) direction (dimensionless)", unit = u"1"]
    end

    @variables begin
        tanϕ(t), [description = "Slope steepness (dimensionless, rise/run)", unit = u"1"]
        slope_aspect(t), [description = "Aspect angle from east (upslope direction)", unit = u"rad"]
        dzdx_out(t), [description = "Elevation gradient in x direction (output)", unit = u"1"]
        dzdy_out(t), [description = "Elevation gradient in y direction (output)", unit = u"1"]
    end

    eqs = [
        tanϕ ~ sqrt(dzdx^2 + dzdy^2) + zero_1,
        slope_aspect ~ atan(dzdy, dzdx) + zero_rad,
        dzdx_out ~ dzdx + zero_1,
        dzdy_out ~ dzdy + zero_1,
    ]

    return System(
        eqs, t; name,
        metadata = Dict(CoupleType => TerrainSlopeCoupler)
    )
end

# ---- MidflameWind component --------------------------------------------------

"""
    MidflameWind(; name = :MidflameWind)

A component that converts horizontal wind components (u, v) to Rothermel's
midflame-height wind speed `U` and wind direction `ω` relative to the
upslope direction.

Uses a wind reduction factor to convert 10m wind to midflame height (~6.1m
for forest fuels). The default reduction factor of 0.4 is typical for
timber fuel types (Baughman & Albini, 1980).
"""
@component function MidflameWind(; name = :MidflameWind)
    @constants begin
        zero_ms = 0.0, [description = "Zero wind speed", unit = u"m/s"]
        zero_rad = 0.0, [description = "Zero angle", unit = u"rad"]
    end

    @parameters begin
        u_wind = 0.0, [description = "Eastward wind component at 10m", unit = u"m/s"]
        v_wind = 0.0, [description = "Northward wind component at 10m", unit = u"m/s"]
        slope_aspect = 0.0, [description = "Upslope direction from east", unit = u"rad"]
        wind_reduction = 0.4, [description = "Wind reduction factor (10m to midflame)", unit = u"1"]
    end

    @variables begin
        U(t), [description = "Wind speed at midflame height", unit = u"m/s"]
        ω(t), [description = "Wind direction relative to upslope", unit = u"rad"]
        u_mf_x(t), [description = "Midflame wind x-component (eastward)", unit = u"m/s"]
        u_mf_y(t), [description = "Midflame wind y-component (northward)", unit = u"m/s"]
    end

    eqs = [
        U ~ wind_reduction * sqrt(u_wind^2 + v_wind^2) + zero_ms,
        ω ~ atan(v_wind, u_wind) - slope_aspect + zero_rad,
        u_mf_x ~ wind_reduction * u_wind + zero_ms,
        u_mf_y ~ wind_reduction * v_wind + zero_ms,
    ]

    return System(
        eqs, t; name,
        metadata = Dict(CoupleType => MidflameWindCoupler)
    )
end

# ---- couple2 methods (inter-component, no EarthSciData dependency) -----------

# FuelModelLookup → RothermelFireSpread (fuel parameters except w0, which is routed
# through FuelConsumption for fuel depletion feedback)
function couple2(fm::FuelModelLookupCoupler, r::RothermelCoupler)
    fm, r = fm.sys, r.sys
    r = param_to_var(r, :σ, :δ, :Mx, :h)
    return ConnectorSystem(
        [
            r.σ ~ fm.σ,
            r.δ ~ fm.δ,
            r.Mx ~ fm.Mx,
            r.h ~ fm.h,
        ], r, fm
    )
end

# FuelModelLookup → FuelConsumption (w0 → w0_initial)
# Routes the fuel load through FuelConsumption so it can be scaled by fuel depletion.
function couple2(fm::FuelModelLookupCoupler, fc::FuelConsumptionCoupler)
    fm, fc = fm.sys, fc.sys
    fc = param_to_var(fc, :w0_initial)
    return ConnectorSystem([fc.w0_initial ~ fm.w0], fc, fm)
end

# TerrainSlope → RothermelFireSpread (slope)
function couple2(ts::TerrainSlopeCoupler, r::RothermelCoupler)
    ts, r = ts.sys, r.sys
    r = param_to_var(r, :tanϕ)
    return ConnectorSystem([r.tanϕ ~ ts.tanϕ], r, ts)
end

# MidflameWind → RothermelFireSpread (wind speed)
function couple2(mw::MidflameWindCoupler, r::RothermelCoupler)
    mw, r = mw.sys, r.sys
    r = param_to_var(r, :U)
    return ConnectorSystem([r.U ~ mw.U], r, mw)
end

# TerrainSlope → MidflameWind (slope aspect for wind direction)
function couple2(ts::TerrainSlopeCoupler, mw::MidflameWindCoupler)
    ts, mw = ts.sys, mw.sys
    mw = param_to_var(mw, :slope_aspect)
    return ConnectorSystem([mw.slope_aspect ~ ts.slope_aspect], mw, ts)
end

# EquilibriumMoistureContent → OneHourFuelMoisture (EMC → EMCPRM)
function couple2(emc::EMCCoupler, fm1::OneHourFuelMoistureCoupler)
    emc, fm1 = emc.sys, fm1.sys
    fm1 = param_to_var(fm1, :EMCPRM)
    return ConnectorSystem([fm1.EMCPRM ~ emc.EMC], fm1, emc)
end

# OneHourFuelMoisture → RothermelFireSpread (MC1 → Mf)
function couple2(fm1::OneHourFuelMoistureCoupler, r::RothermelCoupler)
    fm1, r = fm1.sys, r.sys
    r = param_to_var(r, :Mf)
    return ConnectorSystem([r.Mf ~ fm1.MC1], r, fm1)
end

# RothermelFireSpread → FireSpreadDirection (Rothermel outputs → direction inputs)
function couple2(r::RothermelCoupler, fsd::FireSpreadDirectionCoupler)
    r, fsd = r.sys, fsd.sys
    fsd = param_to_var(fsd, :R0, :φw, :φs, :β_ratio, :C_coeff, :B_coeff, :E_coeff)
    return ConnectorSystem(
        [
            fsd.R0 ~ r.R0,
            fsd.φw ~ r.φw,
            fsd.φs ~ r.φs,
            fsd.β_ratio ~ r.β_ratio,
            fsd.C_coeff ~ r.C_coeff,
            fsd.B_coeff ~ r.B_coeff,
            fsd.E_coeff ~ r.E_coeff,
        ], fsd, r
    )
end

# MidflameWind → FireSpreadDirection (wind direction)
function couple2(mw::MidflameWindCoupler, fsd::FireSpreadDirectionCoupler)
    mw, fsd = mw.sys, fsd.sys
    fsd = param_to_var(fsd, :ω)
    return ConnectorSystem([fsd.ω ~ mw.ω], fsd, mw)
end

# RothermelFireSpread → LevelSetFireSpread (Rothermel coefficients for Mandel normal-projection)
function couple2(r::RothermelCoupler, ls::LevelSetCoupler)
    r, ls = r.sys, ls.sys
    ls = param_to_var(ls, :R_0, :C_wind, :B_wind, :E_wind, :β_ratio, :φs_coeff)
    eq_vars = collect(Symbolics.get_variables(equations(ls)[1]))
    find_sym(name) = only(filter(v -> Symbolics.tosymbol(v, escape = false) == name, eq_vars))
    return ConnectorSystem(
        [
            find_sym(:R_0) ~ r.R0,
            find_sym(:C_wind) ~ r.C_coeff,
            find_sym(:B_wind) ~ r.B_coeff,
            find_sym(:E_wind) ~ r.E_coeff,
            find_sym(:β_ratio) ~ r.β_ratio,
            find_sym(:φs_coeff) ~ r.φs_coeff,
        ],
        ls, r
    )
end

# MidflameWind → LevelSetFireSpread (wind vector components)
function couple2(mw::MidflameWindCoupler, ls::LevelSetCoupler)
    mw, ls = mw.sys, ls.sys
    ls = param_to_var(ls, :u_x, :u_y)
    eq_vars = collect(Symbolics.get_variables(equations(ls)[1]))
    find_sym(name) = only(filter(v -> Symbolics.tosymbol(v, escape = false) == name, eq_vars))
    return ConnectorSystem(
        [find_sym(:u_x) ~ mw.u_mf_x, find_sym(:u_y) ~ mw.u_mf_y],
        ls, mw
    )
end

# TerrainSlope → LevelSetFireSpread (terrain gradient components)
function couple2(ts::TerrainSlopeCoupler, ls::LevelSetCoupler)
    ts, ls = ts.sys, ls.sys
    ls = param_to_var(ls, :dzdx, :dzdy)
    eq_vars = collect(Symbolics.get_variables(equations(ls)[1]))
    find_sym(name) = only(filter(v -> Symbolics.tosymbol(v, escape = false) == name, eq_vars))
    return ConnectorSystem(
        [find_sym(:dzdx) ~ ts.dzdx_out, find_sym(:dzdy) ~ ts.dzdy_out],
        ls, ts
    )
end

# LevelSetFireSpread → FuelConsumption (ψ → is_burning)
# Drives the burning state from the level-set function using a smooth Heaviside
# approximation: is_burning = 0.5 * (1 - tanh(ψ/ε)), where ε is a smoothing width.
# When ψ < 0 (burning region), is_burning ≈ 1; when ψ > 0 (unburned), is_burning ≈ 0.
function couple2(ls::LevelSetCoupler, fc::FuelConsumptionCoupler)
    ls, fc = ls.sys, fc.sys
    fc = param_to_var(fc, :is_burning)
    # Extract ψ from the PDE dependent variables (first dv is ψ(t, x, y))
    ψ_sym = ls.dvs[1]
    # Smooth Heaviside: is_burning = 0.5 * (1 - tanh(ψ/ε))
    # ε_h provides meter units to cancel ψ's meters, making the tanh argument dimensionless
    @constants ε_h = 1.0, [description = "Heaviside smoothing width", unit = u"m"]
    return ConnectorSystem(
        [fc.is_burning ~ 0.5 * (1.0 - tanh(ψ_sym / ε_h))],
        fc, ls
    )
end

# FuelConsumption → RothermelFireSpread (w0_effective → w0)
# Scales the Rothermel fuel load by the remaining fuel fraction, causing the fire
# to stop spreading at locations where fuel has been consumed (F → 0).
function couple2(fc::FuelConsumptionCoupler, r::RothermelCoupler)
    fc, r = fc.sys, r.sys
    r = param_to_var(r, :w0)
    return ConnectorSystem([r.w0 ~ fc.w0_effective], r, fc)
end
