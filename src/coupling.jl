export RothermelCoupler, TerrainSlope, TerrainSlopeCoupler, MidflameWind, MidflameWindCoupler
export FuelModelLookup, FuelModelLookupCoupler
export EMCCoupler, OneHourFuelMoistureCoupler
export LevelSetCoupler
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

struct FireSpreadDirectionCoupler
    sys::Any
end

# ---- Fuel model lookup functions (registered for symbolic use) ---------------

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
fuel_heat(code) = _lookup_fuel(_FUEL_HEAT, code, 8000.0 * 2326.0)

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
    end

    eqs = [
        tanϕ ~ sqrt(dzdx^2 + dzdy^2) + zero_1,
        slope_aspect ~ atan(dzdy, dzdx) + zero_rad,
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
    end

    eqs = [
        U ~ wind_reduction * sqrt(u_wind^2 + v_wind^2) + zero_ms,
        ω ~ atan(v_wind, u_wind) - slope_aspect + zero_rad,
    ]

    return System(
        eqs, t; name,
        metadata = Dict(CoupleType => MidflameWindCoupler)
    )
end

# ---- couple2 methods (inter-component, no EarthSciData dependency) -----------

# FuelModelLookup → RothermelFireSpread (fuel parameters)
function couple2(fm::FuelModelLookupCoupler, r::RothermelCoupler)
    fm, r = fm.sys, r.sys
    r = param_to_var(r, :σ, :w0, :δ, :Mx, :h)
    return ConnectorSystem(
        [
            r.σ ~ fm.σ,
            r.w0 ~ fm.w0,
            r.δ ~ fm.δ,
            r.Mx ~ fm.Mx,
            r.h ~ fm.h,
        ], r, fm
    )
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

# FireSpreadDirection → LevelSetFireSpread (R_H, Z, α → level-set parameters)
function couple2(fsd::FireSpreadDirectionCoupler, ls::LevelSetCoupler)
    fsd, ls = fsd.sys, ls.sys
    ls = param_to_var(ls, :R_H, :Z, :α)
    # Find the promoted variables from the substituted equation
    eq_vars = collect(Symbolics.get_variables(equations(ls)[1]))
    R_H_sym = only(filter(v -> Symbolics.tosymbol(v, escape = false) == :R_H, eq_vars))
    Z_sym = only(filter(v -> Symbolics.tosymbol(v, escape = false) == :Z, eq_vars))
    α_sym = only(filter(v -> Symbolics.tosymbol(v, escape = false) == :α, eq_vars))
    return ConnectorSystem(
        [R_H_sym ~ fsd.R_H, Z_sym ~ fsd.Z, α_sym ~ fsd.α],
        ls, fsd
    )
end
