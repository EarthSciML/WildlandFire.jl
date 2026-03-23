module WildlandFireEarthSciDataExt

using WildlandFire
using EarthSciMLBase
using EarthSciMLBase: ConnectorSystem, param_to_var
import EarthSciMLBase: couple2
using EarthSciData

# LANDFIRE fuel_model → FuelModelLookup
function couple2(lf::EarthSciData.LANDFIRECoupler, fm::WildlandFire.FuelModelLookupCoupler)
    lf, fm = lf.sys, fm.sys
    fm = param_to_var(fm, :fuel_model)
    return ConnectorSystem([fm.fuel_model ~ lf.fuel_model], fm, lf)
end

# USGS3DEP dzdx/dzdy → TerrainSlope
function couple2(dep::EarthSciData.USGS3DEPCoupler, ts::WildlandFire.TerrainSlopeCoupler)
    dep, ts = dep.sys, ts.sys
    ts = param_to_var(ts, :dzdx, :dzdy)
    return ConnectorSystem(
        [
            ts.dzdx ~ dep.dzdx,
            ts.dzdy ~ dep.dzdy,
        ], ts, dep
    )
end

# ERA5 → MidflameWind (wind components)
function couple2(era::EarthSciData.ERA5Coupler, mw::WildlandFire.MidflameWindCoupler)
    era, mw = era.sys, mw.sys
    mw = param_to_var(mw, :u_wind, :v_wind)
    return ConnectorSystem(
        [
            mw.u_wind ~ era.pl₊u,
            mw.v_wind ~ era.pl₊v,
        ], mw, era
    )
end

# ERA5 → EquilibriumMoistureContent (temperature and relative humidity)
function couple2(era::EarthSciData.ERA5Coupler, emc::WildlandFire.EMCCoupler)
    era, emc = era.sys, emc.sys
    emc = param_to_var(emc, :TEMP, :RH)
    return ConnectorSystem(
        [
            emc.TEMP ~ era.pl₊t,
            emc.RH ~ era.pl₊r,
        ], emc, era
    )
end

end # module
