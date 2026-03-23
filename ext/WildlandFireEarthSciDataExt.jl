module WildlandFireEarthSciDataExt

using WildlandFire
using EarthSciMLBase
using EarthSciMLBase: ConnectorSystem, param_to_var
import EarthSciMLBase: couple2
using EarthSciData

# Only define coupling methods for coupler types that exist in the loaded
# version of EarthSciData. This allows the extension to work with both
# released and development versions.

if isdefined(EarthSciData, :LANDFIRECoupler)
    @eval begin
        # LANDFIRE fuel_model → FuelModelLookup
        function couple2(lf::EarthSciData.LANDFIRECoupler, fm::WildlandFire.FuelModelLookupCoupler)
            lf, fm = lf.sys, fm.sys
            fm = param_to_var(fm, :fuel_model)
            return ConnectorSystem([fm.fuel_model ~ lf.fuel_model], fm, lf)
        end
    end
end

if isdefined(EarthSciData, :ERA5Coupler)
    @eval begin
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
    end
end

end # module
