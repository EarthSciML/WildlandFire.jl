"""
    WildlandFire

Wildland fire modeling components for EarthSciML, including:
- Rothermel surface fire spread model
- National Fire Danger Rating System (NFDRS) fuel moisture and fire behavior models
"""
module WildlandFire

using DomainSets
using DynamicQuantities
using MethodOfLines
using ModelingToolkit
using ModelingToolkit: t, D
using OrdinaryDiffEqDefault
using EarthSciMLBase
using EarthSciData

# Include implementations
include("rothermel.jl")
include("fire_spread_direction.jl")
include("nfdrs.jl")
include("level_set_fire_spread.jl")
include("fsim.jl")
include("coupling.jl")

# Export NFDRS components
export EquilibriumMoistureContent
export OneHourFuelMoisture, TenHourFuelMoisture
export HundredHourFuelMoisture, ThousandHourFuelMoisture
export HerbaceousFuelMoisture, WoodyFuelMoisture
export FuelLoadingTransfer
export SpreadComponent, EnergyReleaseComponent
export BurningIndex, IgnitionComponent
export HumanFireOccurrenceIndex, LightningFireOccurrenceIndex
export FireLoadIndex

# Export fuel model utilities
export NFDRSFuelModel, NFDRS_FUEL_MODELS, get_fuel_model
export fuel_loading_to_kg_per_sqm, fuel_loading_to_lb_per_sqft  # latter deprecated

# Export Rothermel model
export RothermelFireSpread

# Export level-set fire spread model
export LevelSetFireSpread, FuelConsumption, FireHeatFlux
export anderson_fuel_coefficients, ANDERSON_FUEL_DATA
# Export FSim components
export FireOccurrenceLogistic, FireContainment
export BurnProbability, ERCTimeSeries, FlameLengthCategory

end # module
