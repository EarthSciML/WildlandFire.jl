"""
    WildlandFire

Wildland fire modeling components for EarthSciML, including:
- Rothermel surface fire spread model
- National Fire Danger Rating System (NFDRS) fuel moisture and fire behavior models
"""
module WildlandFire

using DynamicQuantities
using ModelingToolkit
using ModelingToolkit: t, D
using OrdinaryDiffEqDefault

# Include implementations
include("rothermel.jl")
include("fire_spread_direction.jl")
include("nfdrs.jl")

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

end # module
