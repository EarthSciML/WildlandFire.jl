module WildlandFire

using ModelingToolkit
using ModelingToolkit: t, D
using DynamicQuantities

include("rothermel.jl")
include("fire_spread_direction.jl")

end # module
