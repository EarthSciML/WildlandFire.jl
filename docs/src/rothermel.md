# Rothermel Surface Fire Spread Model

## Overview

The Rothermel surface fire spread model is a semi-empirical model that predicts the rate
of spread of surface fires in wildland fuels. Developed by Richard C. Rothermel in 1972,
it remains one of the most widely used fire behavior models in the United States.

The model calculates fire spread rate based on a heat source/heat sink balance:
- **Heat source**: Energy released by the fire (reaction intensity) modified by wind and slope
- **Heat sink**: Energy required to raise fuel ahead of the fire to ignition temperature

**Reference**: Andrews, Patricia L. 2018. The Rothermel surface fire spread model and associated developments: A comprehensive explanation. Gen. Tech. Rep. RMRS-GTR-371. Fort Collins, CO: U.S. Department of Agriculture, Forest Service, Rocky Mountain Research Station. 121 p.

```@docs
RothermelFireSpread
DynamicFuelLoadTransfer
LiveFuelMoistureExtinction
EffectiveMidflameWindSpeed
WindLimit
```

## Implementation

The Rothermel model is implemented as a ModelingToolkit.jl algebraic system. The main function
`RothermelFireSpread()` creates a system with 26 equations representing all intermediate
calculations and output variables.

### Important Note on Units

This implementation uses US customary units (ft, lb, Btu) internally, as the empirical
coefficients in the Rothermel equations were calibrated in this unit system. Key conversions:

| Quantity | US Customary | SI Equivalent |
|----------|--------------|---------------|
| Load | 1 ton/acre | 0.224 kg/m² |
| Load | 1 lb/ft² | 4.88 kg/m² |
| Wind speed | 1 mi/h | 26.82 m/min |
| Energy | 1 Btu | 1.055 kJ |
| Length | 1 ft | 0.305 m |
| Mass | 1 lb | 0.454 kg |

### State Variables

The model computes the following output variables:

| Variable | Description | Units |
|----------|-------------|-------|
| `R` | Rate of spread | ft/min |
| `R0` | No-wind no-slope rate of spread | ft/min |
| `IR` | Reaction intensity | Btu/ft²/min |
| `t_r` | Flame residence time | min |
| `HA` | Heat per unit area | Btu/ft² |
| `IB` | Fireline intensity (Byram) | Btu/ft/s |
| `F_L` | Flame length (Byram) | ft |

### Parameters

| Parameter | Description | Units |
|-----------|-------------|-------|
| `σ` | Surface-area-to-volume ratio | 1/ft |
| `w0` | Oven-dry fuel load | lb/ft² |
| `δ` | Fuel bed depth | ft |
| `Mx` | Dead fuel moisture of extinction | fraction |
| `Mf` | Fuel moisture content | fraction |
| `U` | Wind velocity at midflame height | ft/min |
| `tanϕ` | Slope steepness (rise/run) | dimensionless |
| `h` | Low heat content | Btu/lb |
| `S_T` | Total mineral content | fraction |
| `S_e` | Effective mineral content | fraction |
| `ρ_p` | Oven-dry particle density | lb/ft³ |

### Main Equations

The fundamental rate of spread equation is:

```math
R = \frac{I_R \cdot \xi \cdot (1 + \phi_w + \phi_s)}{\rho_b \cdot \varepsilon \cdot Q_{ig}}
```

Where:
- ``I_R`` = Reaction intensity
- ``\xi`` = Propagating flux ratio
- ``\phi_w`` = Wind factor
- ``\phi_s`` = Slope factor
- ``\rho_b`` = Bulk density
- ``\varepsilon`` = Effective heating number
- ``Q_{ig}`` = Heat of preignition

## Analysis

### Example: Fire Spread Under Different Conditions

```@example rothermel
using WildlandFire
using ModelingToolkit
using NonlinearSolve
using Plots

# Create the Rothermel system
sys = RothermelFireSpread()
compiled_sys = mtkcompile(sys)

# Fuel Model 1 parameters (short grass)
base_params = Dict(
    compiled_sys.σ => 3500.0,      # SAV ratio (1/ft)
    compiled_sys.w0 => 0.034,      # Fuel load (lb/ft²)
    compiled_sys.δ => 1.0,         # Fuel bed depth (ft)
    compiled_sys.Mx => 0.12,       # Moisture of extinction
    compiled_sys.Mf => 0.05,       # Fuel moisture content
    compiled_sys.tanϕ => 0.0       # Flat terrain
)

# Calculate spread rate for different wind speeds
wind_speeds = 0:50:500  # ft/min (0 to ~6 mi/h)
spread_rates = Float64[]

for U in wind_speeds
    prob = NonlinearProblem(compiled_sys, [], merge(base_params, Dict(compiled_sys.U => U)))
    sol = solve(prob)
    push!(spread_rates, sol[compiled_sys.R])
end

# Plot results
plot(wind_speeds ./ 88,  # Convert ft/min to mi/h
     spread_rates,
     xlabel = "Wind Speed (mi/h)",
     ylabel = "Rate of Spread (ft/min)",
     title = "Fire Spread Rate vs Wind Speed\n(Fuel Model 1, Short Grass)",
     legend = false,
     linewidth = 2,
     marker = :circle)
savefig("wind_effect.png"); nothing # hide
```

![Wind effect on spread rate](wind_effect.png)

### Effect of Fuel Moisture

```@example rothermel
# Calculate spread rate for different fuel moisture contents
moisture_contents = 0.02:0.01:0.12  # 2% to 12% (extinction)
spread_rates_moist = Float64[]

for Mf in moisture_contents
    params = merge(base_params, Dict(compiled_sys.Mf => Mf, compiled_sys.U => 220.0))  # 2.5 mi/h wind
    prob = NonlinearProblem(compiled_sys, [], params)
    sol = solve(prob)
    push!(spread_rates_moist, sol[compiled_sys.R])
end

plot(moisture_contents .* 100,
     spread_rates_moist,
     xlabel = "Fuel Moisture Content (%)",
     ylabel = "Rate of Spread (ft/min)",
     title = "Fire Spread Rate vs Fuel Moisture\n(Fuel Model 1, 2.5 mi/h wind)",
     legend = false,
     linewidth = 2,
     marker = :circle)
savefig("moisture_effect.png"); nothing # hide
```

![Moisture effect on spread rate](moisture_effect.png)

### Flame Length Prediction

```@example rothermel
# Calculate flame length for different wind speeds
flame_lengths = Float64[]

for U in wind_speeds
    prob = NonlinearProblem(compiled_sys, [], merge(base_params, Dict(compiled_sys.U => U)))
    sol = solve(prob)
    push!(flame_lengths, sol[compiled_sys.F_L])
end

plot(wind_speeds ./ 88,
     flame_lengths,
     xlabel = "Wind Speed (mi/h)",
     ylabel = "Flame Length (ft)",
     title = "Flame Length vs Wind Speed\n(Fuel Model 1, Short Grass)",
     legend = false,
     linewidth = 2,
     marker = :circle)
savefig("flame_length.png"); nothing # hide
```

![Flame length](flame_length.png)

## Standard Fuel Models

The Rothermel model uses standardized fuel model parameters. The 13 original fuel models and
40 Scott and Burgan fuel models define typical vegetation types. For Fuel Model 1 (short grass)
used in the examples above:

| Parameter | Value |
|-----------|-------|
| SAV ratio (σ) | 3,500 1/ft |
| Fuel load (w0) | 0.034 lb/ft² (0.74 ton/acre) |
| Fuel bed depth (δ) | 1.0 ft |
| Moisture of extinction (Mx) | 0.12 (12%) |
| Heat content (h) | 8,000 Btu/lb |

## Related Models

### Dynamic Fuel Model Load Transfer

The `DynamicFuelLoadTransfer` component simulates the curing of herbaceous fuels by
transferring load from live to dead fuel categories based on live herbaceous moisture:

```math
T = -1.11 \cdot M_{f,live} + 1.33 \quad (0 \leq T \leq 1)
```

### Live Fuel Moisture of Extinction

The `LiveFuelMoistureExtinction` component calculates the moisture content at which
live fuel will no longer sustain fire spread:

```math
M_{x,live} = \max\left(M_{x,dead}, 2.9 W \left(1 - \frac{M_{f,dead}}{M_{x,dead}}\right) - 0.226\right)
```

### Wind Limit

The `WindLimit` component calculates the maximum effective wind speed that can influence
fire spread rate. The corrected equation (Andrews et al. 2013) is:

```math
U_{limit} = 96.8 \cdot I_R^{1/3}
```
