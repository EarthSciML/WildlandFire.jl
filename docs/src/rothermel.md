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

This implementation uses SI units throughout. The original Rothermel equations were calibrated
in US customary units (ft, lb, Btu), so all empirical coefficients have been converted to SI.
Key conversions from the original formulation:

| Quantity | US Customary | SI Equivalent |
|----------|--------------|---------------|
| Load | 1 lb/ft² | 4.88 kg/m² |
| Load | 1 ton/acre | 0.224 kg/m² |
| Wind speed | 1 mi/h | 0.447 m/s |
| Wind speed | 1 ft/min | 0.00508 m/s |
| Energy | 1 Btu | 1.055 kJ |
| Length | 1 ft | 0.305 m |
| Mass | 1 lb | 0.454 kg |
| SAV ratio | 1 1/ft | 3.28 1/m |
| Density | 1 lb/ft³ | 16.02 kg/m³ |
| Heat content | 1 Btu/lb | 2326 J/kg |
| Power flux | 1 Btu/ft²/min | 189.3 W/m² |

### State Variables

```@example rothermel_impl
using DataFrames, ModelingToolkit, Symbolics, DynamicQuantities
using WildlandFire

sys = RothermelFireSpread()
vars = unknowns(sys)
units = [ModelingToolkit.get_unit(v) for v in vars]
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [isnothing(u) ? "dimensionless" : string(u) for u in units],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

### Parameters

```@example rothermel_impl
params = parameters(sys)
punits = [ModelingToolkit.get_unit(p) for p in params]
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params],
    :Units => [isnothing(u) ? "dimensionless" : string(u) for u in punits],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

### Equations

```@example rothermel_impl
equations(sys)
```

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

# Unit conversion factors
ft_to_m = 0.3048
lb_to_kg = 0.453592
lbft2_to_kgm2 = lb_to_kg / ft_to_m^2
invft_to_invm = 1 / ft_to_m
ftmin_to_ms = ft_to_m / 60.0

# Fuel Model 1 parameters (short grass) converted to SI
# Original US: σ = 3500 1/ft, w0 = 0.034 lb/ft², δ = 1.0 ft, Mx = 0.12
base_params = Dict(
    compiled_sys.σ => 3500.0 * invft_to_invm,      # SAV ratio (1/m)
    compiled_sys.w0 => 0.034 * lbft2_to_kgm2,      # Fuel load (kg/m²)
    compiled_sys.δ => 1.0 * ft_to_m,               # Fuel bed depth (m)
    compiled_sys.Mx => 0.12,                        # Moisture of extinction
    compiled_sys.Mf => 0.05,                        # Fuel moisture content
    compiled_sys.tanϕ => 0.0                        # Flat terrain
)

# Calculate spread rate for different wind speeds
# Wind speeds from 0 to ~6 mi/h in m/s
wind_speeds_ms = [U * ftmin_to_ms for U in 0:50:500]
spread_rates = Float64[]

for U in wind_speeds_ms
    prob = NonlinearProblem(compiled_sys, merge(base_params, Dict(compiled_sys.U => U)))
    sol = solve(prob)
    push!(spread_rates, sol[compiled_sys.R])
end

# Plot results (convert wind to mi/h and spread rate to ft/min for familiar units)
plot(wind_speeds_ms ./ 0.447,  # Convert m/s to mi/h
     spread_rates ./ ftmin_to_ms,  # Convert m/s to ft/min
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

U_wind = 220.0 * ftmin_to_ms  # 2.5 mi/h in m/s

for Mf in moisture_contents
    params = merge(base_params, Dict(compiled_sys.Mf => Mf, compiled_sys.U => U_wind))
    prob = NonlinearProblem(compiled_sys, params)
    sol = solve(prob)
    push!(spread_rates_moist, sol[compiled_sys.R])
end

plot(moisture_contents .* 100,
     spread_rates_moist ./ ftmin_to_ms,
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

for U in wind_speeds_ms
    prob = NonlinearProblem(compiled_sys, merge(base_params, Dict(compiled_sys.U => U)))
    sol = solve(prob)
    push!(flame_lengths, sol[compiled_sys.F_L])
end

plot(wind_speeds_ms ./ 0.447,  # Convert m/s to mi/h
     flame_lengths ./ ft_to_m,  # Convert m to ft
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

| Parameter | US Customary | SI Value |
|-----------|--------------|----------|
| SAV ratio (σ) | 3,500 1/ft | 11,483 1/m |
| Fuel load (w0) | 0.034 lb/ft² (0.74 ton/acre) | 0.166 kg/m² |
| Fuel bed depth (δ) | 1.0 ft | 0.305 m |
| Moisture of extinction (Mx) | 0.12 (12%) | 0.12 (12%) |
| Heat content (h) | 8,000 Btu/lb | 18,608,000 J/kg |
| Particle density (ρ_p) | 32 lb/ft³ | 512.6 kg/m³ |

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
fire spread rate. The corrected equation (Andrews et al. 2013) in SI units is:

```math
U_{limit} = 0.0857 \cdot I_R^{1/3}
```

where ``U_{limit}`` is in m/s and ``I_R`` is in W/m².
