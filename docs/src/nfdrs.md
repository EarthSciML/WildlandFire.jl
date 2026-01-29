# National Fire Danger Rating System (NFDRS)

## Overview

The National Fire Danger Rating System (NFDRS) is a comprehensive system for assessing fire danger across wildlands in the United States. This implementation provides ModelingToolkit.jl components for the NFDRS equations, enabling integration with other Earth system models.

**Reference**: Cohen, Jack D.; Deeming, John E. "The National Fire-Danger Rating System: basic equations." Gen. Tech. Rep. PSW-82. Berkeley, CA: Pacific Southwest Forest and Range Experiment Station, Forest Service, U.S. Department of Agriculture; 1985. 16 p.

Note: The NFDRS equations use imperial units (°F, lb, ft, Btu, etc.) as originally published. All inputs and outputs are in these original units unless otherwise noted.

## Equilibrium Moisture Content

```@docs
EquilibriumMoistureContent
```

### Equations

The EMC is calculated using regression equations (Eq. 1a, 1b, 1c from Cohen & Deeming 1985) based on relative humidity ranges:
- RH < 10%: Linear approximation
- 10% ≤ RH < 50%: Intermediate range equation
- RH ≥ 50%: High humidity equation

## Dead Fuel Moisture Models

### 1-Hour Timelag Fuel Moisture

```@docs
OneHourFuelMoisture
```

### 10-Hour Timelag Fuel Moisture

```@docs
TenHourFuelMoisture
```

### 100-Hour Timelag Fuel Moisture

```@docs
HundredHourFuelMoisture
```

### 1000-Hour Timelag Fuel Moisture

```@docs
ThousandHourFuelMoisture
```

## Live Fuel Moisture Models

### Herbaceous Fuel Moisture

```@docs
HerbaceousFuelMoisture
```

### Woody Fuel Moisture

```@docs
WoodyFuelMoisture
```

## Fuel Loading Transfer

```@docs
FuelLoadingTransfer
```

## Fire Behavior Indices

### Burning Index

```@docs
BurningIndex
```

### Ignition Component

```@docs
IgnitionComponent
```

## Fire Occurrence Indices

### Human-Caused Fire Occurrence Index

```@docs
HumanFireOccurrenceIndex
```

### Fire Load Index

```@docs
FireLoadIndex
```

## Fuel Models

```@docs
NFDRSFuelModel
NFDRS_FUEL_MODELS
get_fuel_model
```

### Available Fuel Models

The NFDRS includes 20 fuel models representing different vegetation types:

```@example fuel_models
using WildlandFire

for (code, model) in sort(collect(NFDRS_FUEL_MODELS), by=x->x[1])
    println("$code: $(model.description) (depth: $(model.DEPTH) ft)")
end
```

## Analysis

### Fuel Model Comparison

Different fuel models have varying fuel bed depths and loadings, affecting fire behavior:

```@example fuel_comparison
using WildlandFire
using Plots

models = collect(values(NFDRS_FUEL_MODELS))
sort!(models, by=m -> m.DEPTH)

bar([string(m.name) for m in models], [m.DEPTH for m in models],
    xlabel = "Fuel Model",
    ylabel = "Fuel Bed Depth (ft)",
    title = "NFDRS Fuel Model Bed Depths",
    legend = false,
    rotation = 45)
savefig("fuel_depths.svg"); nothing # hide
```

![Fuel Bed Depths](fuel_depths.svg)
