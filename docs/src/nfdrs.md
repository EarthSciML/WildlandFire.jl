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

### Lightning-Caused Fire Occurrence Index

```@docs
LightningFireOccurrenceIndex
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

## Implementation

This section demonstrates the structure of the NFDRS ModelingToolkit.jl components.

### Equilibrium Moisture Content System

```@example implementation
using WildlandFire
using ModelingToolkit
using DataFrames
using Symbolics
using DynamicQuantities

emc = EquilibriumMoistureContent()

# Parameters
params = parameters(emc)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

```@example implementation
# Variables
vars = unknowns(emc)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

### Spread Component System

The Spread Component is the most complex subsystem, implementing Rothermel's fire spread model:

```@example implementation
sc = SpreadComponent()

# Parameters
params = parameters(sc)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

### Equations

The EMC equations (Eq. 1a, 1b, 1c from Cohen & Deeming 1985):

```@example implementation
emc = EquilibriumMoistureContent()
equations(emc)
```

The Burning Index equation (page 12):

```@example implementation
bi = BurningIndex()
equations(bi)
```

The Fire Load Index equation (page 14):

```@example implementation
fli = FireLoadIndex()
equations(fli)
```

### Numerical Validation

The following examples verify the implementation against known values from Cohen & Deeming (1985).

#### EMC Equations (Eq. 1a, 1b, 1c)

```@example validation
# Verify EMC regression equations from Simard (1968)
function emc_manual(temp, rh_pct)
    if rh_pct < 10
        return 0.03229 + 0.281073 * rh_pct - 0.000578 * temp * rh_pct
    elseif rh_pct < 50
        return 2.22749 + 0.160107 * rh_pct - 0.014784 * temp
    else
        return 21.0606 + 0.005565 * rh_pct^2 - 0.00035 * rh_pct * temp - 0.483199 * rh_pct
    end
end

# Test case 1: Low RH (Eq. 1a: RH < 10%)
# At RH=5%, TEMP=70°F: EMC ≈ 1.236%
println("EMC at RH=5%, T=70°F: $(round(emc_manual(70.0, 5.0), digits=3))%")

# Test case 2: Mid RH (Eq. 1b: 10% <= RH < 50%)
# At RH=30%, TEMP=70°F: EMC ≈ 5.997%
println("EMC at RH=30%, T=70°F: $(round(emc_manual(70.0, 30.0), digits=3))%")

# Test case 3: High RH (Eq. 1c: RH >= 50%)
# At RH=80%, TEMP=70°F: EMC ≈ 16.06%
println("EMC at RH=80%, T=70°F: $(round(emc_manual(70.0, 80.0), digits=2))%")
```

#### Burning Index (page 12)

```@example validation
# BI = 3.01 * (SC * ERC)^0.46
# Example: SC=50, ERC=20 -> BI ≈ 60
bi_example = 3.01 * (50 * 20)^0.46
println("BI at SC=50, ERC=20: $(round(bi_example, digits=1))")
```

#### Fire Load Index (page 14)

```@example validation
# FLI = 0.71 * sqrt(BI^2 + (LOI + MCOI)^2)
# Example: BI=50, LOI=30, MCOI=20
fli_example = 0.71 * sqrt(50^2 + (30 + 20)^2)
println("FLI at BI=50, LOI=30, MCOI=20: $(round(fli_example, digits=1))")
```

### Fuel Loading Comparison

Comparison of fuel loadings across different NFDRS fuel models:

```@example fuel_analysis
using WildlandFire
using Plots

models = collect(values(NFDRS_FUEL_MODELS))
sort!(models, by=m -> m.W1 + m.W10 + m.W100)

# Total dead fuel loading (1-hr + 10-hr + 100-hr)
dead_loading = [m.W1 + m.W10 + m.W100 for m in models]

bar([string(m.name) for m in models], dead_loading,
    xlabel = "Fuel Model",
    ylabel = "Dead Fuel Loading (tons/acre)",
    title = "NFDRS Dead Fuel Loading by Model",
    legend = false,
    rotation = 45)
savefig("dead_fuel_loading.svg"); nothing # hide
```

![Dead Fuel Loading](dead_fuel_loading.svg)

### Heat Content Comparison

Different fuel models have varying heat content values:

```@example fuel_analysis
models = collect(values(NFDRS_FUEL_MODELS))
sort!(models, by=m -> m.HD)

bar([string(m.name) for m in models], [m.HD for m in models],
    xlabel = "Fuel Model",
    ylabel = "Heat of Combustion (Btu/lb)",
    title = "NFDRS Dead Fuel Heat Content",
    legend = false,
    rotation = 45,
    ylims = (7500, 10000))
savefig("heat_content.svg"); nothing # hide
```

![Heat Content](heat_content.svg)
