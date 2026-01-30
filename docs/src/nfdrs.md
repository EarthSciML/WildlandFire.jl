# National Fire Danger Rating System (NFDRS)

## Overview

The National Fire Danger Rating System (NFDRS) is a comprehensive system for assessing fire danger across wildlands in the United States. This implementation provides ModelingToolkit.jl components for the NFDRS equations, enabling integration with other Earth system models.

**Reference**: Cohen, Jack D.; Deeming, John E. "The National Fire-Danger Rating System: basic equations." Gen. Tech. Rep. PSW-82. Berkeley, CA: Pacific Southwest Forest and Range Experiment Station, Forest Service, U.S. Department of Agriculture; 1985. 16 p.

## Unit Convention

This implementation uses **SI units** throughout for compatibility with other EarthSciML components:

| Quantity | SI Unit | Original NFDRS Unit |
|----------|---------|---------------------|
| Temperature | K (Kelvin) | °F (Fahrenheit) |
| Fuel loading | kg/m² | tons/acre, lb/ft² |
| Fuel bed depth | m | ft |
| SAV ratio | m⁻¹ | ft⁻¹ |
| Heat of combustion | J/kg | Btu/lb |
| Wind speed | m/s | mph |
| Rate of spread | m/s | ft/min |
| Moisture content | dimensionless (0-1) | percent |
| Density | kg/m³ | lb/ft³ |

The original NFDRS equations use imperial units. All empirical coefficients have been converted to SI units, with conversion factors explicitly defined in the source code.

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

### Spread Component

```@docs
SpreadComponent
```

### Energy Release Component

```@docs
EnergyReleaseComponent
```

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
fuel_loading_to_kg_per_sqm
```

### Available Fuel Models

The NFDRS includes 20 fuel models representing different vegetation types:

```@example fuel_models
using WildlandFire

for (code, model) in sort(collect(NFDRS_FUEL_MODELS), by=x->x[1])
    println("$code: $(model.description) (depth: $(round(model.DEPTH, digits=3)) m)")
end
```

### Fuel Model Parameters Table (Original Units)

The following table reproduces the fuel model parameters from the Appendix of Cohen & Deeming (1985), page 15, in their **original imperial units** for reference. The implementation converts these to SI units automatically.

| Model | Description | SG1 (ft⁻¹) | W1 (t/ac) | SG10 | W10 | Depth (ft) | MXD (%) | HD (Btu/lb) | SCM (ft/min) | WNDFC |
|-------|-------------|------------|-----------|------|-----|------------|---------|-------------|--------------|-------|
| A | Western grasses (annual) | 3000 | 0.20 | 109 | 0.0 | 0.80 | 15 | 8000 | 300 | 0.6 |
| B | California chaparral | 700 | 3.50 | 109 | 4.0 | 4.50 | 15 | 9500 | 58 | 0.5 |
| C | Pine-grass savanna | 2000 | 0.40 | 109 | 1.0 | 0.75 | 20 | 8000 | 32 | 0.4 |
| D | Southern rough | 1250 | 2.00 | 109 | 1.0 | 2.00 | 30 | 9000 | 25 | 0.4 |
| E | Hardwood litter (winter) | 2000 | 1.50 | 109 | 0.5 | 0.40 | 25 | 8000 | 25 | 0.4 |
| F | Intermediate brush | 700 | 2.50 | 109 | 2.0 | 4.50 | 15 | 9500 | 24 | 0.5 |
| G | Short needle (heavy dead) | 2000 | 2.50 | 109 | 2.0 | 1.00 | 25 | 8000 | 30 | 0.4 |
| H | Short needle (normal dead) | 2000 | 1.50 | 109 | 1.0 | 0.30 | 20 | 8000 | 8 | 0.4 |
| I | Heavy slash | 1500 | 12.0 | 109 | 12.0 | 2.00 | 25 | 8000 | 65 | 0.5 |
| J | Intermediate slash | 1500 | 7.00 | 109 | 7.0 | 1.30 | 25 | 8000 | 44 | 0.5 |
| K | Light slash | 1500 | 2.50 | 109 | 2.5 | 0.60 | 25 | 8000 | 23 | 0.5 |
| L | Western grasses (perennial) | 2000 | 0.25 | 109 | 1.5 | 1.00 | 15 | 8000 | 178 | 0.6 |
| N | Sawgrass | 1600 | 1.50 | 109 | 3.0 | 3.00 | 25 | 8700 | 167 | 0.6 |
| O | High pocosin | 1500 | 2.00 | 109 | 1.0 | 4.00 | 30 | 9000 | 99 | 0.5 |
| P | Southern pine plantation | 1750 | 1.00 | 109 | 2.5 | 0.40 | 30 | 8000 | 14 | 0.4 |
| Q | Alaskan black spruce | 1500 | 2.00 | 109 | 0.5 | 3.00 | 25 | 8000 | 59 | 0.4 |
| R | Hardwood litter (summer) | 1500 | 0.50 | 109 | 0.5 | 0.25 | 25 | 8000 | 6 | 0.4 |
| S | Tundra | 1500 | 0.50 | 109 | 0.5 | 0.40 | 25 | 8000 | 17 | 0.6 |
| T | Sagebrush-grass | 2500 | 1.00 | 109 | 1.5 | 1.25 | 15 | 8000 | 73 | 0.6 |
| U | Western pines | 1750 | 1.50 | 109 | 1.0 | 0.50 | 20 | 8000 | 16 | 0.4 |

Notes:
- SG1, SG10: Surface-area-to-volume ratio (ft⁻¹) → converted to m⁻¹ by multiplying by 3.28084
- W1, W10: Fuel loading (tons/acre) → converted to kg/m² by multiplying by 0.2241702
- MXD: Dead fuel moisture of extinction (percent) - unchanged
- HD: Heat of combustion (Btu/lb) → converted to J/kg by multiplying by 2326.0
- SCM: Spread component threshold (ft/min) → converted to m/s by multiplying by 0.00508
- WNDFC: Wind reduction factor (dimensionless) - unchanged

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
    ylabel = "Fuel Bed Depth (m)",
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

The following examples verify the implementation by running the actual ModelingToolkit components with known input values from Cohen & Deeming (1985).

#### EMC Equations (Eq. 1a, 1b, 1c)

The implementation uses SI units (temperature in Kelvin, RH as fraction 0-1).

```@example validation
using WildlandFire
using ModelingToolkit
using OrdinaryDiffEqDefault

# Create and compile the EMC system
emc_sys = EquilibriumMoistureContent()
compiled_emc = mtkcompile(emc_sys)

# Helper function to compute EMC for given temperature (°F) and RH (%)
function compute_emc(temp_f, rh_pct)
    # Convert °F to K: K = (°F + 459.67) × 5/9
    temp_k = (temp_f + 459.67) * 5 / 9
    rh_frac = rh_pct / 100.0

    prob = ODEProblem(compiled_emc,
        Dict(compiled_emc.TEMP => temp_k, compiled_emc.RH => rh_frac),
        (0.0, 1.0))
    sol = solve(prob)
    return sol[compiled_emc.EMC][end] * 100  # Convert fraction to %
end

# Test case 1: Low RH (Eq. 1a: RH < 10%)
# At RH=5%, TEMP=70°F: EMC ≈ 1.2%
println("EMC at RH=5%, T=70°F: $(round(compute_emc(70.0, 5.0), digits=2))%")

# Test case 2: Mid RH (Eq. 1b: 10% <= RH < 50%)
# At RH=30%, TEMP=70°F: EMC ≈ 6.0%
println("EMC at RH=30%, T=70°F: $(round(compute_emc(70.0, 30.0), digits=2))%")

# Test case 3: High RH (Eq. 1c: RH >= 50%)
# At RH=80%, TEMP=70°F: EMC ≈ 16.1%
println("EMC at RH=80%, T=70°F: $(round(compute_emc(70.0, 80.0), digits=2))%")
```

#### Burning Index (page 12)

```@example validation
using WildlandFire
using ModelingToolkit
using OrdinaryDiffEqDefault

# Create and compile the Burning Index system
bi_sys = BurningIndex()
compiled_bi = mtkcompile(bi_sys)

# SC=50, ERC=20 -> BI ≈ 60
prob = ODEProblem(compiled_bi,
    Dict(compiled_bi.SC => 50.0, compiled_bi.ERC => 20.0),
    (0.0, 1.0))
sol = solve(prob)
bi_result = sol[compiled_bi.BI][end]
println("BI at SC=50, ERC=20: $(round(bi_result, digits=1))")
```

#### Fire Load Index (page 14)

```@example validation
using WildlandFire
using ModelingToolkit
using OrdinaryDiffEqDefault

# Create and compile the Fire Load Index system
fli_sys = FireLoadIndex()
compiled_fli = mtkcompile(fli_sys)

# BI=50, LOI=30, MCOI=20
prob = ODEProblem(compiled_fli,
    Dict(compiled_fli.BI => 50.0, compiled_fli.LOI => 30.0, compiled_fli.MCOI => 20.0),
    (0.0, 1.0))
sol = solve(prob)
fli_result = sol[compiled_fli.FLI][end]
println("FLI at BI=50, LOI=30, MCOI=20: $(round(fli_result, digits=1))")
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
    ylabel = "Dead Fuel Loading (kg/m²)",
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

bar([string(m.name) for m in models], [m.HD / 1e6 for m in models],
    xlabel = "Fuel Model",
    ylabel = "Heat of Combustion (MJ/kg)",
    title = "NFDRS Dead Fuel Heat Content",
    legend = false,
    rotation = 45,
    ylims = (17, 23))
savefig("heat_content.svg"); nothing # hide
```

![Heat Content](heat_content.svg)
