# National Fire Danger Rating System (NFDRS)

## Overview

The National Fire Danger Rating System (NFDRS) is a comprehensive system for assessing fire danger across wildlands in the United States. This implementation provides ModelingToolkit.jl components for the NFDRS equations, enabling integration with other Earth system models.

**Reference**: Cohen, Jack D.; Deeming, John E. "The National Fire-Danger Rating System: basic equations." Gen. Tech. Rep. PSW-82. Berkeley, CA: Pacific Southwest Forest and Range Experiment Station, Forest Service, U.S. Department of Agriculture; 1985. 16 p.

```@docs
EquilibriumMoistureContent
```

## Implementation

This section demonstrates the structure of the NFDRS ModelingToolkit.jl components.

### Unit Convention

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

### Equilibrium Moisture Content System

The EMC component calculates equilibrium moisture content using regression equations (Eq. 1a, 1b, 1c from Cohen & Deeming 1985):

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
    :Units => [string(dimension(ModelingToolkit.get_unit(p))) for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

```@example implementation
# Variables
vars = unknowns(emc)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [string(dimension(ModelingToolkit.get_unit(v))) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

### Equations

The EMC equations (Eq. 1a, 1b, 1c from Cohen & Deeming 1985):

```@example implementation
equations(emc)
```

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

The Spread Component is the most complex subsystem, implementing Rothermel's fire spread model:

```@docs
SpreadComponent
```

```@example implementation
sc = SpreadComponent()

# Parameters with units
params = parameters(sc)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params],
    :Units => [string(dimension(ModelingToolkit.get_unit(p))) for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

### Energy Release Component

```@docs
EnergyReleaseComponent
```

### Burning Index

```@docs
BurningIndex
```

The Burning Index equation (page 12):

```@example implementation
bi = BurningIndex()
equations(bi)
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

The Fire Load Index equation (page 14):

```@example implementation
fli = FireLoadIndex()
equations(fli)
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

This section validates the implementation against expected values from Cohen & Deeming (1985) and demonstrates the behavior of the NFDRS components.

### EMC Validation (Eq. 1a, 1b, 1c, page 1)

The equilibrium moisture content equations calculate EMC based on temperature and relative humidity. The implementation uses SI units (temperature in Kelvin, RH as fraction 0-1).

```@example validation
using WildlandFire
using ModelingToolkit
using OrdinaryDiffEqDefault
using Plots

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

# Test cases from Cohen & Deeming (1985)
println("EMC Validation at T=70°F:")
println("  RH=5% (Eq. 1a):  EMC = $(round(compute_emc(70.0, 5.0), digits=2))%")
println("  RH=30% (Eq. 1b): EMC = $(round(compute_emc(70.0, 30.0), digits=2))%")
println("  RH=80% (Eq. 1c): EMC = $(round(compute_emc(70.0, 80.0), digits=2))%")
```

The following plot shows EMC as a function of relative humidity at 70°F, demonstrating the three humidity regimes:

```@example validation
rh_range = 1:99
emc_values = [compute_emc(70.0, rh) for rh in rh_range]

p = plot(rh_range, emc_values,
    xlabel = "Relative Humidity (%)",
    ylabel = "EMC (%)",
    title = "Equilibrium Moisture Content at 70°F",
    legend = false,
    linewidth = 2)

# Mark the regime boundaries
vline!(p, [10, 50], linestyle=:dash, color=:gray, label="")
annotate!(p, [(5, 15, text("Eq. 1a", 8)),
              (30, 15, text("Eq. 1b", 8)),
              (75, 15, text("Eq. 1c", 8))])
p
```

### EMC Temperature Sensitivity

EMC decreases with increasing temperature at constant humidity:

```@example validation
temps_f = 40:10:100
rh_values = [30, 50, 70]

p = plot(xlabel = "Temperature (°F)",
         ylabel = "EMC (%)",
         title = "EMC Temperature Sensitivity",
         legend = :topright)

for rh in rh_values
    emc_vals = [compute_emc(t, rh) for t in temps_f]
    plot!(p, temps_f, emc_vals, label="RH=$rh%", linewidth=2)
end
p
```

### Burning Index Validation (page 12)

The Burning Index (BI) is calculated from the Spread Component and Energy Release Component:

BI = 3.01 × (SC × ERC)^0.46

where SC is in ft/min for the original equation.

```@example validation
using WildlandFire
using ModelingToolkit
using OrdinaryDiffEqDefault
using Plots

bi_sys = BurningIndex()
compiled_bi = mtkcompile(bi_sys)

# Conversion factor: m/s to ft/min
FPM_TO_MS = 0.00508

# Calculate BI for various SC and ERC combinations
function compute_bi(sc_fpm, erc)
    sc_ms = sc_fpm * FPM_TO_MS
    prob = ODEProblem(compiled_bi,
        Dict(compiled_bi.SC => sc_ms, compiled_bi.ERC => erc, compiled_bi.fuels_wet => 0.0),
        (0.0, 1.0))
    sol = solve(prob)
    return sol[compiled_bi.BI][end]
end

# Validate at SC=50 ft/min, ERC=20
expected_bi = 3.01 * (50.0 * 20.0)^0.46
computed_bi = compute_bi(50.0, 20.0)
println("BI Validation:")
println("  SC=50 ft/min, ERC=20")
println("  Expected: $(round(expected_bi, digits=1))")
println("  Computed: $(round(computed_bi, digits=1))")
```

BI contour plot showing the relationship between SC, ERC, and fire behavior intensity:

```@example validation
sc_range = 10:10:100
erc_range = 5:5:50
bi_matrix = [compute_bi(sc, erc) for erc in erc_range, sc in sc_range]

contourf(sc_range, erc_range, bi_matrix,
    xlabel = "Spread Component (ft/min)",
    ylabel = "Energy Release Component",
    title = "Burning Index",
    colorbar_title = "BI",
    levels = 10)
```

### Fire Load Index Validation (page 14)

The Fire Load Index combines the Burning Index with fire occurrence indices:

FLI = 0.71 × √(BI² + (LOI + MCOI)²)

```@example validation
using WildlandFire
using ModelingToolkit
using OrdinaryDiffEqDefault

fli_sys = FireLoadIndex()
compiled_fli = mtkcompile(fli_sys)

# Test case from Cohen & Deeming (1985)
prob = ODEProblem(compiled_fli,
    Dict(compiled_fli.BI => 50.0, compiled_fli.LOI => 30.0, compiled_fli.MCOI => 20.0),
    (0.0, 1.0))
sol = solve(prob)
computed_fli = sol[compiled_fli.FLI][end]
expected_fli = 0.71 * sqrt(50.0^2 + (30.0 + 20.0)^2)

println("FLI Validation:")
println("  BI=50, LOI=30, MCOI=20")
println("  Expected: $(round(expected_fli, digits=1))")
println("  Computed: $(round(computed_fli, digits=1))")
```

### Fuel Model Comparison

Different fuel models have varying fuel bed depths and loadings, affecting fire behavior:

```@example fuel_comparison
using WildlandFire
using Plots

models = collect(values(NFDRS_FUEL_MODELS))
sort!(models, by=m -> m.DEPTH)

p = bar([string(m.name) for m in models], [m.DEPTH for m in models],
    xlabel = "Fuel Model",
    ylabel = "Fuel Bed Depth (m)",
    title = "NFDRS Fuel Model Bed Depths (Cohen & Deeming 1985, Appendix)",
    legend = false,
    rotation = 45)
p
```

### Dead Fuel Loading Comparison

Comparison of total dead fuel loadings (1-hr + 10-hr + 100-hr) across different NFDRS fuel models:

```@example fuel_analysis
using WildlandFire
using Plots

models = collect(values(NFDRS_FUEL_MODELS))
sort!(models, by=m -> m.W1 + m.W10 + m.W100)

# Total dead fuel loading (1-hr + 10-hr + 100-hr)
dead_loading = [m.W1 + m.W10 + m.W100 for m in models]

p = bar([string(m.name) for m in models], dead_loading,
    xlabel = "Fuel Model",
    ylabel = "Dead Fuel Loading (kg/m²)",
    title = "NFDRS Dead Fuel Loading by Model",
    legend = false,
    rotation = 45)
p
```

### Heat Content Comparison

Different fuel models have varying heat content values reflecting their chemical composition:

```@example fuel_analysis
models = collect(values(NFDRS_FUEL_MODELS))
sort!(models, by=m -> m.HD)

p = bar([string(m.name) for m in models], [m.HD / 1e6 for m in models],
    xlabel = "Fuel Model",
    ylabel = "Heat of Combustion (MJ/kg)",
    title = "NFDRS Dead Fuel Heat Content",
    legend = false,
    rotation = 45,
    ylims = (17, 23))
p
```

### Fuel Loading Transfer

The fuel loading transfer model (Eq. 5-8) calculates the fraction of herbaceous fuel that behaves as dead fuel based on the herbaceous moisture content:

```@example fuel_transfer
using WildlandFire
using ModelingToolkit
using OrdinaryDiffEqDefault
using Plots

flt_sys = FuelLoadingTransfer()
compiled_flt = mtkcompile(flt_sys)

# Calculate FCTCUR across the moisture range
function compute_fctcur(mcherb_frac)
    prob = ODEProblem(compiled_flt,
        Dict(compiled_flt.MCHERB => mcherb_frac,
             compiled_flt.W1 => 0.05,
             compiled_flt.WHERB => 0.02),
        (0.0, 1.0))
    sol = solve(prob)
    return sol[compiled_flt.FCTCUR][end]
end

mcherb_range = 0.30:0.01:1.30
fctcur_values = [compute_fctcur(mc) for mc in mcherb_range]

plot(mcherb_range .* 100, fctcur_values .* 100,
    xlabel = "Herbaceous Moisture Content (%)",
    ylabel = "Fraction Transferred to Dead (%)",
    title = "Fuel Loading Transfer (Eq. 5)",
    legend = false,
    linewidth = 2)
```
