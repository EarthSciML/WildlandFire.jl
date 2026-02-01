# WildlandFire.jl

```@meta
CurrentModule = WildlandFire
```

WildlandFire.jl provides ModelingToolkit.jl-based equation systems for wildland fire behavior
modeling. The package is part of the [EarthSciML](https://github.com/EarthSciML) ecosystem.

## Features

- **Rothermel Surface Fire Spread Model**: Semi-empirical model for predicting fire spread rates
- **Dynamic Fuel Load Transfer**: Load transfer from live to dead herbaceous fuel based on curing
- **Live Fuel Moisture of Extinction**: Calculation of live fuel extinction moisture
- **Related Fire Behavior Models**: Flame length, fireline intensity, residence time
- **Fire Spread Direction**: Vector addition for calculating direction of maximum spread when wind is not aligned with slope
- **Elliptical Fire Spread**: Fire shape calculations from a single ignition point
- **Fire Perimeter Spread**: Rate of spread normal to the fire perimeter
- **National Fire Danger Rating System (NFDRS)**: Fuel moisture and fire danger assessment
  - Fuel Moisture Models: Dead fuel moisture (1-hr, 10-hr, 100-hr, 1000-hr timelag classes) and live fuel moisture (herbaceous and woody)
  - Fire Behavior Indices: Spread Component, Energy Release Component, Burning Index
  - Fire Occurrence Indices: Ignition Component, Human-Caused Fire Occurrence Index, Fire Load Index
  - Fuel Model Database: All 20 NFDRS fuel models (A-U, excluding M) with vegetation-specific parameters

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/EarthSciML/WildlandFire.jl")
```

## Quick Start

### Rothermel Fire Spread

```julia
using WildlandFire
using ModelingToolkit
using NonlinearSolve

# Create the Rothermel fire spread system
sys = RothermelFireSpread()
compiled_sys = mtkcompile(sys)

# Set parameters for short grass (Fuel Model 1) in SI units
# Original US values: σ=3500 1/ft, w0=0.034 lb/ft², δ=1.0 ft, U=440 ft/min (5 mi/h)
prob = NonlinearProblem(compiled_sys, Dict(
    compiled_sys.σ => 11483.5,     # SAV ratio (1/m), converted from 3500 1/ft
    compiled_sys.w0 => 0.166,      # Fuel load (kg/m²), converted from 0.034 lb/ft²
    compiled_sys.δ => 0.3048,      # Fuel bed depth (m), converted from 1.0 ft
    compiled_sys.Mx => 0.12,       # Moisture of extinction (dimensionless)
    compiled_sys.Mf => 0.05,       # Fuel moisture content (dimensionless)
    compiled_sys.U => 2.235,       # Wind speed (m/s), converted from 5 mi/h
    compiled_sys.tanϕ => 0.0       # Flat terrain
))
sol = solve(prob)

# Get results in SI units
println("Rate of spread: ", sol[compiled_sys.R], " m/s")
println("Flame length: ", sol[compiled_sys.F_L], " m")
```

### NFDRS Fire Danger Rating

```julia
using WildlandFire
using ModelingToolkit

# Create an equilibrium moisture content model
emc = EquilibriumMoistureContent()

# Get fuel model parameters
fuel_model_a = get_fuel_model(:A)  # Western grasses (annual)

# Create a burning index model
bi = BurningIndex()
```

## Contents

```@contents
Pages = ["nfdrs.md", "rothermel.md", "fire_spread_direction.md"]
Depth = 2
```

## Module Documentation

```@docs
WildlandFire.WildlandFire
```

## Index

```@index
```
