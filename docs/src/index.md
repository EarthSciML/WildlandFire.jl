# WildlandFire.jl

WildlandFire.jl provides ModelingToolkit.jl-based equation systems for wildland fire behavior
modeling. The package is part of the [EarthSciML](https://github.com/EarthSciML) ecosystem.

## Features

- **Rothermel Surface Fire Spread Model**: Semi-empirical model for predicting fire spread rates
- **Dynamic Fuel Load Transfer**: Load transfer from live to dead herbaceous fuel based on curing
- **Live Fuel Moisture of Extinction**: Calculation of live fuel extinction moisture
- **Related Fire Behavior Models**: Flame length, fireline intensity, residence time

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/EarthSciML/WildlandFire.jl")
```

## Quick Start

```julia
using WildlandFire
using ModelingToolkit
using NonlinearSolve

# Create the Rothermel fire spread system
sys = RothermelFireSpread()
compiled_sys = mtkcompile(sys)

# Set parameters for short grass (Fuel Model 1) in SI units
# Original US values: σ=3500 1/ft, w0=0.034 lb/ft², δ=1.0 ft, U=440 ft/min (5 mi/h)
prob = NonlinearProblem(compiled_sys, [], [
    compiled_sys.σ => 11483.5,     # SAV ratio (1/m), converted from 3500 1/ft
    compiled_sys.w0 => 0.166,      # Fuel load (kg/m²), converted from 0.034 lb/ft²
    compiled_sys.δ => 0.3048,      # Fuel bed depth (m), converted from 1.0 ft
    compiled_sys.Mx => 0.12,       # Moisture of extinction (dimensionless)
    compiled_sys.Mf => 0.05,       # Fuel moisture content (dimensionless)
    compiled_sys.U => 2.235,       # Wind speed (m/s), converted from 5 mi/h
    compiled_sys.tanϕ => 0.0       # Flat terrain
])
sol = solve(prob)

# Get results in SI units
println("Rate of spread: ", sol[compiled_sys.R], " m/s")
println("Flame length: ", sol[compiled_sys.F_L], " m")
```

## Contents

```@contents
Pages = ["rothermel.md"]
Depth = 2
```

## Index

```@index
```
