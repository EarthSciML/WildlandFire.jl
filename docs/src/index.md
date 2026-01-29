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

# Set parameters for short grass (Fuel Model 1)
prob = NonlinearProblem(compiled_sys, [], [
    compiled_sys.σ => 3500.0,      # SAV ratio (1/ft)
    compiled_sys.w0 => 0.034,      # Fuel load (lb/ft²)
    compiled_sys.δ => 1.0,         # Fuel bed depth (ft)
    compiled_sys.Mx => 0.12,       # Moisture of extinction
    compiled_sys.Mf => 0.05,       # Fuel moisture content
    compiled_sys.U => 440.0,       # Wind speed (5 mi/h)
    compiled_sys.tanϕ => 0.0       # Flat terrain
])
sol = solve(prob)

# Get results
println("Rate of spread: ", sol[compiled_sys.R], " ft/min")
println("Flame length: ", sol[compiled_sys.F_L], " ft")
```

## Contents

```@contents
Pages = ["rothermel.md"]
Depth = 2
```

## Index

```@index
```
