# WildlandFire.jl

WildlandFire.jl provides ModelingToolkit.jl-based implementations of wildland fire behavior models for use in Earth science simulations.

## Overview

This package implements the Rothermel surface fire spread model and associated developments as described in:

> Andrews, Patricia L. 2018. The Rothermel surface fire spread model and associated developments: A comprehensive explanation. Gen. Tech. Rep. RMRS-GTR-371. Fort Collins, CO: U.S. Department of Agriculture, Forest Service, Rocky Mountain Research Station. 121 p.

## Features

- **Rothermel Fire Spread Model**: Core fire spread rate calculations including wind and slope effects
- **Fire Spread Direction**: Vector addition for calculating direction of maximum spread when wind is not aligned with slope
- **Elliptical Fire Spread**: Fire shape calculations from a single ignition point
- **Fire Perimeter Spread**: Rate of spread normal to the fire perimeter

## Quick Start

```julia
using WildlandFire, ModelingToolkit, OrdinaryDiffEqDefault

# Create a fire spread direction system
sys = FireSpreadDirection()
compiled_sys = mtkcompile(sys)

# Set up problem with cross-slope wind
prob = NonlinearProblem(compiled_sys, Dict(
    compiled_sys.R0 => 0.01,      # No-wind no-slope rate of spread (m/s)
    compiled_sys.φw => 5.0,        # Wind factor
    compiled_sys.φs => 2.0,        # Slope factor
    compiled_sys.ω => π/4,         # Wind direction 45° from upslope
    compiled_sys.β_ratio => 0.5,
    compiled_sys.C_coeff => 7.47,
    compiled_sys.B_coeff => 0.5,
    compiled_sys.E_coeff => 0.5,
    compiled_sys.elapsed_time => 60.0
))

sol = solve(prob)
```

## Documentation Contents

```@contents
Pages = ["fire_spread_direction.md"]
```
