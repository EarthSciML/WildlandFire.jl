# WildlandFire.jl

*A Julia package for wildland fire behavior modeling using ModelingToolkit.jl*

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://wildlandfire.earthsci.dev)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://wildlandfire.earthsci.dev/dev)

## Overview

WildlandFire.jl provides a complete, validated implementation of the Rothermel surface fire spread model based on the foundational work of Rothermel (1972) and Albini (1976a). The package leverages [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl) for symbolic modeling, enabling automatic differentiation, optimization, and seamless integration with the broader SciML ecosystem.

## Key Features

- **Complete Rothermel Model**: All equations from Rothermel (1972) and Albini (1976a)
- **Validated Physics**: Thoroughly tested against expected fire behavior
- **Homogeneous & Heterogeneous Fuels**: Support for both single-class and multi-class fuel models
- **Symbolic Framework**: Built on ModelingToolkit.jl for maximum flexibility
- **Production Ready**: Comprehensive test suite and validation
- **Well Documented**: Extensive examples and API documentation

## What is the Rothermel Model?

The Rothermel fire spread model is a semi-empirical mathematical model that predicts the forward rate of spread of a surface fire in wildland fuels. It was developed by Richard C. Rothermel in 1972 and has become the foundation for fire behavior prediction systems used worldwide, including:

- BEHAVE fire behavior prediction system
- BehavePlus fire modeling system
- FlamMap fire behavior mapping system
- FARSITE fire area simulator

The model calculates how fast a fire will spread based on:
- **Fuel properties**: Type, load, moisture, and structure
- **Environmental conditions**: Wind speed and slope
- **Fire physics**: Heat transfer, combustion rates, and energy balance

## Quick Example

```@example index
using WildlandFire
using OrdinaryDiffEqDefault

# Create the Rothermel model
sys = Rothermel()

# Define fuel and environmental parameters for grass
params = [
    h => 8000.0,          # Heat content (Btu/lb)
    σ => 3500.0,          # Surface-area-to-volume ratio (ft²/ft³)
    w_o => 0.138,         # Fuel load (lb/ft²)
    δ => 1.0,             # Fuel bed depth (ft)
    M_f => 0.05,          # Fuel moisture (5%)
    U => 352.0,           # Wind speed (ft/min) ≈ 4 mph
    tan_ϕ => 0.0,         # Slope (flat)
    S_T => 0.0555,        # Total mineral content
    S_e => 0.010,         # Effective mineral content
    ρ_p => 32.0,          # Particle density (lb/ft³)
    M_x => 0.12,          # Moisture of extinction
]

# Solve for rate of spread
prob = ODEProblem(sys, [], (0.0, 1.0), params)
sol = solve(prob)

# Get the rate of spread
R_ftmin = sol[R][end]
println("Rate of spread: $(round(R_ftmin, digits=2)) ft/min")
```

## Installation

```julia
using Pkg
Pkg.add("WildlandFire")
```

## Package Contents

```@contents
Pages = [
    "getting_started.md",
    "model_description.md",
    "usage_guide.md",
    "examples/basic.md",
    "examples/sensitivity.md",
    "examples/fuel_types.md",
    "api.md",
    "references.md",
]
Depth = 1
```

## Citation

If you use this package in research, please cite both the package and the original Rothermel paper:

**Package Citation:**
```
EarthSciML authors and contributors (2026). WildlandFire.jl: Wildland Fire
Behavior Modeling in Julia. https://github.com/EarthSciML/WildlandFire.jl
```

**Original Model Citation:**
```
Rothermel, R.C. (1972). A mathematical model for predicting fire spread in
wildland fuels. USDA Forest Service Research Paper INT-115.
```

## License

MIT License - Copyright (c) 2026 EarthSciML authors and contributors

## Acknowledgments

This implementation is based on the foundational work of:
- **Richard C. Rothermel** - Original fire spread model (1972)
- **Frank A. Albini** - Model extensions and refinements (1976)
- **Patricia L. Andrews** - Comprehensive model documentation (2018)

## Getting Help

- **Documentation**: [https://wildlandfire.earthsci.dev](https://wildlandfire.earthsci.dev)
- **Issues**: [GitHub Issues](https://github.com/EarthSciML/WildlandFire.jl/issues)
- **Discussions**: [GitHub Discussions](https://github.com/EarthSciML/WildlandFire.jl/discussions)
