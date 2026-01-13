# WildlandFire.jl

[![MIT License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

A Julia package for wildland fire behavior modeling using [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl).

## Overview

WildlandFire.jl provides a complete, validated implementation of the Rothermel surface fire spread model based on Rothermel (1972) and Albini (1976a). The package uses symbolic modeling with ModelingToolkit.jl for automatic differentiation, optimization, and seamless integration with the SciML ecosystem.

### Key Features

- ✅ **Complete Rothermel Model**: All equations from the original papers
- ✅ **Validated Physics**: Tested against expected fire behavior
- ✅ **Homogeneous & Heterogeneous Fuels**: Both single-class and multi-class fuel models
- ✅ **Symbolic Framework**: Built on ModelingToolkit.jl for maximum flexibility
- ✅ **Well-Documented**: Comprehensive documentation and examples
- ✅ **Production-Ready**: Validated results and test suite included

## Installation

```julia
using Pkg
Pkg.add("WildlandFire")
```

Or for development:

```julia
using Pkg
Pkg.develop(path="/path/to/WildlandFire.jl")
```

## Quick Start

```julia
using WildlandFire
using OrdinaryDiffEqDefault

# Create the Rothermel model using the factory function
sys = Rothermel()

# Define fuel and environmental parameters
params = [
    h => 8000.0,          # Heat content (Btu/lb)
    S_T => 0.0555,        # Total mineral content
    S_e => 0.010,         # Effective mineral content
    ρ_p => 32.0,          # Particle density (lb/ft³)
    σ => 3500.0,          # Surface-area-to-volume ratio (ft²/ft³)
    w_o => 0.138,         # Fuel load (lb/ft²)
    δ => 1.0,             # Fuel bed depth (ft)
    M_x => 0.12,          # Moisture of extinction
    M_f => 0.05,          # Fuel moisture
    U => 352.0,           # Wind speed (ft/min) ≈ 4 mph
    tan_ϕ => 0.0,         # Slope (fraction)
]

# Solve for steady-state rate of spread
prob = ODEProblem(sys, [], (0.0, 1.0), params)
sol = solve(prob)

# Extract rate of spread
R_ftmin = sol[R][end]      # ft/min
R_mhr = R_ftmin * 18.288   # m/hr
R_chains = R_ftmin * 1.829 # chains/hr

println("Rate of spread: $R_ftmin ft/min")
println("Rate of spread: $R_mhr m/hr")
println("Rate of spread: $R_chains chains/hr")
```

**Alternative:** You can also use the pre-built model for backward compatibility:
```julia
sys = rothermel_simplified  # Pre-built simplified model
```

## Model Description

The Rothermel model calculates the forward rate of spread (R) of a surface fire using:

```
R = I_R·ξ(1 + ϕ_w + ϕ_s) / (ρ_b·ε·Q_ig)
```

Where:
- **I_R**: Reaction intensity (Btu/ft²·min) - heat release rate
- **ξ**: Propagating flux ratio - proportion heating adjacent fuel
- **ϕ_w**: Wind factor - dimensionless wind effect multiplier
- **ϕ_s**: Slope factor - dimensionless slope effect multiplier
- **ρ_b**: Bulk density (lb/ft³)
- **ε**: Effective heating number - particle size effect
- **Q_ig**: Heat of preignition (Btu/lb) - energy to ignite fuel

### Key Components

#### Reaction Intensity (I_R)
Rate of heat release from burning fuel:
- Depends on fuel load, heat content, and damping coefficients
- Typical values: 1000-10000 Btu/ft²·min

#### Wind Factor (ϕ_w)
Exponential increase in spread with wind:
- Formula: C × U^B × (β/β_opt)^(-E)
- Strong effect at higher wind speeds
- Example: 10 mph wind can increase ROS by 50-100x

#### Slope Factor (ϕ_s)
Quadratic increase with slope angle:
- Formula: 5.275β^(-0.3)(tan ϕ)²
- Upslope fires spread much faster
- Example: 30% slope increases ROS by ~3-4x

#### Moisture Damping (η_M)
Reduction in reaction intensity due to fuel moisture:
- Cubic polynomial function approaching zero at extinction moisture
- Wet fuels require more energy to ignite

## Examples

See the [`examples/`](examples/) directory for detailed examples:

- **[basic_examples.jl](examples/basic_examples.jl)**: Comprehensive demonstrations including:
  - Basic fire spread calculation
  - Wind speed sensitivity analysis
  - Fuel moisture sensitivity analysis
  - Slope effect analysis
  - Comparison of different NFFL fuel types

- **[wind_factor_diagnostic.jl](examples/wind_factor_diagnostic.jl)**: Detailed wind factor behavior analysis

### Running Examples

```bash
julia --project=. examples/basic_examples.jl
```

This generates plots showing:
- Wind speed vs rate of spread (exponential relationship)
- Fuel moisture vs rate of spread (inverse exponential)
- Slope vs rate of spread (quadratic relationship)

## Validation

The model has been thoroughly validated against expected fire behavior:

| Test | Result | Status |
|------|--------|--------|
| Wind sensitivity | Exponential increase with wind speed | ✅ Validated |
| Moisture sensitivity | Exponential decay with moisture | ✅ Validated |
| Slope effect | Quadratic increase with slope | ✅ Validated |
| Fuel type differences | Realistic spread rate variations | ✅ Validated |
| Physical bounds | All values positive and reasonable | ✅ Validated |

See [FINAL_VALIDATION.md](FINAL_VALIDATION.md) for detailed validation results.

### Example Results

**Baseline (grass fuel, 4 mph wind, 5% moisture, flat):**
- Rate of spread: 6.34 ft/min (116 m/hr)

**Wind effect (15 mph):**
- Rate of spread: 100 ft/min (1830 m/hr)
- 158x increase compared to calm conditions

**Slope effect (50% upslope):**
- Rate of spread: ~50 ft/min
- 7.8x increase compared to flat ground

## Input Parameters

### Fuel Particle Parameters
- **h**: Low heat content (Btu/lb) - typically ~8000
- **S_T**: Total mineral content (fraction) - typically 0.0555
- **S_e**: Effective mineral content (fraction) - typically 0.010
- **ρ_p**: Particle density (lb/ft³) - typically 32 for wood

### Fuel Array Parameters
- **σ**: Surface-area-to-volume ratio (ft²/ft³)
  - Fine grass: 3000-3500
  - Medium fuels: 1500-2000
  - Coarse fuels: 100-1000
- **w_o**: Oven-dry fuel load (lb/ft²)
- **δ**: Fuel bed depth (ft)
- **M_x**: Dead fuel moisture of extinction - typically 0.12-0.40

### Environmental Parameters
- **M_f**: Fuel moisture content (fraction, dry weight basis)
- **U**: Wind velocity at midflame height (ft/min)
  - Conversion: 1 mph ≈ 88 ft/min
- **tan_ϕ**: Slope steepness (rise/run)
  - Example: 30% slope → tan_ϕ = 0.30

## Unit Conversions

### Rate of Spread
- 1 ft/min = 0.3048 m/min = 18.288 m/hr
- 1 ft/min = 1.829 chains/hr (wildland fire standard)

### Wind Speed
- 1 mph ≈ 88 ft/min
- 1 m/s ≈ 196.85 ft/min

### Fuel Load
- 1 lb/ft² = 4.88 kg/m²
- 1 ton/acre = 0.2242 kg/m²

## Heterogeneous Fuel Model

For fuel beds with multiple size classes:

```julia
using WildlandFire

# Create model with 2 dead classes and 1 live class
sys = create_heterogeneous_rothermel_model(2, 1)

# Define parameters for each size class
params = [
    # Dead fuel class 1 (1-hr timelag)
    h_d[1] => 8000.0,
    σ_d[1] => 3500.0,
    (w_o)_d[1] => 0.10,
    # ... more parameters
]

# Solve as usual
prob = ODEProblem(sys, [], (0.0, 1.0), params)
sol = solve(prob)
```

## API Reference

### Exported Symbols

**Model Systems:**
- `rothermel_simplified`: Basic Rothermel model (simplified ODESystem)
- `rothermel_system`: Basic Rothermel model (full ODESystem)
- `create_heterogeneous_rothermel_model(n_dead, n_live)`: Multi-class fuel model

**Variables (for accessing solution):**
- `R`: Rate of spread (ft/min)
- `I_R`: Reaction intensity (Btu/ft²·min)
- `ϕ_w`: Wind factor
- `ϕ_s`: Slope factor
- `β`: Packing ratio
- `η_M`: Moisture damping coefficient
- And many more...

**Parameters:**
- `h`, `σ`, `w_o`, `δ`, `M_f`, `M_x`, `U`, `tan_ϕ`, etc.

## Testing

Run the test suite:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## References

1. **Rothermel, R.C.** 1972. A mathematical model for predicting fire spread in wildland fuels. USDA Forest Service Research Paper INT-115.

2. **Albini, F.A.** 1976a. Estimating wildfire behavior and effects. USDA Forest Service General Technical Report INT-30.

3. **Andrews, P.L.** 2018. The Rothermel surface fire spread model and associated developments: A comprehensive explanation. USDA Forest Service General Technical Report RMRS-GTR-371.

## Model Limitations

These are inherent limitations of the Rothermel model:

1. **Steady-state**: Does not capture temporal fire dynamics
2. **One-dimensional**: Assumes uniform fire front
3. **Surface fire only**: No crown fire or spotting
4. **No fuel consumption**: Does not track burnout over time
5. **Homogeneous fuel bed** (basic model): Extension available for heterogeneous fuels

## Contributing

Contributions are welcome! Please ensure:
- Code follows Julia style guidelines
- Tests pass and new features include tests
- Documentation is updated
- Physical accuracy is maintained

## License

MIT License - Copyright (c) 2026 EarthSciML authors and contributors

See [LICENSE](LICENSE) for details.

## Acknowledgments

This implementation is based on the foundational work of Richard C. Rothermel and Frank A. Albini on wildland fire behavior modeling.

## Citation

If you use this package in research, please cite:

```bibtex
@software{wildlandfire_jl,
  author = {EarthSciML authors and contributors},
  title = {WildlandFire.jl: Wildland Fire Behavior Modeling in Julia},
  year = {2026},
  url = {https://github.com/your-org/WildlandFire.jl}
}
```

And the original Rothermel paper:

```bibtex
@techreport{rothermel1972,
  author = {Rothermel, Richard C.},
  title = {A mathematical model for predicting fire spread in wildland fuels},
  institution = {USDA Forest Service},
  year = {1972},
  number = {INT-115},
  type = {Research Paper}
}
```
