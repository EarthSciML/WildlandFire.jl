# Model Description

This page provides a detailed description of the Rothermel surface fire spread model as implemented in WildlandFire.jl.

## Overview

The Rothermel model is a semi-empirical mathematical model that predicts the forward rate of spread of a surface fire through wildland fuels. It was developed by Richard C. Rothermel in 1972 at the USDA Forest Service Fire Laboratory in Missoula, Montana.

## Main Equation

The fundamental equation for rate of spread is:

```math
R = \frac{I_R \xi (1 + \phi_w + \phi_s)}{\rho_b \varepsilon Q_{ig}}
```

Where:
- **R**: Rate of spread (ft/min)
- **I_R**: Reaction intensity - heat release rate per unit area (Btu/ft²·min)
- **ξ** (xi): Propagating flux ratio - fraction of heat that advances the fire (dimensionless)
- **φ_w** (phi_w): Wind factor - effect of wind on spread rate (dimensionless)
- **φ_s** (phi_s): Slope factor - effect of slope on spread rate (dimensionless)
- **ρ_b**: Bulk density of the fuel bed (lb/ft³)
- **ε** (epsilon): Effective heating number - accounts for fuel particle size (dimensionless)
- **Q_ig**: Heat of preignition - energy required to ignite fuel (Btu/lb)

This equation represents a balance between the **heat source** (numerator) and the **heat sink** (denominator).

## Physical Interpretation

### Heat Source: I_R ξ (1 + φ_w + φ_s)

The numerator represents the energy available to preheat and ignite unburned fuel ahead of the fire front.

1. **I_R**: Total heat release from burning fuel
2. **ξ**: Fraction of that heat that actually advances the fire (not lost to radiation upward)
3. **(1 + φ_w + φ_s)**: Enhancement due to wind and slope

**Wind effect (φ_w)**: Wind tilts flames forward, preheating fuels ahead of the fire and providing more oxygen.

**Slope effect (φ_s)**: Upslope spread is faster because flames are closer to unburned fuel, increasing radiant heating.

### Heat Sink: ρ_b ε Q_ig

The denominator represents the energy required to bring the fuel to ignition temperature.

1. **ρ_b**: Amount of fuel per unit volume
2. **ε**: Adjustment for how efficiently heat penetrates fuel particles
3. **Q_ig**: Energy needed per unit mass of fuel (increases with moisture)

## Model Components

### 1. Reaction Intensity (I_R)

The rate at which fuel burns and releases heat:

```math
I_R = \Gamma w_n h \eta_M \eta_s
```

Where:
- **Γ** (Gamma): Optimum reaction velocity (min⁻¹)
- **w_n**: Net fuel load, accounting for mineral content (lb/ft²)
- **h**: Heat content of the fuel (Btu/lb)
- **η_M**: Moisture damping coefficient (0-1)
- **η_s**: Mineral damping coefficient (0-1)

#### Reaction Velocity (Γ)

The reaction velocity determines how quickly fuel burns:

```math
\Gamma = \Gamma_{max} \left(\frac{\beta}{\beta_{opt}}\right)^A \exp\left[A\left(1 - \frac{\beta}{\beta_{opt}}\right)\right]
```

Where:
- **Γ_max**: Maximum reaction velocity, function of σ
- **β**: Packing ratio (ρ_b / ρ_p) - how tightly fuel is packed
- **β_opt**: Optimum packing ratio for maximum burning
- **A**: Constant that depends on σ

Reaction velocity peaks at β_opt. Too loose (β << β_opt) means poor heat transfer; too tight (β >> β_opt) means limited oxygen.

#### Damping Coefficients

**Moisture Damping (η_M)**: Reduces reaction intensity as fuel moisture increases:

```math
\eta_M = 1 - 2.59\frac{M_f}{M_x} + 5.11\left(\frac{M_f}{M_x}\right)^2 - 3.52\left(\frac{M_f}{M_x}\right)^3
```

Where:
- **M_f**: Fuel moisture content (fraction of dry weight)
- **M_x**: Moisture of extinction - moisture level at which fire will not spread

When M_f ≥ M_x, η_M = 0 and the fire cannot spread.

**Mineral Damping (η_s)**: Reduces reaction intensity due to inert mineral content:

```math
\eta_s = 0.174 S_e^{-0.19}
```

Where S_e is the effective mineral content (silica-free ash).

### 2. Propagating Flux Ratio (ξ)

The fraction of reaction intensity that heats adjacent unburned fuel:

```math
\xi = \frac{\exp\left[(0.792 + 0.681\sigma^{0.5})(\beta + 0.1)\right]}{192 + 0.2595\sigma}
```

Where σ is the surface-area-to-volume ratio of fuel particles.

Finer fuels (higher σ) have lower ξ because more heat is lost to radiation.

### 3. Wind Factor (φ_w)

The wind factor describes the exponential increase in spread rate with wind speed:

```math
\phi_w = C U^B \left(\frac{\beta}{\beta_{opt}}\right)^{-E}
```

Where:
- **U**: Wind velocity at midflame height (ft/min)
- **C, B, E**: Functions of σ and β

Key features:
- Wind effect is roughly exponential (U^B where B ≈ 1-2)
- Effect is stronger for looser fuel beds (lower β)
- Effect is stronger for finer fuels (higher σ)

### 4. Slope Factor (φ_s)

The slope factor describes the quadratic increase in upslope spread rate:

```math
\phi_s = 5.275 \beta^{-0.3} (\tan \phi)^2
```

Where:
- **tan φ**: Slope steepness (rise/run)

Key features:
- Slope effect is quadratic - doubling slope quadruples the effect
- Effect is stronger for looser fuel beds (lower β)
- Only affects upslope spread; downslope spread is slower

### 5. Heat of Preignition (Q_ig)

Energy required to bring fuel from ambient temperature to ignition:

```math
Q_{ig} = 250 + 1116 M_f
```

Increases linearly with fuel moisture because energy is needed to evaporate water.

### 6. Effective Heating Number (ε)

Accounts for the fact that only the surface of fuel particles is initially heated:

```math
\varepsilon = \exp\left(-\frac{138}{\sigma}\right)
```

- **ε → 1** for very fine fuels (high σ): entire particle heats quickly
- **ε → 0** for coarse fuels (low σ): only surface heats initially

## Fuel Bed Characteristics

### Packing Ratio (β)

```math
\beta = \frac{\rho_b}{\rho_p} = \frac{w_o}{\delta \rho_p}
```

The packing ratio describes how tightly fuel is packed:
- **Low β** (0.001-0.01): Loose, airy fuel beds (grass)
- **Medium β** (0.01-0.05): Moderate density (brush, litter)
- **High β** (0.05-0.15): Tightly packed (compacted fuels)

### Optimum Packing Ratio (β_opt)

```math
\beta_{opt} = 3.348 \sigma^{-0.8189}
```

The packing ratio that gives maximum reaction velocity:
- Finer fuels (high σ) have lower β_opt
- Coarser fuels (low σ) have higher β_opt

### Surface-Area-to-Volume Ratio (σ)

Describes fuel particle fineness (ft²/ft³):
- **Fine grass**: 3000-4000 ft²/ft³
- **Leaves, pine needles**: 1500-2500 ft²/ft³
- **Small twigs**: 500-1500 ft²/ft³
- **Large branches**: 100-500 ft²/ft³

Finer fuels ignite more easily but burn faster.

## Model Assumptions

The Rothermel model makes several important assumptions:

1. **Steady-state**: Fire has reached an equilibrium spread rate
2. **Homogeneous fuel**: Fuel properties are uniform across the landscape
3. **One-dimensional spread**: Fire spreads as a uniform front in one direction
4. **Surface fire**: Model does not account for crown fire or spotting
5. **No fuel consumption**: Does not track fuel burnout or residual burning
6. **Continuous fuel bed**: No gaps or discontinuities in fuel

## Model Limitations

Understanding the limitations is important for proper application:

### Physical Limitations
- Does not predict fire intensity or flame length directly
- No temporal dynamics (transient behavior)
- No interaction between multiple fires
- No spotting or ember transport
- No crown fire or torching

### Fuel Limitations
- Best suited for relatively homogeneous fuel beds
- Heterogeneous fuel extension is complex
- Does not account for fuel moisture gradients
- Does not model fuel consumption or burnout

### Environmental Limitations
- Wind should be at midflame height (requires adjustment from standard measurements)
- Does not account for atmospheric instability
- No feedback from fire to atmosphere
- Slope effect is empirical and may not apply to very steep slopes (>60%)

### Operational Limitations
- Requires detailed fuel data
- Sensitive to parameter uncertainty
- Best used for relative comparisons rather than absolute predictions
- Should be validated against local conditions when possible

## Validation and Accuracy

The Rothermel model has been extensively validated:

### Strengths
- Widely validated against experimental burns
- Foundation of operational fire behavior systems (BEHAVE, BehavePlus, FARSITE)
- Captures key physical processes
- Generally accurate for low to moderate intensity fires

### Known Issues
- May underpredict spread rate in very high winds (>25 mph)
- Slope effect becomes less reliable on very steep slopes
- Wind-slope interaction is additive (may not be physically accurate)
- Does not account for fire-atmosphere coupling in extreme fires

### Recommended Use
- Relative comparisons between scenarios
- Fire danger rating
- Fuel management planning
- Initial attack planning
- Training and education

The model should be used as one tool among many, with local knowledge and experience informing decisions.

## References

1. **Rothermel, R.C.** 1972. A mathematical model for predicting fire spread in wildland fuels. USDA Forest Service Research Paper INT-115.

2. **Albini, F.A.** 1976a. Estimating wildfire behavior and effects. USDA Forest Service General Technical Report INT-30.

3. **Andrews, P.L.** 2018. The Rothermel surface fire spread model and associated developments: A comprehensive explanation. USDA Forest Service General Technical Report RMRS-GTR-371.

4. **Anderson, H.E.** 1982. Aids to determining fuel models for estimating fire behavior. USDA Forest Service General Technical Report INT-122.
