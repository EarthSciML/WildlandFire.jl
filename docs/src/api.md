# API Reference

This page provides detailed documentation for all exported functions, variables, and parameters in WildlandFire.jl.

## Model Creation Functions

### Rothermel

```julia
Rothermel(; name=:rothermel, simplify=true)
```

Factory function to create a Rothermel surface fire spread model.

**Arguments:**
- `name::Symbol=:rothermel`: Name for the ODESystem
- `simplify::Bool=true`: Whether to apply structural_simplify to the system

**Returns:**
- `ODESystem`: A ModelingToolkit ODESystem representing the Rothermel fire spread model

**Description:**

Creates a complete Rothermel surface fire spread model including all equations from Rothermel (1972) and Albini (1976a). The model calculates the forward rate of spread based on fuel properties, environmental conditions, and topography.

The main equation is:
```
R = I_R·ξ(1 + ϕ_w + ϕ_s) / (ρ_b·ε·Q_ig)
```

**Example:**
```@example api1
using WildlandFire
using OrdinaryDiffEqDefault

# Create model
sys = Rothermel()

# Define parameters
params = [
    h => 8000.0,
    σ => 3500.0,
    w_o => 0.138,
    δ => 1.0,
    M_f => 0.05,
    U => 352.0,
    tan_ϕ => 0.0,
    S_T => 0.0555,
    S_e => 0.010,
    ρ_p => 32.0,
    M_x => 0.12,
]

# Solve
prob = ODEProblem(sys, [], (0.0, 1.0), params)
sol = solve(prob)

# Extract result
R_ftmin = sol[R][end]
```

---

### rothermel_simplified

```julia
rothermel_simplified::ODESystem
```

Pre-built simplified Rothermel model (backward compatibility).

**Description:**

This is a pre-instantiated version of `Rothermel(name=:rothermel, simplify=true)` provided for backward compatibility. Use this if you want to use the model without explicitly calling the factory function.

**Example:**
```julia
using WildlandFire
using OrdinaryDiffEqDefault

sys = rothermel_simplified  # Use pre-built model

params = [h => 8000.0, σ => 3500.0, ...]  # ... represents other parameters
prob = ODEProblem(sys, [], (0.0, 1.0), params)
sol = solve(prob)
```

---

### create_heterogeneous_rothermel_model

```julia
create_heterogeneous_rothermel_model(n_dead::Int, n_live::Int; name=:rothermel_hetero)
```

Creates a Rothermel model with multiple fuel size classes.

**Arguments:**
- `n_dead::Int`: Number of dead fuel size classes
- `n_live::Int`: Number of live fuel size classes
- `name::Symbol=:rothermel_hetero`: Name for the ODESystem

**Returns:**
- `ODESystem`: Extended Rothermel model with heterogeneous fuels

**Description:**

Creates an extended Rothermel model that handles fuel beds with multiple size classes. Each size class has its own fuel properties (heat content, SAV ratio, fuel load, moisture, etc.). The model computes weighted averages based on the contribution of each size class to the total reaction intensity.

**Example:**
```julia
using WildlandFire
using OrdinaryDiffEqDefault

# Create model with 2 dead classes and 1 live class
sys = create_heterogeneous_rothermel_model(2, 1)

# Define parameters for each size class
params = [
    # Dead fuel class 1 (1-hr timelag)
    h_d[1] => 8000.0,
    σ_d[1] => 3500.0,
    w_o_d[1] => 0.10,
    M_f_d[1] => 0.05,

    # Dead fuel class 2 (10-hr timelag)
    h_d[2] => 8000.0,
    σ_d[2] => 109.0,
    w_o_d[2] => 0.05,
    M_f_d[2] => 0.07,

    # Live fuel class 1
    h_l[1] => 8000.0,
    σ_l[1] => 1500.0,
    w_o_l[1] => 0.15,
    M_f_l[1] => 1.0,

    # Common parameters
    δ => 1.0,
    U => 352.0,
    tan_ϕ => 0.0,
]

prob = ODEProblem(sys, [], (0.0, 1.0), params)
sol = solve(prob)
```

---

## Model Variables

These variables are available in the solution after solving the ODE system.

### Primary Output

#### R
```julia
R(t)
```
**Rate of spread** (ft/min) - The forward rate of spread of the fire front.

**Access in solution:**
```julia
sol[R][end]  # Get final value
sol[R]       # Get all values over time
```

**Typical range:** 0.1 - 100 ft/min (0 - 1830 m/hr)

---

### Intermediate Calculations

#### I_R
```julia
I_R(t)
```
**Reaction intensity** (Btu/ft²·min) - Rate of heat release per unit area of the fire front.

**Formula:** `I_R = Γ·w_n·h·η_M·η_s`

**Typical range:** 1,000 - 10,000 Btu/ft²·min

---

#### Γ
```julia
Γ(t)
```
**Optimum reaction velocity** (min⁻¹) - Rate at which fuel is consumed by combustion.

**Formula:** `Γ = Γ_max(β/β_opt)^A * exp[A(1 - β/β_opt)]`

**Typical range:** 1 - 20 min⁻¹

---

#### Γ_max
```julia
Γ_max(t)
```
**Maximum reaction velocity** (min⁻¹) - Maximum possible reaction velocity for the fuel.

**Formula:** `Γ_max = σ^1.5 / (495 + 0.0594σ^1.5)`

---

#### β
```julia
β(t)
```
**Packing ratio** (dimensionless) - Ratio of bulk density to particle density.

**Formula:** `β = ρ_b / ρ_p = w_o / (δ·ρ_p)`

**Typical range:** 0.001 - 0.01

---

#### β_opt
```julia
β_opt(t)
```
**Optimum packing ratio** (dimensionless) - Packing ratio at which reaction velocity is maximized.

**Formula:** `β_opt = 3.348σ^(-0.8189)`

---

#### β_ratio
```julia
β_ratio(t)
```
**Relative packing ratio** (dimensionless) - Ratio of actual to optimum packing ratio.

**Formula:** `β_ratio = β / β_opt`

---

#### ρ_b
```julia
ρ_b(t)
```
**Bulk density** (lb/ft³) - Mass of fuel per unit volume of fuel bed.

**Formula:** `ρ_b = w_o / δ`

**Typical range:** 0.1 - 10 lb/ft³

---

#### w_n
```julia
w_n(t)
```
**Net fuel load** (lb/ft²) - Oven-dry fuel load with mineral content removed.

**Formula:** `w_n = w_o(1 - S_T)`

---

#### η_M
```julia
η_M(t)
```
**Moisture damping coefficient** (dimensionless) - Reduction in reaction intensity due to fuel moisture.

**Formula:** `η_M = 1 - 2.59(M_f/M_x) + 5.11(M_f/M_x)² - 3.52(M_f/M_x)³`

**Range:** 0 to 1 (1 = dry, 0 = extinction)

---

#### η_s
```julia
η_s(t)
```
**Mineral damping coefficient** (dimensionless) - Reduction in reaction intensity due to mineral content.

**Formula:** `η_s = 0.174S_e^(-0.19)` (max = 1.0)

**Range:** 0 to 1

---

#### ξ
```julia
ξ(t)
```
**Propagating flux ratio** (dimensionless) - Proportion of reaction intensity that heats adjacent unburned fuel.

**Formula:** `ξ = exp[(0.792 + 0.681σ^0.5)(β + 0.1)] / (192 + 0.2595σ)`

**Typical range:** 0.01 - 0.5

---

#### ϕ_w
```julia
ϕ_w(t)
```
**Wind factor** (dimensionless) - Dimensionless multiplier for wind effect on spread rate.

**Formula:** `ϕ_w = C·U^B·(β/β_opt)^(-E)` where:
- `C = 7.47exp(-0.133σ^0.55)`
- `B = 0.02526σ^0.54`
- `E = 0.715exp(-3.59×10^(-4)σ)`

**Range:** 0 to 100+ (can be very large at high winds)

---

#### ϕ_s
```julia
ϕ_s(t)
```
**Slope factor** (dimensionless) - Dimensionless multiplier for slope effect on upslope spread rate.

**Formula:** `ϕ_s = 5.275β^(-0.3)(tan ϕ)²`

**Range:** 0 to 10+ (quadratic with slope)

---

#### ε
```julia
ε(t)
```
**Effective heating number** (dimensionless) - Proportion of fuel particles that are effectively heated by flames.

**Formula:** `ε = exp(-138/σ)`

**Range:** 0 to 1 (approaches 1 for fine fuels)

---

#### Q_ig
```julia
Q_ig(t)
```
**Heat of preignition** (Btu/lb) - Energy required to ignite fuel per unit mass.

**Formula:** `Q_ig = 250 + 1116M_f`

**Typical range:** 250 - 500 Btu/lb

---

#### R_0
```julia
R_0(t)
```
**No-wind, no-slope rate of spread** (ft/min) - Rate of spread on flat ground with no wind.

**Formula:** `R_0 = I_R·ξ / (ρ_b·ε·Q_ig)`

---

## Model Parameters

These parameters must be provided when solving the model.

### Fuel Particle Parameters

#### h
```julia
h
```
**Low heat content** (Btu/lb) - Energy released per unit mass of fuel burned.

**Typical value:** 8000 Btu/lb (for most wood fuels)

**Range:** 7500 - 9500 Btu/lb

---

#### S_T
```julia
S_T
```
**Total mineral content** (fraction) - Total mineral fraction in fuel.

**Typical value:** 0.0555 (5.55%)

**Range:** 0.01 - 0.10

---

#### S_e
```julia
S_e
```
**Effective mineral content** (fraction) - Mineral fraction that affects combustion.

**Typical value:** 0.010 (1.0%)

**Range:** 0.001 - 0.05

---

#### ρ_p
```julia
ρ_p
```
**Particle density** (lb/ft³) - Oven-dry density of fuel particles.

**Typical value:** 32 lb/ft³ (for wood)

**Range:** 25 - 40 lb/ft³

---

### Fuel Array Parameters

#### σ
```julia
σ
```
**Surface-area-to-volume ratio** (ft²/ft³) - Surface area of fuel particles per unit volume.

**Typical values:**
- Fine grass: 3000 - 3500 ft²/ft³
- Medium fuels: 1500 - 2000 ft²/ft³
- Coarse fuels: 100 - 1000 ft²/ft³

**Critical parameter:** Dominates fire behavior characteristics

---

#### w_o
```julia
w_o
```
**Oven-dry fuel load** (lb/ft²) - Mass of fuel per unit area of ground.

**Typical range:** 0.01 - 0.50 lb/ft²

**Conversions:**
- 1 lb/ft² = 4.88 kg/m²
- 1 ton/acre = 0.2242 kg/m² = 0.046 lb/ft²

---

#### δ
```julia
δ
```
**Fuel bed depth** (ft) - Vertical depth of the fuel layer.

**Typical range:** 0.5 - 10 ft

**Examples:**
- Short grass: 1 ft
- Tall grass: 2.5 ft
- Brush: 2 - 6 ft
- Timber litter: 0.5 - 1 ft

---

#### M_x
```julia
M_x
```
**Moisture of extinction** (fraction) - Fuel moisture above which fire will not spread.

**Typical values:**
- Grass: 0.12 (12%)
- Timber: 0.15 - 0.25 (15-25%)
- Chaparral: 0.20 - 0.40 (20-40%)

**Range:** 0.10 - 0.40

---

### Environmental Parameters

#### M_f
```julia
M_f
```
**Fuel moisture content** (fraction, dry weight basis) - Mass of water per unit dry mass of fuel.

**Formula:** `M_f = (wet_mass - dry_mass) / dry_mass`

**Critical thresholds:**
- < 0.05 (5%): Extreme fire danger
- 0.05 - 0.10: High fire danger
- 0.10 - 0.20: Moderate fire danger
- > M_x: No spread

---

#### U
```julia
U
```
**Wind velocity at midflame height** (ft/min) - Wind speed at the height of the flame midpoint.

**Conversions:**
- 1 mph ≈ 88 ft/min
- 1 m/s ≈ 196.85 ft/min

**Typical values:**
- Calm: 0 - 88 ft/min (0-1 mph)
- Light: 88 - 440 ft/min (1-5 mph)
- Moderate: 440 - 880 ft/min (5-10 mph)
- Strong: 880 - 1760 ft/min (10-20 mph)

**Critical parameter:** Exponential effect on spread rate

---

#### tan_ϕ
```julia
tan_ϕ
```
**Slope steepness** (fraction) - Tangent of slope angle (rise/run).

**Formula:** `tan_ϕ = vertical_rise / horizontal_distance`

**Conversions:**
- Percent slope: `tan_ϕ = slope% / 100`
- Degrees: `tan_ϕ = tan(angle_degrees * π/180)`

**Examples:**
- 10% slope → tan_ϕ = 0.10
- 30% slope → tan_ϕ = 0.30
- 45° slope → tan_ϕ = 1.00

**Range:** 0 - 1.0 (0-100% slope or 0-45°)

---

### Wind Factor Constants

#### C_const
```julia
C_const = 7.47
```
Constant used in wind factor calculation. Default value is 7.47.

#### B_const
```julia
B_const = 0.02526
```
Constant used in wind factor calculation. Default value is 0.02526.

#### E_const
```julia
E_const = 0.715
```
Constant used in wind factor calculation. Default value is 0.715.

---

## Usage Patterns

### Basic Simulation

```@example apiusage
using WildlandFire
using OrdinaryDiffEqDefault

sys = Rothermel()
params = [
    h => 8000.0, σ => 3500.0, w_o => 0.138, δ => 1.0,
    M_f => 0.05, U => 352.0, tan_ϕ => 0.0,
    S_T => 0.0555, S_e => 0.010, ρ_p => 32.0, M_x => 0.12
]
prob = ODEProblem(sys, [], (0.0, 1.0), params)
sol = solve(prob)
R_ftmin = sol[R][end]
println("Rate of spread: $(round(R_ftmin, digits=2)) ft/min")
```

### Accessing Multiple Variables

```@example apiusage
ros = sol[R][end]           # Rate of spread
intensity = sol[I_R][end]   # Reaction intensity
wind_effect = sol[ϕ_w][end] # Wind factor
slope_effect = sol[ϕ_s][end] # Slope factor

println("ROS: $(round(ros, digits=2)) ft/min")
println("Intensity: $(round(intensity, digits=2)) Btu/ft²·min")
println("Wind factor: $(round(wind_effect, digits=2))")
println("Slope factor: $(round(slope_effect, digits=2))")
```

### Parametric Studies

```@example apiusage
wind_speeds = 0:88:880  # 0-10 mph
results = map(wind_speeds) do U_val
    prob = ODEProblem(sys, [], (0.0, 1.0), [params..., U => U_val])
    sol = solve(prob)
    sol[R][end]
end
println("Wind speed study: $(length(results)) simulations completed")
```

---

## Unit Conversions

### Rate of Spread

```@example apiusage
R_ftmin = sol[R][end]               # ft/min (native units)
R_mhr = R_ftmin * 18.288            # m/hr
R_chains = R_ftmin * 1.829          # chains/hr (wildland standard)
R_kmhr = R_ftmin * 0.018288         # km/hr

println("Rate of spread conversions:")
println("  $(round(R_ftmin, digits=2)) ft/min")
println("  $(round(R_mhr, digits=2)) m/hr")
println("  $(round(R_chains, digits=2)) chains/hr")
println("  $(round(R_kmhr, digits=2)) km/hr")
```

### Wind Speed

```@example apiusage
wind_mph = 10
wind_ftmin_from_mph = wind_mph * 88.0        # Convert to ft/min

wind_ms = 5
wind_ftmin_from_ms = wind_ms * 196.85        # Convert to ft/min

println("Wind speed conversions:")
println("  $(wind_mph) mph = $(round(wind_ftmin_from_mph, digits=2)) ft/min")
println("  $(wind_ms) m/s = $(round(wind_ftmin_from_ms, digits=2)) ft/min")
```

### Fuel Load

```@example apiusage
load_lb_ft2 = 0.138                 # lb/ft²
load_kg_m2 = load_lb_ft2 * 4.88    # kg/m²

load_ton_acre = 1.0
load_lb_ft2_from_ton = load_ton_acre * 0.046  # Convert to lb/ft²

println("Fuel load conversions:")
println("  $(load_lb_ft2) lb/ft² = $(round(load_kg_m2, digits=2)) kg/m²")
println("  $(load_ton_acre) ton/acre = $(load_lb_ft2_from_ton) lb/ft²")
```

---

## See Also

- [Getting Started](getting_started.md) - First simulation tutorial
- [Usage Guide](usage_guide.md) - Comprehensive usage examples
- [Examples](examples/basic.md) - Detailed example scenarios
