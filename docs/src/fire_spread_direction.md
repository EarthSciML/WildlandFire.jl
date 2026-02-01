# Fire Spread Direction

## Overview

When wind blows at an angle to the slope, vector addition is used to find the direction and rate of maximum fire spread. This module implements the equations from Section 6 of Andrews (2018), which describe fire behavior in the direction of maximum spread, elliptical fire shape from a single ignition point, and fire spread normal to the perimeter.

**Reference**: Andrews, Patricia L. 2018. The Rothermel surface fire spread model and associated developments: A comprehensive explanation. Gen. Tech. Rep. RMRS-GTR-371. Fort Collins, CO: U.S. Department of Agriculture, Forest Service, Rocky Mountain Research Station. Pages 84-92, Table 26.

Additional reference for fire perimeter equations: Catchpole, E. A.; deMestre, N. J.; Gill, A. M. 1982. Intensity of fire at its perimeter. Australian Forest Research. 12: 47-54.

## FireSpreadDirection

Calculates fire behavior in the direction of maximum spread when wind is not blowing directly upslope.

```@docs
FireSpreadDirection
```

### State Variables

```@example fire_spread
using DataFrames, ModelingToolkit, Symbolics, DynamicQuantities
using WildlandFire

sys = FireSpreadDirection()

vars = unknowns(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [string(dimension(ModelingToolkit.get_unit(v))) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

### Parameters

```@example fire_spread
params = parameters(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params],
    :Units => [string(dimension(ModelingToolkit.get_unit(p))) for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

### Equations

```@example fire_spread
eqs = equations(sys)
```

## EllipticalFireSpread

Calculates fire spread from a single ignition point assuming an elliptical fire shape.

```@docs
EllipticalFireSpread
```

### State Variables

```@example elliptical
using DataFrames, ModelingToolkit, Symbolics, DynamicQuantities
using WildlandFire

sys = EllipticalFireSpread()

vars = unknowns(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [string(dimension(ModelingToolkit.get_unit(v))) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

### Parameters

```@example elliptical
params = parameters(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params],
    :Units => [string(dimension(ModelingToolkit.get_unit(p))) for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

### Equations

```@example elliptical
eqs = equations(sys)
```

## FirePerimeterSpread

Calculates fire spread rate normal to the fire perimeter based on Catchpole et al. (1982).

```@docs
FirePerimeterSpread
```

### State Variables

```@example perimeter
using DataFrames, ModelingToolkit, Symbolics, DynamicQuantities
using WildlandFire

sys = FirePerimeterSpread()

vars = unknowns(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [string(dimension(ModelingToolkit.get_unit(v))) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

### Parameters

```@example perimeter
params = parameters(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params],
    :Units => [string(dimension(ModelingToolkit.get_unit(p))) for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

### Equations

```@example perimeter
eqs = equations(sys)
```

## Analysis

### Effect of Wind Direction on Fire Spread

The following example demonstrates how wind direction affects the direction of maximum fire spread. When wind blows upslope (omega = 0), the fire spreads directly uphill. As the wind direction rotates, the direction of maximum spread shifts.

```@example analysis
using WildlandFire, ModelingToolkit, OrdinaryDiffEqDefault
using Plots, DynamicQuantities

# Create and compile the system
sys = FireSpreadDirection()
compiled_sys = mtkcompile(sys)

# Common parameters
R0_val = 0.01  # m/s
φw_val = 5.0   # wind factor
φs_val = 2.0   # slope factor
t_val = 60.0   # elapsed time (s)

# Sweep wind direction from 0 to 2pi
ω_values = range(0, 2π, length=37)
α_results = Float64[]
R_H_results = Float64[]

for ω_val in ω_values
    prob = NonlinearProblem(compiled_sys, Dict(
        compiled_sys.R0 => R0_val,
        compiled_sys.φw => φw_val,
        compiled_sys.φs => φs_val,
        compiled_sys.ω => ω_val,
        compiled_sys.elapsed_time => t_val,
        compiled_sys.β_ratio => 0.5,
        compiled_sys.C_coeff => 7.47,
        compiled_sys.B_coeff => 0.15566,
        compiled_sys.E_coeff => 0.715,
    ))
    sol = solve(prob)
    push!(α_results, ustrip(sol[compiled_sys.α]))
    push!(R_H_results, ustrip(sol[compiled_sys.R_H]))
end

p1 = plot(rad2deg.(ω_values), rad2deg.(α_results),
    xlabel="Wind direction from upslope (degrees)",
    ylabel="Direction of max spread (degrees)",
    title="Direction of Maximum Spread vs Wind Direction",
    legend=false,
    linewidth=2)

p2 = plot(rad2deg.(ω_values), R_H_results .* 60,  # Convert to m/min
    xlabel="Wind direction from upslope (degrees)",
    ylabel="Head fire rate of spread (m/min)",
    title="Head Fire Rate of Spread vs Wind Direction",
    legend=false,
    linewidth=2)

plot(p1, p2, layout=(2,1), size=(600, 500))
```

### Elliptical Fire Shape

This example shows the elliptical fire perimeter for different length-to-width ratios.

```@example analysis
using WildlandFire, ModelingToolkit, OrdinaryDiffEqDefault
using Plots, DynamicQuantities

sys = EllipticalFireSpread()
compiled_sys = mtkcompile(sys)

R_H_val = 0.1  # m/s
t_val = 3600.0  # 1 hour

p = plot(xlabel="Crosswind distance (m)", ylabel="Upwind/Downwind distance (m)",
         title="Fire Ellipse Shapes for Different Z Values",
         aspect_ratio=:equal, legend=:topleft)

for Z_val in [1.5, 2.0, 3.0, 5.0]
    prob = NonlinearProblem(compiled_sys, Dict(
        compiled_sys.R_H => R_H_val,
        compiled_sys.Z => Z_val,
        compiled_sys.γ => 0.0,
        compiled_sys.elapsed_time => t_val,
    ))
    sol = solve(prob)

    # Get ellipse dimensions
    L = ustrip(sol[compiled_sys.L])
    W = ustrip(sol[compiled_sys.W])
    D_H = ustrip(sol[compiled_sys.D_H])
    D_B = ustrip(sol[compiled_sys.D_B])

    # Draw ellipse centered at ignition point (0, 0)
    # Semi-axes: a = L/2 (upwind/downwind), b = W/2 (crosswind)
    θ = range(0, 2π, length=100)
    a = L / 2  # semi-major axis
    b = W / 2  # semi-minor axis
    # Center offset: ignition point is at focus, center is at (D_H - a, 0)
    center_offset = D_H - a

    x_ellipse = b .* sin.(θ)
    y_ellipse = a .* cos.(θ) .+ center_offset

    plot!(p, x_ellipse, y_ellipse, label="Z = $(Z_val)", linewidth=2)
end

# Mark ignition point
scatter!(p, [0], [0], label="Ignition", markersize=8, color=:red)

p
```

### Fire Intensity Around the Perimeter

This example shows how fireline intensity varies around the fire perimeter.

```@example analysis
using WildlandFire, ModelingToolkit, OrdinaryDiffEqDefault
using Plots, DynamicQuantities

sys = FirePerimeterSpread()
compiled_sys = mtkcompile(sys)

R_H_val = 0.1  # m/s
f_val = 180.0  # semi-major axis (m)
h_val = 90.0   # semi-minor axis (m)
g_val = sqrt(f_val^2 - h_val^2)  # focal distance
H_A_val = 5000.0  # Heat per unit area (J/m²)

γ_values = range(0, 2π, length=73)
R_ψ_results = Float64[]
I_B_results = Float64[]

for γ_val in γ_values
    prob = NonlinearProblem(compiled_sys, Dict(
        compiled_sys.R_H => R_H_val,
        compiled_sys.f => f_val,
        compiled_sys.g => g_val,
        compiled_sys.h => h_val,
        compiled_sys.γ => γ_val,
        compiled_sys.H_A => H_A_val,
    ))
    sol = solve(prob)
    push!(R_ψ_results, ustrip(sol[compiled_sys.R_ψ]))
    push!(I_B_results, ustrip(sol[compiled_sys.I_B]))
end

p1 = plot(rad2deg.(γ_values), R_ψ_results .* 60,  # Convert to m/min
    xlabel="Direction from head fire (degrees)",
    ylabel="Rate of spread normal to perimeter (m/min)",
    title="Perimeter Spread Rate Around Fire",
    legend=false,
    linewidth=2)

p2 = plot(rad2deg.(γ_values), I_B_results ./ 1000,  # Convert to kW/m
    xlabel="Direction from head fire (degrees)",
    ylabel="Fireline intensity (kW/m)",
    title="Fireline Intensity Around Fire Perimeter",
    legend=false,
    linewidth=2)

plot(p1, p2, layout=(2,1), size=(600, 500))
```
