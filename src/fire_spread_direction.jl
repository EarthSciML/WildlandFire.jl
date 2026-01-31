"""
    FireSpreadDirection(; name=:FireSpreadDirection)

Create a system for calculating fire behavior in the direction of maximum spread when
wind is not blowing directly upslope.

The Rothermel surface fire spread model assumes wind and slope effects are in the same direction.
When wind blows at an angle to the slope, vector addition is used to find the direction and
rate of maximum fire spread.

# Model Description

The slope vector has magnitude D_S and direction 0 (upslope). The wind vector has magnitude
D_W in direction ω from upslope. The resultant vector gives the direction of maximum spread (α)
and the head fire rate of spread (R_H).

Key outputs:
- `R_H`: Rate of spread in direction of maximum spread (head fire)
- `α`: Direction of maximum spread relative to upslope
- `φ_E`: Effective wind factor in direction of maximum spread
- `U_E`: Effective wind speed in direction of maximum spread
- `Z`: Fire length-to-width ratio

# Reference

Andrews, Patricia L. 2018. The Rothermel surface fire spread model and associated
developments: A comprehensive explanation. Gen. Tech. Rep. RMRS-GTR-371. Fort Collins,
CO: U.S. Department of Agriculture, Forest Service, Rocky Mountain Research Station.
Table 26, pages 86-87.

# Example

```julia
using ModelingToolkit, WildlandFire, OrdinaryDiffEqDefault

sys = FireSpreadDirection()
compiled_sys = mtkcompile(sys)

# Calculate with cross-slope wind at 45 degrees from upslope
prob = NonlinearProblem(compiled_sys, Dict(
    compiled_sys.R0 => 0.01,      # No-wind no-slope rate of spread (m/s)
    compiled_sys.φw => 5.0,        # Wind factor (dimensionless)
    compiled_sys.φs => 2.0,        # Slope factor (dimensionless)
    compiled_sys.ω => π/4,         # Wind direction 45° from upslope (rad)
    compiled_sys.β_ratio => 0.5,   # Relative packing ratio
    compiled_sys.C_coeff => 7.47,  # Wind coefficient C
    compiled_sys.B_coeff => 0.5,   # Wind coefficient B
    compiled_sys.E_coeff => 0.5,   # Wind coefficient E
    compiled_sys.elapsed_time => 60.0  # Elapsed time (s)
))
sol = solve(prob)
```
"""
@component function FireSpreadDirection(; name=:FireSpreadDirection)
    @constants begin
        one = 1.0, [description = "Dimensionless one for unit balancing", unit = u"1"]
        zero_rad = 0.0, [description = "Zero angle", unit = u"rad"]
        # Conversion factor for length-to-width ratio equation
        # Z = 1 + 0.25*U_E where U_E is in mi/h
        # In SI: Z = 1 + 0.25 * U_E * 60 / 1609.34 = 1 + 0.00932 * U_E (U_E in m/s)
        c_Z = 0.00932, [description = "Length-to-width ratio coefficient (SI)", unit = u"s/m"]
        U_ref = 1.0, [description = "Reference wind speed", unit = u"m/s"]
        D_min = 1e-10, [description = "Minimum distance to avoid division by zero", unit = u"m"]
    end

    # Input parameters from RothermelFireSpread
    @parameters begin
        R0, [description = "No-wind no-slope rate of spread", unit = u"m/s"]
        φw, [description = "Wind factor (dimensionless)", unit = u"1"]
        φs, [description = "Slope factor (dimensionless)", unit = u"1"]
        ω, [description = "Wind direction relative to upslope", unit = u"rad"]
        elapsed_time, [description = "Elapsed time for vector calculations", unit = u"s"]
        # Parameters needed for effective wind speed calculation
        β_ratio, [description = "Relative packing ratio β/β_op (dimensionless)", unit = u"1"]
        C_coeff, [description = "Wind coefficient C (dimensionless)", unit = u"1"]
        B_coeff, [description = "Wind coefficient B (dimensionless)", unit = u"1"]
        E_coeff, [description = "Wind coefficient E (dimensionless)", unit = u"1"]
    end

    # Intermediate variables - vector components
    @variables begin
        D_S(t), [description = "Slope vector magnitude", unit = u"m"]
        D_W(t), [description = "Wind vector magnitude", unit = u"m"]
        X(t), [description = "X-component of resultant vector", unit = u"m"]
        Y(t), [description = "Y-component of resultant vector", unit = u"m"]
        D_H(t), [description = "Head fire distance magnitude", unit = u"m"]
    end

    # Output variables
    @variables begin
        R_H(t), [description = "Rate of spread in direction of maximum spread (head fire)", unit = u"m/s"]
        α(t), [description = "Direction of maximum spread relative to upslope", unit = u"rad"]
        φ_E(t), [description = "Effective wind factor in direction of maximum spread (dimensionless)", unit = u"1"]
        U_E(t), [description = "Effective wind speed in direction of maximum spread", unit = u"m/s"]
        Z(t), [description = "Fire length-to-width ratio (dimensionless)", unit = u"1"]
    end

    eqs = [
        # Slope and wind vector magnitudes - Table 26, Andrews (2018)
        D_S ~ R0 * φs * elapsed_time,                    # Eq. D_S
        D_W ~ R0 * φw * elapsed_time,                    # Eq. D_W

        # Resultant vector components - Table 26, Andrews (2018)
        X ~ D_S + D_W * cos(ω),                          # Eq. X
        Y ~ D_W * sin(ω),                                # Eq. Y

        # Head fire distance magnitude - Table 26, Andrews (2018)
        D_H ~ sqrt(X^2 + Y^2),                           # Eq. D_H

        # Rate of spread in direction of maximum spread - Table 26, Andrews (2018)
        R_H ~ R0 + D_H / elapsed_time,                   # Eq. R_H

        # Direction of maximum spread relative to upslope - Table 26, Andrews (2018)
        # Note: Using atan2-like behavior with abs(Y) to handle quadrants
        α ~ asin(abs(Y) / max(D_H, D_min)),              # Eq. α

        # Effective wind factor - Table 26, Andrews (2018)
        φ_E ~ R_H / R0 - one,                            # Eq. φ_E

        # Effective wind speed - Table 26, Andrews (2018)
        # U_E = [φ_E * (β/β_op)^E / C]^(1/B)
        U_E ~ U_ref * (φ_E * β_ratio^E_coeff / C_coeff)^(one / B_coeff),  # Eq. U_E

        # Length-to-width ratio - Table 26, Andrews (2018)
        # Z = L/W = 1 + 0.25*U_E (with U_E in mi/h)
        # In SI: Z = 1 + c_Z * U_E
        Z ~ one + c_Z * U_E,                             # Eq. Z
    ]

    return System(eqs, t; name)
end


"""
    EllipticalFireSpread(; name=:EllipticalFireSpread)

Create a system for calculating fire spread from a single ignition point assuming
an elliptical fire shape.

When a fire starts from a single ignition point, it is assumed to be elliptically shaped
with the fire spreading at a steady rate throughout the projection time. The length-to-width
ratio is determined from the effective wind speed in the direction of maximum spread.

# Key Outputs
- `e`: Eccentricity of the fire ellipse
- `R_γ`: Rate of spread from ignition point in direction γ
- `R_B`: Rate of spread of backing fire
- `L`: Fire length
- `W`: Fire width
- `D_F`: Flanking spread distance

# Reference

Andrews, Patricia L. 2018. The Rothermel surface fire spread model and associated
developments: A comprehensive explanation. Gen. Tech. Rep. RMRS-GTR-371. Fort Collins,
CO: U.S. Department of Agriculture, Forest Service, Rocky Mountain Research Station.
Table 26, pages 86-89.

Catchpole, E. A.; deMestre, N. J.; Gill, A. M. 1982. Intensity of fire at its perimeter.
Australian Forest Research. 12: 47–54.

# Example

```julia
using ModelingToolkit, WildlandFire, OrdinaryDiffEqDefault

sys = EllipticalFireSpread()
compiled_sys = mtkcompile(sys)

prob = NonlinearProblem(compiled_sys, Dict(
    compiled_sys.R_H => 0.1,       # Head fire rate of spread (m/s)
    compiled_sys.Z => 2.0,         # Length-to-width ratio
    compiled_sys.γ => π/4,         # Direction 45° from head fire direction (rad)
    compiled_sys.elapsed_time => 3600.0  # Elapsed time (s)
))
sol = solve(prob)
```
"""
@component function EllipticalFireSpread(; name=:EllipticalFireSpread)
    @constants begin
        one = 1.0, [description = "Dimensionless one for unit balancing", unit = u"1"]
        two = 2.0, [description = "Dimensionless two", unit = u"1"]
    end

    # Input parameters
    @parameters begin
        R_H, [description = "Rate of spread in direction of maximum spread (head fire)", unit = u"m/s"]
        Z, [description = "Fire length-to-width ratio (dimensionless)", unit = u"1"]
        γ, [description = "Direction from ignition point relative to max spread direction", unit = u"rad"]
        elapsed_time, [description = "Elapsed time", unit = u"s"]
    end

    # Intermediate variables
    @variables begin
        e(t), [description = "Eccentricity of fire ellipse (dimensionless)", unit = u"1"]
        D_H(t), [description = "Heading spread distance", unit = u"m"]
        D_B(t), [description = "Backing spread distance", unit = u"m"]
        L(t), [description = "Fire length", unit = u"m"]
        W(t), [description = "Fire width", unit = u"m"]
    end

    # Output variables
    @variables begin
        R_γ(t), [description = "Rate of spread from ignition point in direction γ", unit = u"m/s"]
        R_B(t), [description = "Rate of spread of backing fire", unit = u"m/s"]
        D_F(t), [description = "Flanking spread distance", unit = u"m"]
        f(t), [description = "Ellipse semi-major axis", unit = u"m"]
        g(t), [description = "Distance from focus to center", unit = u"m"]
        h(t), [description = "Ellipse semi-minor axis", unit = u"m"]
    end

    eqs = [
        # Eccentricity of ellipse - Table 26, Andrews (2018)
        e ~ sqrt(Z^2 - one) / Z,                         # Eq. e

        # Rate of spread from ignition point in direction γ - Table 26, Andrews (2018)
        R_γ ~ R_H * (one - e) / (one - e * cos(γ)),      # Eq. R_γ

        # Rate of spread of backing fire - Table 26, Andrews (2018)
        R_B ~ R_H * (one - e) / (one + e),               # Eq. R_B

        # Heading and backing spread distances - Table 26, Andrews (2018)
        D_H ~ R_H * elapsed_time,                        # Eq. D_H
        D_B ~ R_B * elapsed_time,                        # Eq. D_B

        # Fire length and width - Table 26, Andrews (2018)
        L ~ D_H + D_B,                                   # Eq. L
        W ~ L / Z,                                       # Eq. W

        # Flanking spread distance - Table 26, Andrews (2018)
        D_F ~ W / two,                                   # Eq. D_F

        # Ellipse dimensions for fire perimeter spread calculations - Table 26, Andrews (2018)
        f ~ L / two,                                     # Eq. f (semi-major axis)
        g ~ D_H - f,                                     # Eq. g (focus to center)
        h ~ D_F,                                         # Eq. h (semi-minor axis = D_F)
    ]

    return System(eqs, t; name)
end


"""
    FirePerimeterSpread(; name=:FirePerimeterSpread)

Create a system for calculating fire spread rate normal to the fire perimeter.

Fire spread rate and intensity on the perimeter of a fire are described by
Catchpole et al. (1982). Given the direction from the ignition point (γ), this
system calculates the associated direction normal to the perimeter (ψ) and the
rate of spread in that direction (R_ψ).

# Key Outputs
- `R_ψ`: Rate of spread normal to the fire perimeter
- `ψ`: Direction normal to the perimeter
- `θ`: Intermediate angle from center of ellipse
- `I_B`: Fireline intensity (Byram's)

# Reference

Andrews, Patricia L. 2018. The Rothermel surface fire spread model and associated
developments: A comprehensive explanation. Gen. Tech. Rep. RMRS-GTR-371. Fort Collins,
CO: U.S. Department of Agriculture, Forest Service, Rocky Mountain Research Station.
Table 26, pages 89-91.

Catchpole, E. A.; deMestre, N. J.; Gill, A. M. 1982. Intensity of fire at its perimeter.
Australian Forest Research. 12: 47–54.
"""
@component function FirePerimeterSpread(; name=:FirePerimeterSpread)
    @constants begin
        one = 1.0, [description = "Dimensionless one", unit = u"1"]
    end

    # Input parameters - ellipse dimensions and fire behavior
    @parameters begin
        R_H, [description = "Rate of spread in direction of maximum spread (head fire)", unit = u"m/s"]
        f, [description = "Ellipse semi-major axis", unit = u"m"]
        g, [description = "Distance from focus to center", unit = u"m"]
        h, [description = "Ellipse semi-minor axis", unit = u"m"]
        γ, [description = "Direction from ignition point relative to max spread direction", unit = u"rad"]
        H_A, [description = "Heat per unit area", unit = u"J/m^2"]
    end

    # Intermediate variables
    @variables begin
        cos_θ(t), [description = "Cosine of intermediate angle θ (dimensionless)", unit = u"1"]
        sin_θ(t), [description = "Sine of intermediate angle θ (dimensionless)", unit = u"1"]
        θ(t), [description = "Intermediate angle from center of ellipse", unit = u"rad"]
    end

    # Output variables
    @variables begin
        R_ψ(t), [description = "Rate of spread normal to fire perimeter", unit = u"m/s"]
        ψ(t), [description = "Direction normal to perimeter relative to max spread", unit = u"rad"]
        I_B(t), [description = "Fireline intensity (Byram)", unit = u"W/m"]
    end

    eqs = [
        # Intermediate angle θ - Eq. [5] in Catchpole et al. (1982), Table 26 Andrews (2018)
        # cos(θ) = [h·cos(γ)·(h²cos²γ + (f² - g²)sin²γ)^0.5 - f·g·sin²γ] / (h²cos²γ + f²sin²γ)
        cos_θ ~ (h * cos(γ) * sqrt(h^2 * cos(γ)^2 + (f^2 - g^2) * sin(γ)^2)
                 - f * g * sin(γ)^2) / (h^2 * cos(γ)^2 + f^2 * sin(γ)^2),

        # Calculate sin_θ from cos_θ
        sin_θ ~ sqrt(one - cos_θ^2),

        # θ from cos_θ
        θ ~ acos(cos_θ),

        # Rate of spread normal to perimeter - Eq. [4] in Catchpole et al. (1982)
        # R_ψ = R_H·(g·cos(θ) + f) / sqrt((f/h)²sin²θ + cos²θ) / f
        # Simplified: R_ψ = R_H·h·(g·cos(θ) + f) / (f·sqrt(h²cos²θ + f²sin²θ))
        R_ψ ~ R_H * h * (g * cos_θ + f) / (f * sqrt(h^2 * cos_θ^2 + f^2 * sin_θ^2)),

        # Direction ψ normal to perimeter - Table 26, Andrews (2018)
        # tan(ψ) = f/h · tan(θ)
        ψ ~ atan(f / h * sin_θ / cos_θ),

        # Fireline intensity - Catchpole et al. (1982), Table 26 Andrews (2018)
        # Uses R_ψ, not R_γ
        I_B ~ H_A * R_ψ,
    ]

    return System(eqs, t; name)
end


export FireSpreadDirection, EllipticalFireSpread, FirePerimeterSpread
