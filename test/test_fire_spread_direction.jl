using WildlandFire
using Test
using ModelingToolkit
using OrdinaryDiffEqDefault
using DynamicQuantities

# ==============================================================================
# Test Setup Snippet
# ==============================================================================

@testsnippet FireSpreadDirectionSetup begin
    using WildlandFire
    using Test
    using ModelingToolkit
    using OrdinaryDiffEqDefault
    using DynamicQuantities
    using DynamicQuantities: ustrip
end

# ==============================================================================
# FireSpreadDirection Tests
# ==============================================================================

@testitem "FireSpreadDirection - Structural Verification" setup = [FireSpreadDirectionSetup] tags = [:fire_spread_direction] begin
    sys = FireSpreadDirection()

    # Verify system structure
    @test sys isa System

    # Check that all expected variables exist
    var_names = [string(Symbolics.tosymbol(v, escape = false)) for v in unknowns(sys)]
    expected_vars = ["φ_combined", "R_H", "α", "φ_E", "U_E", "Z"]
    for v in expected_vars
        @test v in var_names
    end

    # Check that all expected parameters exist
    param_names = [string(Symbolics.tosymbol(p, escape = false)) for p in parameters(sys)]
    expected_params = ["R0", "φw", "φs", "ω", "β_ratio", "C_coeff", "B_coeff", "E_coeff"]
    for p in expected_params
        @test p in param_names
    end

    # Verify equation count (6 equations for 6 unknowns)
    @test length(equations(sys)) == 6
end

@testitem "FireSpreadDirection - Wind Aligned with Slope" setup = [FireSpreadDirectionSetup] tags = [:fire_spread_direction] begin
    # When wind direction ω = 0 (aligned with slope), the resultant should be
    # the simple sum of slope and wind effects

    sys = FireSpreadDirection()
    compiled_sys = mtkcompile(sys)

    # Test parameters
    R0_val = 0.01  # m/s - no-wind no-slope rate of spread
    φw_val = 5.0   # wind factor
    φs_val = 2.0   # slope factor
    ω_val = 0.0    # wind aligned with slope (upslope)
    β_ratio_val = 0.5
    C_val = 7.47
    B_val = 0.15566
    E_val = 0.715

    prob = NonlinearProblem(
        compiled_sys, Dict(
            compiled_sys.R0 => R0_val,
            compiled_sys.φw => φw_val,
            compiled_sys.φs => φs_val,
            compiled_sys.ω => ω_val,
            compiled_sys.β_ratio => β_ratio_val,
            compiled_sys.C_coeff => C_val,
            compiled_sys.B_coeff => B_val,
            compiled_sys.E_coeff => E_val,
        )
    )
    sol = solve(prob)

    # Direction of max spread should be 0 (aligned with slope)
    @test ustrip(sol[compiled_sys.α]) ≈ 0.0 atol = 1.0e-10

    # When ω = 0: φ_combined = sqrt((φs + φw)² + 0²) = φs + φw
    expected_φ_combined = φw_val + φs_val
    @test ustrip(sol[compiled_sys.φ_combined]) ≈ expected_φ_combined rtol = 1.0e-6

    # R_H = R0 * (1 + φ_combined) = R0 * (1 + φw + φs)
    expected_R_H = R0_val * (1 + expected_φ_combined)
    @test ustrip(sol[compiled_sys.R_H]) ≈ expected_R_H rtol = 1.0e-6

    # φ_E = φ_combined = (φw + φs)
    expected_φ_E = φw_val + φs_val
    @test ustrip(sol[compiled_sys.φ_E]) ≈ expected_φ_E rtol = 1.0e-6
end

@testitem "FireSpreadDirection - Wind Perpendicular to Slope" setup = [FireSpreadDirectionSetup] tags = [:fire_spread_direction] begin
    # When wind direction ω = π/2 (perpendicular to slope)

    sys = FireSpreadDirection()
    compiled_sys = mtkcompile(sys)

    R0_val = 0.01
    φw_val = 5.0
    φs_val = 2.0
    ω_val = π / 2  # perpendicular
    β_ratio_val = 0.5
    C_val = 7.47
    B_val = 0.15566
    E_val = 0.715

    prob = NonlinearProblem(
        compiled_sys, Dict(
            compiled_sys.R0 => R0_val,
            compiled_sys.φw => φw_val,
            compiled_sys.φs => φs_val,
            compiled_sys.ω => ω_val,
            compiled_sys.β_ratio => β_ratio_val,
            compiled_sys.C_coeff => C_val,
            compiled_sys.B_coeff => B_val,
            compiled_sys.E_coeff => E_val,
        )
    )
    sol = solve(prob)

    # When ω = π/2:
    # φ_combined = sqrt(φs² + φw²)
    expected_φ_combined = sqrt(φs_val^2 + φw_val^2)
    @test ustrip(sol[compiled_sys.φ_combined]) ≈ expected_φ_combined rtol = 1.0e-6

    # α = asin(|φw * sin(ω)| / φ_combined) = asin(φw / sqrt(φs² + φw²))
    expected_α = asin(φw_val / expected_φ_combined)
    @test ustrip(sol[compiled_sys.α]) ≈ expected_α rtol = 1.0e-6

    # R_H = R0 * (1 + φ_combined)
    expected_R_H = R0_val * (1 + expected_φ_combined)
    @test ustrip(sol[compiled_sys.R_H]) ≈ expected_R_H rtol = 1.0e-6
end

@testitem "FireSpreadDirection - Length-to-Width Ratio" setup = [FireSpreadDirectionSetup] tags = [:fire_spread_direction] begin
    # Test that Z ≥ 1 always (ellipse cannot be narrower than wide)

    sys = FireSpreadDirection()
    compiled_sys = mtkcompile(sys)

    # Low wind scenario
    prob = NonlinearProblem(
        compiled_sys, Dict(
            compiled_sys.R0 => 0.01,
            compiled_sys.φw => 0.1,  # very low wind
            compiled_sys.φs => 0.1,  # very low slope
            compiled_sys.ω => 0.0,
            compiled_sys.β_ratio => 0.5,
            compiled_sys.C_coeff => 7.47,
            compiled_sys.B_coeff => 0.15566,
            compiled_sys.E_coeff => 0.715,
        )
    )
    sol = solve(prob)

    @test ustrip(sol[compiled_sys.Z]) ≥ 1.0
    @test ustrip(sol[compiled_sys.U_E]) ≥ 0.0
end

@testitem "FireSpreadDirection - Realistic U_E and Z values" setup = [FireSpreadDirectionSetup] tags = [:fire_spread_direction] begin
    # Regression test for issue #50: U_ref must be 0.3048/60 m/s (1 ft/min)
    # With typical fire conditions, U_E and Z should be physically realistic.
    # Before the fix, U_E was ~197x too large and Z was ~496 instead of ~3.5.

    sys = FireSpreadDirection()
    compiled_sys = mtkcompile(sys)

    # Typical conditions: moderate wind factor, no slope, wind aligned upslope
    prob = NonlinearProblem(
        compiled_sys, Dict(
            compiled_sys.R0 => 0.01,
            compiled_sys.φw => 5.0,
            compiled_sys.φs => 0.0,
            compiled_sys.ω => 0.0,
            compiled_sys.β_ratio => 0.5,
            compiled_sys.C_coeff => 7.47,
            compiled_sys.B_coeff => 0.15566,
            compiled_sys.E_coeff => 0.715,
        )
    )
    sol = solve(prob)

    U_E_val = ustrip(sol[compiled_sys.U_E])
    Z_val = ustrip(sol[compiled_sys.Z])

    # U_E should be on the order of a few m/s, not hundreds
    @test U_E_val < 50.0  # m/s; before fix this was ~885 m/s
    @test U_E_val > 0.0

    # Z should be a reasonable fire length-to-width ratio (typically 1-10)
    @test Z_val < 30.0  # before fix this was ~496
    @test Z_val > 1.0

    # Verify U_E is consistent with the Rothermel model's U_ref = 0.3048/60
    # U_E = U_ref * (φ_E * β_ratio^E / C)^(1/B)
    U_ref = 0.3048 / 60  # must match the implementation
    φ_E_val = ustrip(sol[compiled_sys.φ_E])
    expected_U_E = U_ref * (φ_E_val * 0.5^0.715 / 7.47)^(1.0 / 0.15566)
    @test U_E_val ≈ expected_U_E rtol = 1.0e-6
end

# ==============================================================================
# EllipticalFireSpread Tests
# ==============================================================================

@testitem "EllipticalFireSpread - Structural Verification" setup = [FireSpreadDirectionSetup] tags = [:fire_spread_direction] begin
    sys = EllipticalFireSpread()

    @test sys isa System

    var_names = [string(Symbolics.tosymbol(v, escape = false)) for v in unknowns(sys)]
    expected_vars = ["e", "D_H", "D_B", "L", "W", "R_γ", "R_B", "D_F", "f", "g", "h"]
    for v in expected_vars
        @test v in var_names
    end

    param_names = [string(Symbolics.tosymbol(p, escape = false)) for p in parameters(sys)]
    expected_params = ["R_H", "Z", "γ", "elapsed_time"]
    for p in expected_params
        @test p in param_names
    end
end

@testitem "EllipticalFireSpread - Eccentricity Bounds" setup = [FireSpreadDirectionSetup] tags = [:fire_spread_direction] begin
    # For Z > 1, eccentricity should be 0 < e < 1

    sys = EllipticalFireSpread()
    compiled_sys = mtkcompile(sys)

    for Z_val in [1.5, 2.0, 3.0, 5.0]
        prob = NonlinearProblem(
            compiled_sys, Dict(
                compiled_sys.R_H => 0.1,
                compiled_sys.Z => Z_val,
                compiled_sys.γ => 0.0,
                compiled_sys.elapsed_time => 3600.0,
            )
        )
        sol = solve(prob)

        # e = sqrt(Z² - 1) / Z
        expected_e = sqrt(Z_val^2 - 1) / Z_val
        @test ustrip(sol[compiled_sys.e]) ≈ expected_e rtol = 1.0e-6
        @test 0.0 < ustrip(sol[compiled_sys.e]) < 1.0
    end
end

@testitem "EllipticalFireSpread - Head Fire Direction" setup = [FireSpreadDirectionSetup] tags = [:fire_spread_direction] begin
    # When γ = 0 (head fire direction), R_γ should equal R_H

    sys = EllipticalFireSpread()
    compiled_sys = mtkcompile(sys)

    R_H_val = 0.1
    Z_val = 2.0

    prob = NonlinearProblem(
        compiled_sys, Dict(
            compiled_sys.R_H => R_H_val,
            compiled_sys.Z => Z_val,
            compiled_sys.γ => 0.0,  # head fire direction
            compiled_sys.elapsed_time => 3600.0,
        )
    )
    sol = solve(prob)

    # At γ = 0: R_γ = R_H*(1-e)/(1-e*cos(0)) = R_H*(1-e)/(1-e) = R_H
    @test ustrip(sol[compiled_sys.R_γ]) ≈ R_H_val rtol = 1.0e-6
end

@testitem "EllipticalFireSpread - Back Fire Direction" setup = [FireSpreadDirectionSetup] tags = [:fire_spread_direction] begin
    # When γ = π (back fire direction), R_γ should equal R_B

    sys = EllipticalFireSpread()
    compiled_sys = mtkcompile(sys)

    R_H_val = 0.1
    Z_val = 2.0

    prob = NonlinearProblem(
        compiled_sys, Dict(
            compiled_sys.R_H => R_H_val,
            compiled_sys.Z => Z_val,
            compiled_sys.γ => π,  # back fire direction
            compiled_sys.elapsed_time => 3600.0,
        )
    )
    sol = solve(prob)

    # At γ = π: R_γ = R_H*(1-e)/(1-e*cos(π)) = R_H*(1-e)/(1+e) = R_B
    @test ustrip(sol[compiled_sys.R_γ]) ≈ ustrip(sol[compiled_sys.R_B]) rtol = 1.0e-6
end

@testitem "EllipticalFireSpread - Fire Dimensions" setup = [FireSpreadDirectionSetup] tags = [:fire_spread_direction] begin
    sys = EllipticalFireSpread()
    compiled_sys = mtkcompile(sys)

    R_H_val = 0.1  # m/s
    Z_val = 2.0
    t_val = 3600.0  # 1 hour

    prob = NonlinearProblem(
        compiled_sys, Dict(
            compiled_sys.R_H => R_H_val,
            compiled_sys.Z => Z_val,
            compiled_sys.γ => 0.0,
            compiled_sys.elapsed_time => t_val,
        )
    )
    sol = solve(prob)

    # D_H = R_H * t
    expected_D_H = R_H_val * t_val
    @test ustrip(sol[compiled_sys.D_H]) ≈ expected_D_H rtol = 1.0e-6

    # L = D_H + D_B
    @test ustrip(sol[compiled_sys.L]) ≈ ustrip(sol[compiled_sys.D_H]) + ustrip(sol[compiled_sys.D_B]) rtol = 1.0e-6

    # W = L / Z
    @test ustrip(sol[compiled_sys.W]) ≈ ustrip(sol[compiled_sys.L]) / Z_val rtol = 1.0e-6

    # D_F = W / 2
    @test ustrip(sol[compiled_sys.D_F]) ≈ ustrip(sol[compiled_sys.W]) / 2.0 rtol = 1.0e-6

    # Ellipse dimensions
    @test ustrip(sol[compiled_sys.f]) ≈ ustrip(sol[compiled_sys.L]) / 2.0 rtol = 1.0e-6
    @test ustrip(sol[compiled_sys.h]) ≈ ustrip(sol[compiled_sys.D_F]) rtol = 1.0e-6
end

# ==============================================================================
# FirePerimeterSpread Tests
# ==============================================================================

@testitem "FirePerimeterSpread - Structural Verification" setup = [FireSpreadDirectionSetup] tags = [:fire_spread_direction] begin
    sys = FirePerimeterSpread()

    @test sys isa System

    var_names = [string(Symbolics.tosymbol(v, escape = false)) for v in unknowns(sys)]
    expected_vars = ["cos_θ", "sin_θ", "θ", "R_ψ", "ψ", "I_B"]
    for v in expected_vars
        @test v in var_names
    end

    param_names = [string(Symbolics.tosymbol(p, escape = false)) for p in parameters(sys)]
    expected_params = ["R_H", "f", "g", "h", "γ", "H_A"]
    for p in expected_params
        @test p in param_names
    end
end

@testitem "FirePerimeterSpread - Head Fire Perimeter" setup = [FireSpreadDirectionSetup] tags = [:fire_spread_direction] begin
    # At γ = 0 (head fire direction), θ should be 0 and R_ψ can be calculated from the equation

    sys = FirePerimeterSpread()
    compiled_sys = mtkcompile(sys)

    R_H_val = 0.1
    # Use ellipse dimensions consistent with fire geometry
    # For a fire ellipse: f = L/2, g = D_H - f, h = W/2
    # Example: Z = 2 (length-to-width ratio)
    # If L = 360, W = 180, then f = 180, h = 90
    # D_H = R_H * t, D_B = R_B * t, L = D_H + D_B
    # For Z = 2, e = sqrt(3)/2 ≈ 0.866
    # D_H/(D_H + D_B) = (1+e)/(2) ≈ 0.933
    f_val = 180.0   # semi-major axis
    h_val = 90.0    # semi-minor axis
    g_val = sqrt(f_val^2 - h_val^2)  # focal distance from center
    H_A_val = 1000.0  # J/m²

    prob = NonlinearProblem(
        compiled_sys, Dict(
            compiled_sys.R_H => R_H_val,
            compiled_sys.f => f_val,
            compiled_sys.g => g_val,
            compiled_sys.h => h_val,
            compiled_sys.γ => 0.0,
            compiled_sys.H_A => H_A_val,
        )
    )
    sol = solve(prob)

    # At head fire (γ = 0), θ should be 0
    @test ustrip(sol[compiled_sys.θ]) ≈ 0.0 atol = 1.0e-6

    # Verify R_ψ follows the equation at θ = 0:
    # R_ψ = R_H * h * (g*cos(θ) + f) / (f * sqrt(h²cos²θ + f²sin²θ))
    # At θ = 0: R_ψ = R_H * h * (g + f) / (f * h) = R_H * (g + f) / f
    expected_R_ψ = R_H_val * (g_val + f_val) / f_val
    @test ustrip(sol[compiled_sys.R_ψ]) ≈ expected_R_ψ rtol = 1.0e-3

    # R_ψ should be positive
    @test sol[compiled_sys.R_ψ] > 0
end

@testitem "FirePerimeterSpread - Fireline Intensity" setup = [FireSpreadDirectionSetup] tags = [:fire_spread_direction] begin
    sys = FirePerimeterSpread()
    compiled_sys = mtkcompile(sys)

    R_H_val = 0.1
    f_val = 180.0
    h_val = 90.0
    g_val = sqrt(f_val^2 - h_val^2)  # proper focal distance
    H_A_val = 5000.0  # J/m²

    prob = NonlinearProblem(
        compiled_sys, Dict(
            compiled_sys.R_H => R_H_val,
            compiled_sys.f => f_val,
            compiled_sys.g => g_val,
            compiled_sys.h => h_val,
            compiled_sys.γ => 0.0,
            compiled_sys.H_A => H_A_val,
        )
    )
    sol = solve(prob)

    # I_B = H_A * R_ψ - verify the relationship holds
    @test ustrip(sol[compiled_sys.I_B]) ≈ H_A_val * ustrip(sol[compiled_sys.R_ψ]) rtol = 1.0e-6

    # Fireline intensity should be positive
    @test ustrip(sol[compiled_sys.I_B]) > 0.0
end

@testitem "FirePerimeterSpread - Backing Fire Perimeter" setup = [FireSpreadDirectionSetup] tags = [:fire_spread_direction] begin
    # At γ = π (back fire), R_ψ should be less than at head fire

    sys = FirePerimeterSpread()
    compiled_sys = mtkcompile(sys)

    R_H_val = 0.1
    f_val = 180.0
    h_val = 90.0
    g_val = sqrt(f_val^2 - h_val^2)  # proper focal distance

    # Head fire case
    prob_head = NonlinearProblem(
        compiled_sys, Dict(
            compiled_sys.R_H => R_H_val,
            compiled_sys.f => f_val,
            compiled_sys.g => g_val,
            compiled_sys.h => h_val,
            compiled_sys.γ => 0.0,
            compiled_sys.H_A => 1000.0,
        )
    )
    sol_head = solve(prob_head)

    # Back fire case
    prob_back = NonlinearProblem(
        compiled_sys, Dict(
            compiled_sys.R_H => R_H_val,
            compiled_sys.f => f_val,
            compiled_sys.g => g_val,
            compiled_sys.h => h_val,
            compiled_sys.γ => π,
            compiled_sys.H_A => 1000.0,
        )
    )
    sol_back = solve(prob_back)

    # Back fire spread rate should be less than head fire spread rate
    @test ustrip(sol_back[compiled_sys.R_ψ]) < ustrip(sol_head[compiled_sys.R_ψ])
end

# ==============================================================================
# Integration Tests - Combined System Behavior
# ==============================================================================

@testitem "Integration - Fire Spread Direction to Elliptical Shape" setup = [FireSpreadDirectionSetup] tags = [:fire_spread_direction] begin
    # Test that outputs from FireSpreadDirection can be used as inputs to EllipticalFireSpread

    dir_sys = FireSpreadDirection()
    compiled_dir = mtkcompile(dir_sys)

    prob_dir = NonlinearProblem(
        compiled_dir, Dict(
            compiled_dir.R0 => 0.01,
            compiled_dir.φw => 5.0,
            compiled_dir.φs => 2.0,
            compiled_dir.ω => π / 4,
            compiled_dir.β_ratio => 0.5,
            compiled_dir.C_coeff => 7.47,
            compiled_dir.B_coeff => 0.15566,
            compiled_dir.E_coeff => 0.715,
        )
    )
    sol_dir = solve(prob_dir)

    # Use outputs as inputs to EllipticalFireSpread
    ellipse_sys = EllipticalFireSpread()
    compiled_ellipse = mtkcompile(ellipse_sys)

    prob_ellipse = NonlinearProblem(
        compiled_ellipse, Dict(
            compiled_ellipse.R_H => ustrip(sol_dir[compiled_dir.R_H]),
            compiled_ellipse.Z => ustrip(sol_dir[compiled_dir.Z]),
            compiled_ellipse.γ => π / 4,  # 45 degrees from head fire direction
            compiled_ellipse.elapsed_time => 3600.0,
        )
    )
    sol_ellipse = solve(prob_ellipse)

    # Verify reasonable outputs
    @test ustrip(sol_ellipse[compiled_ellipse.L]) > 0
    @test ustrip(sol_ellipse[compiled_ellipse.W]) > 0
    @test ustrip(sol_ellipse[compiled_ellipse.L]) > ustrip(sol_ellipse[compiled_ellipse.W])  # ellipse is longer than wide
end
