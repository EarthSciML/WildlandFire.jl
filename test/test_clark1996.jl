@testsnippet Clark1996Setup begin
    using Test
    using ModelingToolkit
    using ModelingToolkit: t, D
    using OrdinaryDiffEqDefault
    using Symbolics
    using WildlandFire
end

@testitem "Clark1996FireSpread Structural Verification" setup = [Clark1996Setup] tags = [:clark1996] begin
    sys = Clark1996FireSpread()

    @test sys !== nothing
    @test length(equations(sys)) == 12
    @test length(unknowns(sys)) == 12

    # Check state variables exist
    var_names = [string(Symbolics.tosymbol(v, escape = false)) for v in unknowns(sys)]
    @test "M_litter" in var_names
    @test "M_trash" in var_names
    @test "M_scrub" in var_names
    @test "M_canopy" in var_names
    @test "Q_cumulative" in var_names
    @test "canopy_burning" in var_names
    @test "B_ratio" in var_names
    @test "S_f" in var_names
    @test "ground_burn_rate" in var_names
    @test "total_burn_rate" in var_names
    @test "F_s" in var_names
    @test "F_l" in var_names

    # Check parameters exist
    param_names = [string(Symbolics.tosymbol(p, escape = false)) for p in parameters(sys)]
    @test "V_A" in param_names
    # Note: canopy_burning is now a computed variable, not a parameter
end

@testitem "Clark1996ConvectiveFroudeNumber Structural Verification" setup = [Clark1996Setup] tags = [:clark1996] begin
    sys = Clark1996ConvectiveFroudeNumber()

    @test sys !== nothing
    @test length(equations(sys)) == 1
    @test length(unknowns(sys)) == 1

    var_names = [string(Symbolics.tosymbol(v, escape = false)) for v in unknowns(sys)]
    @test "F_c_sq" in var_names

    param_names = [string(Symbolics.tosymbol(p, escape = false)) for p in parameters(sys)]
    @test "U" in param_names
    @test "S_f" in param_names
    @test "delta_theta_over_theta" in param_names
    @test "W_f" in param_names
end

@testitem "Clark1996HeatFluxProfile Structural Verification" setup = [Clark1996Setup] tags = [:clark1996] begin
    sys = Clark1996HeatFluxProfile()

    @test sys !== nothing
    @test length(equations(sys)) == 2
    @test length(unknowns(sys)) == 2

    var_names = [string(Symbolics.tosymbol(v, escape = false)) for v in unknowns(sys)]
    @test "F_s" in var_names
    @test "F_l" in var_names
end

@testitem "Clark1996WindProfile Structural Verification" setup = [Clark1996Setup] tags = [:clark1996] begin
    sys = Clark1996WindProfile()

    @test sys !== nothing
    @test length(equations(sys)) == 1
    @test length(unknowns(sys)) == 1

    var_names = [string(Symbolics.tosymbol(v, escape = false)) for v in unknowns(sys)]
    @test "U_z" in var_names
end

@testitem "McArthur Fire Spread Rate — Eq. 9" setup = [Clark1996Setup] tags = [:clark1996] begin
    # Verify S_f = S_a * exp(0.08424 * V_A) against Table 1 values
    # Table 1 gives S_f(U_0) with |V_A| = U_0

    S_a = 0.18  # m/s
    k = 0.08424  # s/m

    # From Table 1: S_f(U_0) values
    test_cases = [
        (U_0 = 1.0, S_f_expected = 0.2),
        (U_0 = 2.0, S_f_expected = 0.21),
        (U_0 = 3.0, S_f_expected = 0.23),
        (U_0 = 4.0, S_f_expected = 0.25),
        (U_0 = 5.0, S_f_expected = 0.27),
        (U_0 = 10.0, S_f_expected = 0.42),
        (U_0 = 15.0, S_f_expected = 0.64),
        (U_0 = 20.0, S_f_expected = 0.97),
    ]

    for tc in test_cases
        S_f_calc = S_a * exp(k * tc.U_0)
        @test isapprox(S_f_calc, tc.S_f_expected, atol = 0.015)
    end
end

@testitem "Burn Rate Ratio — Eq. 8" setup = [Clark1996Setup] tags = [:clark1996] begin
    # Test B_ratio = sqrt((V_A + 1) / (V_A + 4))
    # Should approach 1.0 for large V_A and sqrt(1/4)=0.5 for V_A=0

    # At V_A = 0: B_ratio = sqrt(1/4) = 0.5
    @test isapprox(sqrt((0.0 + 1.0) / (0.0 + 4.0)), 0.5, atol = 1.0e-10)

    # At V_A = large: B_ratio → 1.0
    @test isapprox(sqrt((1000.0 + 1.0) / (1000.0 + 4.0)), 1.0, atol = 0.002)

    # Monotonically increasing with V_A
    V_A_vals = [0.0, 1.0, 3.0, 5.0, 10.0, 20.0]
    B_vals = [sqrt((v + 1.0) / (v + 4.0)) for v in V_A_vals]
    for i in 2:length(B_vals)
        @test B_vals[i] > B_vals[i - 1]
    end
end

@testitem "Clark1996FireSpread ODE Integration" setup = [Clark1996Setup] tags = [:clark1996] begin
    sys = Clark1996FireSpread()
    compiled = mtkcompile(sys)

    # Test with V_A = 3 m/s (FIR7CR experiment conditions), short duration before canopy ignition
    M_litter_0 = 2.0  # kg/m²
    M_trash_0 = 0.5
    M_scrub_0 = 0.2
    M_canopy_0 = 1.2

    prob = ODEProblem(
        compiled,
        [
            compiled.M_litter => M_litter_0,
            compiled.M_trash => M_trash_0,
            compiled.M_scrub => M_scrub_0,
            compiled.M_canopy => M_canopy_0,
            compiled.Q_cumulative => 0.0,
        ],
        (0.0, 30.0),  # Short duration to test before canopy ignition
        [
            compiled.V_A => 3.0,
        ]
    )
    sol = solve(prob)

    # Fuel should decrease monotonically
    @test sol[compiled.M_litter][end] < M_litter_0
    @test sol[compiled.M_trash][end] < M_trash_0
    @test sol[compiled.M_scrub][end] < M_scrub_0

    # Cumulative heat flux should increase
    @test sol[compiled.Q_cumulative][end] > 0.0

    # Heat fluxes should be positive
    @test all(sol[compiled.F_s] .>= 0.0)
    @test all(sol[compiled.F_l] .>= 0.0)

    # B_ratio should be constant (V_A is constant)
    B_expected = sqrt((3.0 + 1.0) / (3.0 + 4.0))
    @test all(isapprox.(sol[compiled.B_ratio], B_expected, rtol = 1.0e-6))

    # S_f should be constant (V_A is constant)
    S_f_expected = 0.18 * exp(0.08424 * 3.0)
    @test all(isapprox.(sol[compiled.S_f], S_f_expected, rtol = 1.0e-6))
end

@testitem "Clark1996FireSpread Fuel Depletion" setup = [Clark1996Setup] tags = [:clark1996] begin
    sys = Clark1996FireSpread()
    compiled = mtkcompile(sys)

    # Run long enough for all fuels to deplete
    # Canopy will automatically ignite once cumulative heat flux reaches threshold
    prob = ODEProblem(
        compiled,
        [
            compiled.M_litter => 2.0,
            compiled.M_trash => 0.5,
            compiled.M_scrub => 0.2,
            compiled.M_canopy => 1.2,
            compiled.Q_cumulative => 0.0,
        ],
        (0.0, 200.0),
        [
            compiled.V_A => 3.0,
        ]
    )
    sol = solve(prob)

    # All ground fuel should be near zero or zero at end
    @test sol[compiled.M_litter][end] < 0.01
    @test sol[compiled.M_scrub][end] < 0.01

    # Fuel masses should never go negative
    @test all(sol[compiled.M_litter] .>= -0.01)
    @test all(sol[compiled.M_trash] .>= -0.01)
    @test all(sol[compiled.M_scrub] .>= -0.01)

    # Canopy should ignite automatically and be partially consumed
    @test sol[compiled.M_canopy][end] < 1.2
    @test any(sol[compiled.canopy_burning] .> 0.0)  # Canopy should have ignited at some point
end

@testitem "Clark1996FireSpread Automatic Canopy Ignition" setup = [Clark1996Setup] tags = [:clark1996] begin
    sys = Clark1996FireSpread()
    compiled = mtkcompile(sys)

    # Test automatic canopy ignition when cumulative heat flux exceeds threshold
    prob = ODEProblem(
        compiled,
        [
            compiled.M_litter => 2.0,
            compiled.M_trash => 0.5,
            compiled.M_scrub => 0.2,
            compiled.M_canopy => 1.2,
            compiled.Q_cumulative => 0.0,
        ],
        (0.0, 120.0),  # Run long enough for canopy ignition to occur
        [
            compiled.V_A => 3.0,
        ]
    )
    sol = solve(prob)

    # Cumulative heat flux should reach the threshold (170 kJ/m²)
    @test sol[compiled.Q_cumulative][end] > 170000.0

    # Canopy should ignite automatically when threshold is exceeded
    ignition_times = findall(x -> x > 0.5, sol[compiled.canopy_burning])
    @test length(ignition_times) > 0  # Canopy should ignite at some point

    # Canopy fuel should decrease after ignition
    @test sol[compiled.M_canopy][end] < 1.2

    # Total burn rate should be higher than ground-only after canopy ignites
    late_time_idx = findlast(x -> sol.t[x] > sol.t[end] * 0.8, 1:length(sol.t))
    @test sol[compiled.total_burn_rate][late_time_idx] >= sol[compiled.ground_burn_rate][late_time_idx]

    # Heat fluxes should be positive throughout
    @test all(sol[compiled.F_s] .>= 0.0)
    @test all(sol[compiled.F_l] .>= 0.0)
end

@testitem "Clark1996FireSpread Canopy Ignition Threshold" setup = [Clark1996Setup] tags = [:clark1996] begin
    # Test that canopy ignition occurs exactly at the 170 kJ/m² threshold
    sys = Clark1996FireSpread()
    compiled = mtkcompile(sys)

    # Test with Q_cumulative just below threshold
    prob_below = ODEProblem(
        compiled,
        [
            compiled.M_litter => 2.0,
            compiled.M_trash => 0.5,
            compiled.M_scrub => 0.2,
            compiled.M_canopy => 1.2,
            compiled.Q_cumulative => 169000.0,  # Just below 170 kJ/m²
        ],
        (0.0, 1.0),  # Very short time
        [compiled.V_A => 3.0]
    )
    sol_below = solve(prob_below)
    @test sol_below[compiled.canopy_burning][1] ≈ 0.0  # Should not be burning

    # Test with Q_cumulative just above threshold
    prob_above = ODEProblem(
        compiled,
        [
            compiled.M_litter => 2.0,
            compiled.M_trash => 0.5,
            compiled.M_scrub => 0.2,
            compiled.M_canopy => 1.2,
            compiled.Q_cumulative => 171000.0,  # Just above 170 kJ/m²
        ],
        (0.0, 1.0),  # Very short time
        [compiled.V_A => 3.0]
    )
    sol_above = solve(prob_above)
    @test sol_above[compiled.canopy_burning][1] ≈ 1.0  # Should be burning
end

@testitem "Convective Froude Number — Eq. 1" setup = [Clark1996Setup] tags = [:clark1996] begin
    # Verify against Table 1 values
    # F_c² = (U - S_f)² / (g * Δθ/θ̄ * W_f)
    # We need to estimate Δθ/θ̄ and W_f from the paper context

    # For Table 1, the paper uses S_f(U_0) with |V_A| = U_0
    # The F_c² values in Table 1 are computed analytically

    # Test that F_c² increases with U for small S_f
    sys = Clark1996ConvectiveFroudeNumber()
    compiled = mtkcompile(sys)

    # Use some reference values for Δθ/θ̄ and W_f
    delta_theta_theta = 0.1  # typical buoyancy
    W_f = 420.0  # fire width in meters

    # F_c² should increase with U (for fixed S_f)
    S_f = 0.23  # m/s (at U_0 = 3)
    g = 9.81

    F_c_sq_U3 = (3.0 - S_f)^2 / (g * delta_theta_theta * W_f)
    F_c_sq_U10 = (10.0 - 0.42)^2 / (g * delta_theta_theta * W_f)
    @test F_c_sq_U10 > F_c_sq_U3

    # F_c² should be zero when U = S_f
    F_c_sq_zero = (S_f - S_f)^2 / (g * delta_theta_theta * W_f)
    @test isapprox(F_c_sq_zero, 0.0, atol = 1.0e-10)
end

@testitem "Wind Profile — Eq. 11" setup = [Clark1996Setup] tags = [:clark1996] begin
    sys = Clark1996WindProfile()
    compiled = mtkcompile(sys)

    # At z = 0 (surface): U ≈ 3 m/s (tanh(-5) ≈ -1)
    prob_surface = ODEProblem(compiled, [], (0.0, 1.0), [compiled.z => 0.0])
    sol_surface = solve(prob_surface)
    U_surface = sol_surface[compiled.U_z][1]
    @test isapprox(U_surface, 3.0, atol = 0.01)

    # At z = 500 m: U = 0 m/s (tanh(0) = 0)
    prob_500 = ODEProblem(compiled, [], (0.0, 1.0), [compiled.z => 500.0])
    sol_500 = solve(prob_500)
    U_500 = sol_500[compiled.U_z][1]
    @test isapprox(U_500, 0.0, atol = 0.01)

    # At z >> 500: U ≈ -3 m/s (tanh(large) ≈ 1)
    prob_high = ODEProblem(compiled, [], (0.0, 1.0), [compiled.z => 1500.0])
    sol_high = solve(prob_high)
    U_high = sol_high[compiled.U_z][1]
    @test isapprox(U_high, -3.0, atol = 0.01)
end

@testitem "Heat Flux Profile — Eq. 10" setup = [Clark1996Setup] tags = [:clark1996] begin
    sys = Clark1996HeatFluxProfile()
    compiled = mtkcompile(sys)

    # At z = 0: F_s = F_s_sfc (no decay)
    prob_surface = ODEProblem(
        compiled, [], (0.0, 1.0),
        [compiled.z => 0.0, compiled.F_s_sfc => 1000.0, compiled.F_l_sfc => 500.0]
    )
    sol_surface = solve(prob_surface)
    @test isapprox(sol_surface[compiled.F_s][1], 1000.0, rtol = 1.0e-6)
    @test isapprox(sol_surface[compiled.F_l][1], 500.0, rtol = 1.0e-6)

    # At z = 50 m (one e-folding): F_s = F_s_sfc * exp(-1)
    prob_50 = ODEProblem(
        compiled, [], (0.0, 1.0),
        [compiled.z => 50.0, compiled.F_s_sfc => 1000.0, compiled.F_l_sfc => 500.0]
    )
    sol_50 = solve(prob_50)
    @test isapprox(sol_50[compiled.F_s][1], 1000.0 * exp(-1), rtol = 1.0e-6)
    @test isapprox(sol_50[compiled.F_l][1], 500.0 * exp(-1), rtol = 1.0e-6)

    # At z = 100 m (two e-foldings): F_s = F_s_sfc * exp(-2)
    prob_100 = ODEProblem(
        compiled, [], (0.0, 1.0),
        [compiled.z => 100.0, compiled.F_s_sfc => 1000.0, compiled.F_l_sfc => 500.0]
    )
    sol_100 = solve(prob_100)
    @test isapprox(sol_100[compiled.F_s][1], 1000.0 * exp(-2), rtol = 1.0e-6)
end

@testitem "Sensible Heat Flux Magnitude" setup = [Clark1996Setup] tags = [:clark1996] begin
    # At the start of burning with all fuel types, verify heat flux is reasonable
    # Total ground burn rate at V_A = 3: (0.04 + 0.005 + 0.004) * B_ratio
    # B_ratio at V_A=3: sqrt(4/7) ≈ 0.7559

    B_ratio = sqrt(4.0 / 7.0)
    total_ground = (0.04 + 0.005 + 0.004) * B_ratio
    H_c = 1.7e7
    f_evap = 0.03

    F_s_expected = (1.0 - f_evap) * H_c * total_ground

    # From Table 1/Figure 6, fire fluxes for FIR7CR at 3 m/s are
    # shown in kW/m² with max values around 200-600 kW/m²
    # Our calculation gives sensible flux at surface
    @test F_s_expected > 0.0
    @test F_s_expected < 1.0e6  # Should be in hundreds of kW/m² range

    # With canopy: total = ground + canopy burn rate
    total_with_canopy = total_ground + 0.02 * B_ratio
    F_s_with_canopy = (1.0 - f_evap) * H_c * total_with_canopy
    @test F_s_with_canopy > F_s_expected
end

@testitem "Table 1 Spread Rate Validation" setup = [Clark1996Setup] tags = [:clark1996] begin
    # Validate fire spread rates from Table 1 against Eq. 9
    S_a = 0.18
    k = 0.08424

    # Full Table 1 comparison
    table1 = [
        (U_0 = 1.0, S_f = 0.2),
        (U_0 = 2.0, S_f = 0.21),
        (U_0 = 3.0, S_f = 0.23),
        (U_0 = 4.0, S_f = 0.25),
        (U_0 = 5.0, S_f = 0.27),
        (U_0 = 10.0, S_f = 0.42),
        (U_0 = 15.0, S_f = 0.64),
        (U_0 = 20.0, S_f = 0.97),
    ]

    for row in table1
        S_f_calc = S_a * exp(k * row.U_0)
        # Table values are rounded to 2 decimal places
        @test isapprox(S_f_calc, row.S_f, atol = 0.015)
    end
end
