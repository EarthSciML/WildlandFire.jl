using WildlandFire
using Test
using OrdinaryDiffEqDefault
using SciMLBase

@testset "WildlandFire.jl" begin
    @testset "Basic Rothermel Model" begin
        # Test that the model is available
        @test isdefined(WildlandFire, :rothermel_simplified)

        # Test basic fire spread calculation
        params = [
            h => 8000.0,
            S_T => 0.0555,
            S_e => 0.010,
            ρ_p => 32.0,
            σ => 3500.0,
            w_o => 0.138,
            δ => 1.0,
            M_x => 0.12,
            M_f => 0.05,
            U => 352.0,
            tan_ϕ => 0.0,
        ]

        prob = ODEProblem(rothermel_simplified, [], (0.0, 1.0), params)
        sol = solve(prob)

        # Test that solution converges
        @test SciMLBase.successful_retcode(sol)

        # Test that rate of spread is positive and reasonable
        R_val = sol[R][end]
        @test R_val > 0.0
        @test R_val < 1000.0  # ft/min (upper bound check)

        # Test that intermediate values are in expected ranges
        @test sol[β][end] > 0.0
        @test sol[β][end] < 1.0
        @test sol[I_R][end] > 0.0
        @test sol[ξ][end] > 0.0
        @test sol[ξ][end] < 1.0
    end

    @testset "Wind Effect" begin
        # Test with no wind
        params_no_wind = [
            h => 8000.0, S_T => 0.0555, S_e => 0.010, ρ_p => 32.0,
            σ => 3500.0, w_o => 0.138, δ => 1.0, M_x => 0.12,
            M_f => 0.05, U => 0.0, tan_ϕ => 0.0,
        ]

        prob = ODEProblem(rothermel_simplified, [], (0.0, 1.0), params_no_wind)
        sol_no_wind = solve(prob)
        R_no_wind = sol_no_wind[R][end]

        # Test with wind
        params_wind = [
            h => 8000.0, S_T => 0.0555, S_e => 0.010, ρ_p => 32.0,
            σ => 3500.0, w_o => 0.138, δ => 1.0, M_x => 0.12,
            M_f => 0.05, U => 880.0, tan_ϕ => 0.0,  # 10 mph
        ]

        prob = ODEProblem(rothermel_simplified, [], (0.0, 1.0), params_wind)
        sol_wind = solve(prob)
        R_wind = sol_wind[R][end]

        # Wind should increase rate of spread
        @test R_wind > R_no_wind

        # Wind factor should be positive
        @test sol_wind[ϕ_w][end] > 0.0
    end

    @testset "Moisture Effect" begin
        # Dry fuel
        params_dry = [
            h => 8000.0, S_T => 0.0555, S_e => 0.010, ρ_p => 32.0,
            σ => 3500.0, w_o => 0.138, δ => 1.0, M_x => 0.12,
            M_f => 0.01, U => 0.0, tan_ϕ => 0.0,
        ]

        prob = ODEProblem(rothermel_simplified, [], (0.0, 1.0), params_dry)
        sol_dry = solve(prob)
        R_dry = sol_dry[R][end]

        # Wet fuel
        params_wet = [
            h => 8000.0, S_T => 0.0555, S_e => 0.010, ρ_p => 32.0,
            σ => 3500.0, w_o => 0.138, δ => 1.0, M_x => 0.12,
            M_f => 0.10, U => 0.0, tan_ϕ => 0.0,
        ]

        prob = ODEProblem(rothermel_simplified, [], (0.0, 1.0), params_wet)
        sol_wet = solve(prob)
        R_wet = sol_wet[R][end]

        # Higher moisture should decrease rate of spread
        @test R_dry > R_wet

        # Moisture damping coefficient should be less than 1
        @test sol_wet[η_M][end] < 1.0
        @test sol_wet[η_M][end] >= 0.0
    end

    @testset "Slope Effect" begin
        # Flat ground
        params_flat = [
            h => 8000.0, S_T => 0.0555, S_e => 0.010, ρ_p => 32.0,
            σ => 3500.0, w_o => 0.138, δ => 1.0, M_x => 0.12,
            M_f => 0.05, U => 0.0, tan_ϕ => 0.0,
        ]

        prob = ODEProblem(rothermel_simplified, [], (0.0, 1.0), params_flat)
        sol_flat = solve(prob)
        R_flat = sol_flat[R][end]

        # Upslope
        params_slope = [
            h => 8000.0, S_T => 0.0555, S_e => 0.010, ρ_p => 32.0,
            σ => 3500.0, w_o => 0.138, δ => 1.0, M_x => 0.12,
            M_f => 0.05, U => 0.0, tan_ϕ => 0.3,  # 30% slope
        ]

        prob = ODEProblem(rothermel_simplified, [], (0.0, 1.0), params_slope)
        sol_slope = solve(prob)
        R_slope = sol_slope[R][end]

        # Slope should increase rate of spread
        @test R_slope > R_flat

        # Slope factor should be positive
        @test sol_slope[ϕ_s][end] > 0.0
    end

    @testset "Heterogeneous Fuel Model" begin
        # Test that the heterogeneous model creation function exists
        @test isdefined(WildlandFire, :create_heterogeneous_rothermel_model)

        # Note: Full heterogeneous model testing requires ModelingToolkit
        # array handling improvements. The basic model is fully functional.
        # See examples/heterogeneous_example.jl for usage.
    end

    @testset "Physical Bounds" begin
        # Test that all physical quantities stay within reasonable bounds
        params = [
            h => 8000.0, S_T => 0.0555, S_e => 0.010, ρ_p => 32.0,
            σ => 3500.0, w_o => 0.138, δ => 1.0, M_x => 0.12,
            M_f => 0.05, U => 352.0, tan_ϕ => 0.0,
        ]

        prob = ODEProblem(rothermel_simplified, [], (0.0, 1.0), params)
        sol = solve(prob)

        # All rates should be positive
        @test sol[R][end] > 0
        @test sol[I_R][end] > 0
        @test sol[Γ][end] > 0

        # Ratios should be between 0 and 1
        @test 0 < sol[ξ][end] < 1
        @test 0 < sol[β][end] < 1

        # Damping coefficients should be between 0 and 1
        @test 0 <= sol[η_M][end] <= 1
        @test 0 <= sol[η_s][end] <= 1

        # Factors should be non-negative
        @test sol[ϕ_w][end] >= 0
        @test sol[ϕ_s][end] >= 0

        # Check for NaN or Inf
        @test !isnan(sol[R][end])
        @test !isinf(sol[R][end])
    end
end
