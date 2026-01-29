@testsnippet RothermelSetup begin
    using Test
    using ModelingToolkit
    using NonlinearSolve
    using Symbolics
    using WildlandFire
end

@testitem "Structural Verification" setup=[RothermelSetup] tags=[:rothermel] begin
    # Test that RothermelFireSpread creates a valid system
    sys = RothermelFireSpread()

    # Verify system structure
    @test sys !== nothing
    @test length(equations(sys)) == 26
    @test length(unknowns(sys)) == 26

    # Check key output variables exist
    var_names = [string(Symbolics.tosymbol(v, escape=false)) for v in unknowns(sys)]
    @test "R" in var_names      # Rate of spread
    @test "R0" in var_names     # No-wind no-slope rate of spread
    @test "IR" in var_names     # Reaction intensity
    @test "F_L" in var_names    # Flame length
    @test "IB" in var_names     # Fireline intensity

    # Check key parameters exist
    param_names = [string(Symbolics.tosymbol(p, escape=false)) for p in parameters(sys)]
    @test "σ" in param_names    # Surface-area-to-volume ratio
    @test "w0" in param_names   # Fuel load
    @test "δ" in param_names    # Fuel bed depth
    @test "Mx" in param_names   # Moisture of extinction
    @test "Mf" in param_names   # Fuel moisture
    @test "U" in param_names    # Wind speed
    @test "tanϕ" in param_names # Slope
end

@testitem "DynamicFuelLoadTransfer Structure" setup=[RothermelSetup] tags=[:rothermel] begin
    sys = DynamicFuelLoadTransfer()

    @test sys !== nothing
    @test length(equations(sys)) == 3
    @test length(unknowns(sys)) == 3

    var_names = [string(Symbolics.tosymbol(v, escape=false)) for v in unknowns(sys)]
    @test "T_fraction" in var_names
    @test "w0_dead_herb" in var_names
    @test "w0_live_herb_remaining" in var_names
end

@testitem "LiveFuelMoistureExtinction Structure" setup=[RothermelSetup] tags=[:rothermel] begin
    sys = LiveFuelMoistureExtinction()

    @test sys !== nothing
    @test length(equations(sys)) == 1
    @test length(unknowns(sys)) == 1
end

@testitem "Fuel Model 1 - Short Grass" setup=[RothermelSetup] tags=[:rothermel] begin
    # Test with Fuel Model 1 (Short grass) parameters from Table 8
    # Expected values computed from the Rothermel equations

    sys = RothermelFireSpread()
    compiled_sys = mtkcompile(sys)

    # Fuel Model 1 parameters (Short grass, 1 ft depth)
    # From Table 8: σ = 3500 1/ft, w0 = 0.034 lb/ft², δ = 1.0 ft, Mx = 0.12
    prob = NonlinearProblem(compiled_sys, [], [
        compiled_sys.σ => 3500.0,      # SAV ratio (1/ft)
        compiled_sys.w0 => 0.034,      # Fuel load (lb/ft²) = 0.74 ton/acre
        compiled_sys.δ => 1.0,         # Fuel bed depth (ft)
        compiled_sys.Mx => 0.12,       # Moisture of extinction
        compiled_sys.Mf => 0.05,       # 5% moisture content (dry conditions)
        compiled_sys.U => 0.0,         # No wind
        compiled_sys.tanϕ => 0.0       # Flat terrain
    ])

    sol = solve(prob)

    # Get solution values
    R0 = sol[compiled_sys.R0]
    R = sol[compiled_sys.R]
    IR = sol[compiled_sys.IR]
    ρb = sol[compiled_sys.ρb]
    β = sol[compiled_sys.β]
    β_op = sol[compiled_sys.β_op]
    ε = sol[compiled_sys.ε]
    ξ = sol[compiled_sys.ξ]
    η_M = sol[compiled_sys.η_M]
    η_s = sol[compiled_sys.η_s]

    # Verify intermediate calculations
    # Bulk density: ρb = w0 / δ = 0.034 / 1.0 = 0.034 lb/ft³
    @test ρb ≈ 0.034 rtol=1e-6

    # Packing ratio: β = ρb / ρp = 0.034 / 32 ≈ 0.001063
    @test β ≈ 0.034 / 32.0 rtol=1e-6

    # Optimum packing ratio: βop = 3.348 * σ^(-0.8189)
    β_op_expected = 3.348 * 3500.0^(-0.8189)
    @test β_op ≈ β_op_expected rtol=1e-6

    # Effective heating number: ε = exp(-138/σ)
    ε_expected = exp(-138.0 / 3500.0)
    @test ε ≈ ε_expected rtol=1e-6

    # Moisture damping coefficient
    rM = 0.05 / 0.12
    η_M_expected = 1.0 - 2.59*rM + 5.11*rM^2 - 3.52*rM^3
    @test η_M ≈ η_M_expected rtol=1e-6

    # Mineral damping coefficient
    η_s_expected = min(0.174 * 0.010^(-0.19), 1.0)
    @test η_s ≈ η_s_expected rtol=1e-6

    # With no wind and slope, R should equal R0
    @test R ≈ R0 rtol=1e-6

    # R0 should be positive for valid fire spread
    @test R0 > 0

    # IR should be positive
    @test IR > 0
end

@testitem "Wind Effect" setup=[RothermelSetup] tags=[:rothermel] begin
    # Test that wind increases the rate of spread

    sys = RothermelFireSpread()
    compiled_sys = mtkcompile(sys)

    # Base case: no wind
    prob_no_wind = NonlinearProblem(compiled_sys, [], [
        compiled_sys.σ => 3500.0,
        compiled_sys.w0 => 0.034,
        compiled_sys.δ => 1.0,
        compiled_sys.Mx => 0.12,
        compiled_sys.Mf => 0.05,
        compiled_sys.U => 0.0,
        compiled_sys.tanϕ => 0.0
    ])
    sol_no_wind = solve(prob_no_wind)
    R_no_wind = sol_no_wind[compiled_sys.R]

    # With wind: 5 mi/h = 440 ft/min
    prob_wind = NonlinearProblem(compiled_sys, [], [
        compiled_sys.σ => 3500.0,
        compiled_sys.w0 => 0.034,
        compiled_sys.δ => 1.0,
        compiled_sys.Mx => 0.12,
        compiled_sys.Mf => 0.05,
        compiled_sys.U => 440.0,
        compiled_sys.tanϕ => 0.0
    ])
    sol_wind = solve(prob_wind)
    R_wind = sol_wind[compiled_sys.R]

    # Wind should increase rate of spread
    @test R_wind > R_no_wind

    # Wind factor should be positive when wind > 0
    φw = sol_wind[compiled_sys.φw]
    @test φw > 0
end

@testitem "Slope Effect" setup=[RothermelSetup] tags=[:rothermel] begin
    # Test that uphill slope increases the rate of spread

    sys = RothermelFireSpread()
    compiled_sys = mtkcompile(sys)

    # Base case: flat terrain
    prob_flat = NonlinearProblem(compiled_sys, [], [
        compiled_sys.σ => 3500.0,
        compiled_sys.w0 => 0.034,
        compiled_sys.δ => 1.0,
        compiled_sys.Mx => 0.12,
        compiled_sys.Mf => 0.05,
        compiled_sys.U => 0.0,
        compiled_sys.tanϕ => 0.0
    ])
    sol_flat = solve(prob_flat)
    R_flat = sol_flat[compiled_sys.R]

    # With 30% slope (tan(16.7°) ≈ 0.3)
    prob_slope = NonlinearProblem(compiled_sys, [], [
        compiled_sys.σ => 3500.0,
        compiled_sys.w0 => 0.034,
        compiled_sys.δ => 1.0,
        compiled_sys.Mx => 0.12,
        compiled_sys.Mf => 0.05,
        compiled_sys.U => 0.0,
        compiled_sys.tanϕ => 0.3
    ])
    sol_slope = solve(prob_slope)
    R_slope = sol_slope[compiled_sys.R]

    # Slope should increase rate of spread
    @test R_slope > R_flat

    # Slope factor should be positive when slope > 0
    φs = sol_slope[compiled_sys.φs]
    @test φs > 0
end

@testitem "Moisture Effect" setup=[RothermelSetup] tags=[:rothermel] begin
    # Test that increased moisture decreases rate of spread

    sys = RothermelFireSpread()
    compiled_sys = mtkcompile(sys)

    # Low moisture: 5%
    prob_dry = NonlinearProblem(compiled_sys, [], [
        compiled_sys.σ => 3500.0,
        compiled_sys.w0 => 0.034,
        compiled_sys.δ => 1.0,
        compiled_sys.Mx => 0.12,
        compiled_sys.Mf => 0.05,
        compiled_sys.U => 0.0,
        compiled_sys.tanϕ => 0.0
    ])
    sol_dry = solve(prob_dry)
    R_dry = sol_dry[compiled_sys.R]

    # Higher moisture: 10%
    prob_wet = NonlinearProblem(compiled_sys, [], [
        compiled_sys.σ => 3500.0,
        compiled_sys.w0 => 0.034,
        compiled_sys.δ => 1.0,
        compiled_sys.Mx => 0.12,
        compiled_sys.Mf => 0.10,
        compiled_sys.U => 0.0,
        compiled_sys.tanϕ => 0.0
    ])
    sol_wet = solve(prob_wet)
    R_wet = sol_wet[compiled_sys.R]

    # Higher moisture should decrease rate of spread
    @test R_wet < R_dry

    # Moisture damping coefficient should be lower at higher moisture
    η_M_dry = sol_dry[compiled_sys.η_M]
    η_M_wet = sol_wet[compiled_sys.η_M]
    @test η_M_wet < η_M_dry
end

@testitem "Extinction Moisture" setup=[RothermelSetup] tags=[:rothermel] begin
    # Test behavior near extinction moisture

    sys = RothermelFireSpread()
    compiled_sys = mtkcompile(sys)

    # At extinction moisture, moisture damping should approach zero
    prob = NonlinearProblem(compiled_sys, [], [
        compiled_sys.σ => 3500.0,
        compiled_sys.w0 => 0.034,
        compiled_sys.δ => 1.0,
        compiled_sys.Mx => 0.12,
        compiled_sys.Mf => 0.12,  # At extinction moisture
        compiled_sys.U => 0.0,
        compiled_sys.tanϕ => 0.0
    ])
    sol = solve(prob)

    # At extinction moisture, η_M should be approximately 0
    η_M = sol[compiled_sys.η_M]
    @test η_M ≈ 0.0 atol=1e-10

    # Consequently, rate of spread should be approximately 0
    R = sol[compiled_sys.R]
    @test R ≈ 0.0 atol=1e-10
end

@testitem "DynamicFuelLoadTransfer Calculations" setup=[RothermelSetup] tags=[:rothermel] begin
    # Test the dynamic fuel load transfer equations

    sys = DynamicFuelLoadTransfer()
    compiled_sys = mtkcompile(sys)

    # Test at low live herbaceous moisture (fully cured)
    # When Mf_live_herb = 0.3 (30%), T = -1.11*0.3 + 1.33 = 1.0 (capped)
    prob_cured = NonlinearProblem(compiled_sys, [], [
        compiled_sys.w0_live_herb => 0.05,
        compiled_sys.Mf_live_herb => 0.3
    ])
    sol_cured = solve(prob_cured)

    T_cured = sol_cured[compiled_sys.T_fraction]
    @test T_cured ≈ 1.0 rtol=1e-6

    # All load transferred to dead
    w0_dead = sol_cured[compiled_sys.w0_dead_herb]
    @test w0_dead ≈ 0.05 rtol=1e-6

    # Test at high live herbaceous moisture (green)
    # When Mf_live_herb = 1.2 (120%), T = -1.11*1.2 + 1.33 = 0.0 (capped at 0)
    prob_green = NonlinearProblem(compiled_sys, [], [
        compiled_sys.w0_live_herb => 0.05,
        compiled_sys.Mf_live_herb => 1.2
    ])
    sol_green = solve(prob_green)

    T_green = sol_green[compiled_sys.T_fraction]
    @test T_green ≈ 0.0 atol=1e-6

    # No load transferred to dead
    w0_dead_green = sol_green[compiled_sys.w0_dead_herb]
    @test w0_dead_green ≈ 0.0 atol=1e-6

    # Test intermediate case
    # When Mf_live_herb = 0.7 (70%), T = -1.11*0.7 + 1.33 = 0.553
    prob_mid = NonlinearProblem(compiled_sys, [], [
        compiled_sys.w0_live_herb => 0.10,
        compiled_sys.Mf_live_herb => 0.7
    ])
    sol_mid = solve(prob_mid)

    T_mid = sol_mid[compiled_sys.T_fraction]
    T_expected = -1.11 * 0.7 + 1.33
    @test T_mid ≈ T_expected rtol=1e-6

    w0_dead_mid = sol_mid[compiled_sys.w0_dead_herb]
    @test w0_dead_mid ≈ T_expected * 0.10 rtol=1e-6
end

@testitem "LiveFuelMoistureExtinction Calculations" setup=[RothermelSetup] tags=[:rothermel] begin
    sys = LiveFuelMoistureExtinction()
    compiled_sys = mtkcompile(sys)

    # Test live fuel moisture of extinction calculation
    # Mx_live = 2.9*W*(1 - Mf_dead/Mx_dead) - 0.226 (min = Mx_dead)

    prob = NonlinearProblem(compiled_sys, [], [
        compiled_sys.Mx_dead => 0.25,
        compiled_sys.W_ratio => 1.5,
        compiled_sys.Mf_dead => 0.10
    ])
    sol = solve(prob)

    Mx_live = sol[compiled_sys.Mx_live]
    Mx_live_expected = 2.9 * 1.5 * (1 - 0.10/0.25) - 0.226
    @test Mx_live ≈ Mx_live_expected rtol=1e-6

    # Test minimum constraint
    prob_low = NonlinearProblem(compiled_sys, [], [
        compiled_sys.Mx_dead => 0.25,
        compiled_sys.W_ratio => 0.1,  # Low ratio should result in low Mx_live
        compiled_sys.Mf_dead => 0.20
    ])
    sol_low = solve(prob_low)

    Mx_live_low = sol_low[compiled_sys.Mx_live]
    # Should be bounded by Mx_dead
    @test Mx_live_low >= 0.25
end

@testitem "Flame Length Relationship" setup=[RothermelSetup] tags=[:rothermel] begin
    # Verify the Byram flame length equation: F_L = 0.45 * IB^0.46

    sys = RothermelFireSpread()
    compiled_sys = mtkcompile(sys)

    prob = NonlinearProblem(compiled_sys, [], [
        compiled_sys.σ => 3500.0,
        compiled_sys.w0 => 0.034,
        compiled_sys.δ => 1.0,
        compiled_sys.Mx => 0.12,
        compiled_sys.Mf => 0.05,
        compiled_sys.U => 440.0,
        compiled_sys.tanϕ => 0.0
    ])
    sol = solve(prob)

    IB = sol[compiled_sys.IB]
    F_L = sol[compiled_sys.F_L]

    # Verify flame length equation
    F_L_expected = 0.45 * IB^0.46
    @test F_L ≈ F_L_expected rtol=1e-6

    # Flame length should be positive
    @test F_L > 0
end

@testitem "Residence Time Relationship" setup=[RothermelSetup] tags=[:rothermel] begin
    # Verify residence time equation: t_r = 384 / σ

    sys = RothermelFireSpread()
    compiled_sys = mtkcompile(sys)

    prob = NonlinearProblem(compiled_sys, [], [
        compiled_sys.σ => 3500.0,
        compiled_sys.w0 => 0.034,
        compiled_sys.δ => 1.0,
        compiled_sys.Mx => 0.12,
        compiled_sys.Mf => 0.05,
        compiled_sys.U => 0.0,
        compiled_sys.tanϕ => 0.0
    ])
    sol = solve(prob)

    t_r = sol[compiled_sys.t_r]
    t_r_expected = 384.0 / 3500.0

    @test t_r ≈ t_r_expected rtol=1e-6
end

@testitem "Heat Per Unit Area" setup=[RothermelSetup] tags=[:rothermel] begin
    # Verify heat per unit area equation: HA = IR * t_r

    sys = RothermelFireSpread()
    compiled_sys = mtkcompile(sys)

    prob = NonlinearProblem(compiled_sys, [], [
        compiled_sys.σ => 3500.0,
        compiled_sys.w0 => 0.034,
        compiled_sys.δ => 1.0,
        compiled_sys.Mx => 0.12,
        compiled_sys.Mf => 0.05,
        compiled_sys.U => 0.0,
        compiled_sys.tanϕ => 0.0
    ])
    sol = solve(prob)

    IR = sol[compiled_sys.IR]
    t_r = sol[compiled_sys.t_r]
    HA = sol[compiled_sys.HA]

    HA_expected = IR * t_r
    @test HA ≈ HA_expected rtol=1e-6
end

@testitem "WindLimit Equations" setup=[RothermelSetup] tags=[:rothermel] begin
    # Test both the corrected and original wind limit equations

    # Test corrected equation (default)
    sys_corrected = WindLimit(use_corrected=true)
    compiled_corrected = mtkcompile(sys_corrected)

    prob_corrected = NonlinearProblem(compiled_corrected, [], [
        compiled_corrected.IR => 1000.0
    ])
    sol_corrected = solve(prob_corrected)

    U_limit_corrected = sol_corrected[compiled_corrected.U_limit]
    U_limit_expected_corrected = 96.8 * 1000.0^(1/3)
    @test U_limit_corrected ≈ U_limit_expected_corrected rtol=1e-6

    # Test original equation
    sys_original = WindLimit(use_corrected=false)
    compiled_original = mtkcompile(sys_original)

    prob_original = NonlinearProblem(compiled_original, [], [
        compiled_original.IR => 1000.0
    ])
    sol_original = solve(prob_original)

    U_limit_original = sol_original[compiled_original.U_limit]
    U_limit_expected_original = 0.9 * 1000.0
    @test U_limit_original ≈ U_limit_expected_original rtol=1e-6
end
