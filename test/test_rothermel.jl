@testsnippet RothermelSetup begin
    using Test
    using ModelingToolkit
    using NonlinearSolve
    using Symbolics
    using WildlandFire

    # Unit conversion factors (US customary to SI)
    const ft_to_m = 0.3048
    const lb_to_kg = 0.453592
    const Btu_to_J = 1055.06
    const min_to_s = 60.0
    const lbft2_to_kgm2 = lb_to_kg / ft_to_m^2  # 4.88243
    const lbft3_to_kgm3 = lb_to_kg / ft_to_m^3  # 16.0185
    const invft_to_invm = 1 / ft_to_m           # 3.28084
    const ftmin_to_ms = ft_to_m / min_to_s      # 0.00508
    const Btuft2min_to_Wm2 = Btu_to_J / (ft_to_m^2 * min_to_s)  # 189.27

    # Fuel Model 1 parameters in SI units (converted from US customary)
    # Original US: σ = 3500 1/ft, w0 = 0.034 lb/ft², δ = 1.0 ft, Mx = 0.12
    const σ_SI = 3500.0 * invft_to_invm         # 11483.5 1/m
    const w0_SI = 0.034 * lbft2_to_kgm2         # 0.166 kg/m²
    const δ_SI = 1.0 * ft_to_m                  # 0.3048 m
    const ρ_p_SI = 32.0 * lbft3_to_kgm3         # 512.6 kg/m³
end

@testitem "Structural Verification" setup = [RothermelSetup] tags = [:rothermel] begin
    # Test that RothermelFireSpread creates a valid system
    sys = RothermelFireSpread()

    # Verify system structure
    @test sys !== nothing
    @test length(equations(sys)) == 26
    @test length(unknowns(sys)) == 26

    # Check key output variables exist
    var_names = [string(Symbolics.tosymbol(v, escape = false)) for v in unknowns(sys)]
    @test "R" in var_names      # Rate of spread
    @test "R0" in var_names     # No-wind no-slope rate of spread
    @test "IR" in var_names     # Reaction intensity
    @test "F_L" in var_names    # Flame length
    @test "IB" in var_names     # Fireline intensity

    # Check key parameters exist
    param_names = [string(Symbolics.tosymbol(p, escape = false)) for p in parameters(sys)]
    @test "σ" in param_names    # Surface-area-to-volume ratio
    @test "w0" in param_names   # Fuel load
    @test "δ" in param_names    # Fuel bed depth
    @test "Mx" in param_names   # Moisture of extinction
    @test "Mf" in param_names   # Fuel moisture
    @test "U" in param_names    # Wind speed
    @test "tanϕ" in param_names # Slope
end

@testitem "DynamicFuelLoadTransfer Structure" setup = [RothermelSetup] tags = [:rothermel] begin
    sys = DynamicFuelLoadTransfer()

    @test sys !== nothing
    @test length(equations(sys)) == 3
    @test length(unknowns(sys)) == 3

    var_names = [string(Symbolics.tosymbol(v, escape = false)) for v in unknowns(sys)]
    @test "T_fraction" in var_names
    @test "w0_dead_herb" in var_names
    @test "w0_live_herb_remaining" in var_names
end

@testitem "LiveFuelMoistureExtinction Structure" setup = [RothermelSetup] tags = [:rothermel] begin
    sys = LiveFuelMoistureExtinction()

    @test sys !== nothing
    @test length(equations(sys)) == 1
    @test length(unknowns(sys)) == 1
end

@testitem "Fuel Model 1 - Short Grass" setup = [RothermelSetup] tags = [:rothermel] begin
    # Test with Fuel Model 1 (Short grass) parameters from Table 8
    # Using SI units

    sys = RothermelFireSpread()
    compiled_sys = mtkcompile(sys)

    # Fuel Model 1 parameters in SI units
    prob = NonlinearProblem(
        compiled_sys, Dict(
            compiled_sys.σ => σ_SI,            # SAV ratio (1/m)
            compiled_sys.w0 => w0_SI,          # Fuel load (kg/m²)
            compiled_sys.δ => δ_SI,            # Fuel bed depth (m)
            compiled_sys.Mx => 0.12,           # Moisture of extinction
            compiled_sys.Mf => 0.05,           # 5% moisture content (dry conditions)
            compiled_sys.U => 0.0,             # No wind (m/s)
            compiled_sys.tanϕ => 0.0           # Flat terrain
        )
    )

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

    # Verify intermediate calculations in SI units
    # Bulk density: ρb = w0 / δ = 0.166 / 0.3048 = 0.544 kg/m³
    @test ρb ≈ w0_SI / δ_SI rtol = 1.0e-6

    # Packing ratio: β = ρb / ρp = 0.544 / 512.6 ≈ 0.001063
    @test β ≈ (w0_SI / δ_SI) / ρ_p_SI rtol = 1.0e-4

    # Optimum packing ratio (SI): βop = 8.858 * σ^(-0.8189)
    β_op_expected = 8.858 * σ_SI^(-0.8189)
    @test β_op ≈ β_op_expected rtol = 1.0e-6

    # Effective heating number (SI): ε = exp(-452.7/σ)
    ε_expected = exp(-452.7 / σ_SI)
    @test ε ≈ ε_expected rtol = 1.0e-6

    # Moisture damping coefficient (same as US, dimensionless)
    rM = 0.05 / 0.12
    η_M_expected = 1.0 - 2.59 * rM + 5.11 * rM^2 - 3.52 * rM^3
    @test η_M ≈ η_M_expected rtol = 1.0e-6

    # Mineral damping coefficient (same as US, dimensionless)
    η_s_expected = min(0.174 * 0.01^(-0.19), 1.0)
    @test η_s ≈ η_s_expected rtol = 1.0e-6

    # With no wind and slope, R should equal R0
    @test R ≈ R0 rtol = 1.0e-6

    # R0 should be positive for valid fire spread
    @test R0 > 0

    # IR should be positive
    @test IR > 0
end

@testitem "Wind Effect" setup = [RothermelSetup] tags = [:rothermel] begin
    # Test that wind increases the rate of spread

    sys = RothermelFireSpread()
    compiled_sys = mtkcompile(sys)

    # Base case: no wind
    prob_no_wind = NonlinearProblem(
        compiled_sys, Dict(
            compiled_sys.σ => σ_SI,
            compiled_sys.w0 => w0_SI,
            compiled_sys.δ => δ_SI,
            compiled_sys.Mx => 0.12,
            compiled_sys.Mf => 0.05,
            compiled_sys.U => 0.0,
            compiled_sys.tanϕ => 0.0
        )
    )
    sol_no_wind = solve(prob_no_wind)
    R_no_wind = sol_no_wind[compiled_sys.R]

    # With wind: 5 mi/h = 440 ft/min = 2.235 m/s
    U_wind = 440.0 * ftmin_to_ms
    prob_wind = NonlinearProblem(
        compiled_sys, Dict(
            compiled_sys.σ => σ_SI,
            compiled_sys.w0 => w0_SI,
            compiled_sys.δ => δ_SI,
            compiled_sys.Mx => 0.12,
            compiled_sys.Mf => 0.05,
            compiled_sys.U => U_wind,
            compiled_sys.tanϕ => 0.0
        )
    )
    sol_wind = solve(prob_wind)
    R_wind = sol_wind[compiled_sys.R]

    # Wind should increase rate of spread
    @test R_wind > R_no_wind

    # Wind factor should be positive when wind > 0
    φw = sol_wind[compiled_sys.φw]
    @test φw > 0
end

@testitem "Slope Effect" setup = [RothermelSetup] tags = [:rothermel] begin
    # Test that uphill slope increases the rate of spread

    sys = RothermelFireSpread()
    compiled_sys = mtkcompile(sys)

    # Base case: flat terrain
    prob_flat = NonlinearProblem(
        compiled_sys, Dict(
            compiled_sys.σ => σ_SI,
            compiled_sys.w0 => w0_SI,
            compiled_sys.δ => δ_SI,
            compiled_sys.Mx => 0.12,
            compiled_sys.Mf => 0.05,
            compiled_sys.U => 0.0,
            compiled_sys.tanϕ => 0.0
        )
    )
    sol_flat = solve(prob_flat)
    R_flat = sol_flat[compiled_sys.R]

    # With 30% slope (tan(16.7°) ≈ 0.3)
    prob_slope = NonlinearProblem(
        compiled_sys, Dict(
            compiled_sys.σ => σ_SI,
            compiled_sys.w0 => w0_SI,
            compiled_sys.δ => δ_SI,
            compiled_sys.Mx => 0.12,
            compiled_sys.Mf => 0.05,
            compiled_sys.U => 0.0,
            compiled_sys.tanϕ => 0.3
        )
    )
    sol_slope = solve(prob_slope)
    R_slope = sol_slope[compiled_sys.R]

    # Slope should increase rate of spread
    @test R_slope > R_flat

    # Slope factor should be positive when slope > 0
    φs = sol_slope[compiled_sys.φs]
    @test φs > 0
end

@testitem "Moisture Effect" setup = [RothermelSetup] tags = [:rothermel] begin
    # Test that increased moisture decreases rate of spread

    sys = RothermelFireSpread()
    compiled_sys = mtkcompile(sys)

    # Low moisture: 5%
    prob_dry = NonlinearProblem(
        compiled_sys, Dict(
            compiled_sys.σ => σ_SI,
            compiled_sys.w0 => w0_SI,
            compiled_sys.δ => δ_SI,
            compiled_sys.Mx => 0.12,
            compiled_sys.Mf => 0.05,
            compiled_sys.U => 0.0,
            compiled_sys.tanϕ => 0.0
        )
    )
    sol_dry = solve(prob_dry)
    R_dry = sol_dry[compiled_sys.R]

    # Higher moisture: 10%
    prob_wet = NonlinearProblem(
        compiled_sys, Dict(
            compiled_sys.σ => σ_SI,
            compiled_sys.w0 => w0_SI,
            compiled_sys.δ => δ_SI,
            compiled_sys.Mx => 0.12,
            compiled_sys.Mf => 0.1,
            compiled_sys.U => 0.0,
            compiled_sys.tanϕ => 0.0
        )
    )
    sol_wet = solve(prob_wet)
    R_wet = sol_wet[compiled_sys.R]

    # Higher moisture should decrease rate of spread
    @test R_wet < R_dry

    # Moisture damping coefficient should be lower at higher moisture
    η_M_dry = sol_dry[compiled_sys.η_M]
    η_M_wet = sol_wet[compiled_sys.η_M]
    @test η_M_wet < η_M_dry
end

@testitem "Extinction Moisture" setup = [RothermelSetup] tags = [:rothermel] begin
    # Test behavior near extinction moisture

    sys = RothermelFireSpread()
    compiled_sys = mtkcompile(sys)

    # At extinction moisture, moisture damping should approach zero
    prob = NonlinearProblem(
        compiled_sys, Dict(
            compiled_sys.σ => σ_SI,
            compiled_sys.w0 => w0_SI,
            compiled_sys.δ => δ_SI,
            compiled_sys.Mx => 0.12,
            compiled_sys.Mf => 0.12,  # At extinction moisture
            compiled_sys.U => 0.0,
            compiled_sys.tanϕ => 0.0
        )
    )
    sol = solve(prob)

    # At extinction moisture, η_M should be approximately 0
    η_M = sol[compiled_sys.η_M]
    @test η_M ≈ 0.0 atol = 1.0e-10

    # Consequently, rate of spread should be approximately 0
    R = sol[compiled_sys.R]
    @test R ≈ 0.0 atol = 1.0e-10
end

@testitem "DynamicFuelLoadTransfer Calculations" setup = [RothermelSetup] tags = [:rothermel] begin
    # Test the dynamic fuel load transfer equations

    sys = DynamicFuelLoadTransfer()
    compiled_sys = mtkcompile(sys)

    # Test fuel load in SI units: 0.05 lb/ft² = 0.244 kg/m²
    w0_test_SI = 0.05 * lbft2_to_kgm2

    # Test at low live herbaceous moisture (fully cured)
    # When Mf_live_herb = 0.3 (30%), T = -1.11*0.3 + 1.33 = 0.997
    # (not quite capped at 1.0, need even lower moisture)
    prob_cured = NonlinearProblem(
        compiled_sys, Dict(
            compiled_sys.w0_live_herb => w0_test_SI,
            compiled_sys.Mf_live_herb => 0.3
        )
    )
    sol_cured = solve(prob_cured)

    T_cured = sol_cured[compiled_sys.T_fraction]
    T_expected_cured = -1.11 * 0.3 + 1.33  # = 0.997
    @test T_cured ≈ T_expected_cured rtol = 1.0e-6

    # Most load transferred to dead (0.997 of it)
    w0_dead = sol_cured[compiled_sys.w0_dead_herb]
    @test w0_dead ≈ T_expected_cured * w0_test_SI rtol = 1.0e-6

    # Test at high live herbaceous moisture (green)
    # When Mf_live_herb = 1.2 (120%), T = -1.11*1.2 + 1.33 = 0.0 (capped at 0)
    prob_green = NonlinearProblem(
        compiled_sys, Dict(
            compiled_sys.w0_live_herb => w0_test_SI,
            compiled_sys.Mf_live_herb => 1.2
        )
    )
    sol_green = solve(prob_green)

    T_green = sol_green[compiled_sys.T_fraction]
    @test T_green ≈ 0.0 atol = 1.0e-6

    # No load transferred to dead
    w0_dead_green = sol_green[compiled_sys.w0_dead_herb]
    @test w0_dead_green ≈ 0.0 atol = 1.0e-6

    # Test intermediate case
    # When Mf_live_herb = 0.7 (70%), T = -1.11*0.7 + 1.33 = 0.553
    w0_test2_SI = 0.1 * lbft2_to_kgm2
    prob_mid = NonlinearProblem(
        compiled_sys, Dict(
            compiled_sys.w0_live_herb => w0_test2_SI,
            compiled_sys.Mf_live_herb => 0.7
        )
    )
    sol_mid = solve(prob_mid)

    T_mid = sol_mid[compiled_sys.T_fraction]
    T_expected = -1.11 * 0.7 + 1.33
    @test T_mid ≈ T_expected rtol = 1.0e-6

    w0_dead_mid = sol_mid[compiled_sys.w0_dead_herb]
    @test w0_dead_mid ≈ T_expected * w0_test2_SI rtol = 1.0e-6
end

@testitem "LiveFuelMoistureExtinction Calculations" setup = [RothermelSetup] tags = [:rothermel] begin
    sys = LiveFuelMoistureExtinction()
    compiled_sys = mtkcompile(sys)

    # Test live fuel moisture of extinction calculation
    # Mx_live = 2.9*W*(1 - Mf_dead/Mx_dead) - 0.226 (min = Mx_dead)

    prob = NonlinearProblem(
        compiled_sys, Dict(
            compiled_sys.Mx_dead => 0.25,
            compiled_sys.W_ratio => 1.5,
            compiled_sys.Mf_dead => 0.1
        )
    )
    sol = solve(prob)

    Mx_live = sol[compiled_sys.Mx_live]
    Mx_live_expected = 2.9 * 1.5 * (1 - 0.1 / 0.25) - 0.226
    @test Mx_live ≈ Mx_live_expected rtol = 1.0e-6

    # Test minimum constraint
    prob_low = NonlinearProblem(
        compiled_sys, Dict(
            compiled_sys.Mx_dead => 0.25,
            compiled_sys.W_ratio => 0.1,  # Low ratio should result in low Mx_live
            compiled_sys.Mf_dead => 0.2
        )
    )
    sol_low = solve(prob_low)

    Mx_live_low = sol_low[compiled_sys.Mx_live]
    # Should be bounded by Mx_dead
    @test Mx_live_low >= 0.25
end

@testitem "Flame Length Relationship" setup = [RothermelSetup] tags = [:rothermel] begin
    # Verify the Byram flame length equation in SI units

    sys = RothermelFireSpread()
    compiled_sys = mtkcompile(sys)

    U_wind = 440.0 * ftmin_to_ms  # 5 mi/h in m/s

    prob = NonlinearProblem(
        compiled_sys, Dict(
            compiled_sys.σ => σ_SI,
            compiled_sys.w0 => w0_SI,
            compiled_sys.δ => δ_SI,
            compiled_sys.Mx => 0.12,
            compiled_sys.Mf => 0.05,
            compiled_sys.U => U_wind,
            compiled_sys.tanϕ => 0.0
        )
    )
    sol = solve(prob)

    IB = sol[compiled_sys.IB]
    F_L = sol[compiled_sys.F_L]

    # Verify flame length equation (SI): F_L = 0.00323 * IB^0.46
    F_L_expected = 0.00323 * IB^0.46
    @test F_L ≈ F_L_expected rtol = 1.0e-6

    # Flame length should be positive
    @test F_L > 0
end

@testitem "Residence Time Relationship" setup = [RothermelSetup] tags = [:rothermel] begin
    # Verify residence time equation (SI): t_r = 75590.6 / σ
    # This matches the US equation: t_r = 384/σ (min) where σ in 1/ft
    # Converted: 384 * 60 / 0.3048 = 75590.6 (s, σ in 1/m)

    sys = RothermelFireSpread()
    compiled_sys = mtkcompile(sys)

    prob = NonlinearProblem(
        compiled_sys, Dict(
            compiled_sys.σ => σ_SI,
            compiled_sys.w0 => w0_SI,
            compiled_sys.δ => δ_SI,
            compiled_sys.Mx => 0.12,
            compiled_sys.Mf => 0.05,
            compiled_sys.U => 0.0,
            compiled_sys.tanϕ => 0.0
        )
    )
    sol = solve(prob)

    t_r = sol[compiled_sys.t_r]
    t_r_expected = 75590.6 / σ_SI

    @test t_r ≈ t_r_expected rtol = 1.0e-6

    # Cross-validate with US equation
    σ_US = 3500.0  # 1/ft
    t_r_US_min = 384 / σ_US  # minutes
    t_r_SI_s = t_r  # seconds
    @test t_r_SI_s / 60 ≈ t_r_US_min rtol = 1.0e-4  # Compare in minutes
end

@testitem "Heat Per Unit Area" setup = [RothermelSetup] tags = [:rothermel] begin
    # Verify heat per unit area equation: HA = IR * t_r

    sys = RothermelFireSpread()
    compiled_sys = mtkcompile(sys)

    prob = NonlinearProblem(
        compiled_sys, Dict(
            compiled_sys.σ => σ_SI,
            compiled_sys.w0 => w0_SI,
            compiled_sys.δ => δ_SI,
            compiled_sys.Mx => 0.12,
            compiled_sys.Mf => 0.05,
            compiled_sys.U => 0.0,
            compiled_sys.tanϕ => 0.0
        )
    )
    sol = solve(prob)

    IR = sol[compiled_sys.IR]
    t_r = sol[compiled_sys.t_r]
    HA = sol[compiled_sys.HA]

    HA_expected = IR * t_r
    @test HA ≈ HA_expected rtol = 1.0e-6
end

@testitem "WindLimit Equations" setup = [RothermelSetup] tags = [:rothermel] begin
    # Test both the corrected and original wind limit equations in SI units

    # Test corrected equation (default)
    sys_corrected = WindLimit(use_corrected = true)
    compiled_corrected = mtkcompile(sys_corrected)

    # IR in SI units: 1000 Btu/ft²/min = 189270 W/m²
    IR_SI = 1000.0 * Btuft2min_to_Wm2

    prob_corrected = NonlinearProblem(
        compiled_corrected, Dict(
            compiled_corrected.IR => IR_SI
        )
    )
    sol_corrected = solve(prob_corrected)

    U_limit_corrected = sol_corrected[compiled_corrected.U_limit]
    # SI: U_limit = 0.0857 * IR^(1/3)
    U_limit_expected_corrected = 0.0857 * IR_SI^(1 / 3)
    @test U_limit_corrected ≈ U_limit_expected_corrected rtol = 1.0e-6

    # Test original equation
    sys_original = WindLimit(use_corrected = false)
    compiled_original = mtkcompile(sys_original)

    prob_original = NonlinearProblem(
        compiled_original, Dict(
            compiled_original.IR => IR_SI
        )
    )
    sol_original = solve(prob_original)

    U_limit_original = sol_original[compiled_original.U_limit]
    # SI: U_limit = 2.42e-5 * IR
    U_limit_expected_original = 2.42e-5 * IR_SI
    @test U_limit_original ≈ U_limit_expected_original rtol = 1.0e-6
end

@testitem "EffectiveMidflameWindSpeed Calculations" setup = [RothermelSetup] tags = [:rothermel] begin
    # Test the effective midflame wind speed calculation
    # This component converts combined wind and slope factors back to an equivalent wind speed

    sys = EffectiveMidflameWindSpeed()
    compiled_sys = mtkcompile(sys)

    # Test with typical values from a fire spread calculation
    # C = 7.47, B = 0.15, E = 0.5, β_ratio = 0.5, φw = 10.0, φs = 2.0
    prob = NonlinearProblem(
        compiled_sys, Dict(
            compiled_sys.C_coeff => 7.47,
            compiled_sys.B_coeff => 0.15,
            compiled_sys.E_coeff => 0.5,
            compiled_sys.β_ratio => 0.5,
            compiled_sys.φw => 10.0,
            compiled_sys.φs => 2.0
        )
    )
    sol = solve(prob)

    # Combined factor should be sum of wind and slope factors
    φE = sol[compiled_sys.φE]
    @test φE ≈ 12.0 rtol = 1.0e-6

    # Effective wind speed should be positive
    UE = sol[compiled_sys.UE]
    @test UE > 0

    # Test the inverse relationship: UE = [(φE * β_ratio^E) / C]^(1/B)
    UE_expected = (φE * 0.5^0.5 / 7.47)^(1 / 0.15)
    @test UE ≈ UE_expected rtol = 1.0e-6

    # Test with no slope factor
    prob_no_slope = NonlinearProblem(
        compiled_sys, Dict(
            compiled_sys.C_coeff => 7.47,
            compiled_sys.B_coeff => 0.15,
            compiled_sys.E_coeff => 0.5,
            compiled_sys.β_ratio => 0.5,
            compiled_sys.φw => 10.0,
            compiled_sys.φs => 0.0
        )
    )
    sol_no_slope = solve(prob_no_slope)

    UE_no_slope = sol_no_slope[compiled_sys.UE]
    # With less combined factor, effective wind should be lower
    @test UE_no_slope < UE
end

@testitem "Cross-validation with US Customary Units" setup = [RothermelSetup] tags = [:rothermel] begin
    # Cross-validate that SI implementation gives equivalent results to US customary
    # by comparing spread rates after unit conversion

    sys = RothermelFireSpread()
    compiled_sys = mtkcompile(sys)

    # Fuel Model 1 in SI
    prob = NonlinearProblem(
        compiled_sys, Dict(
            compiled_sys.σ => σ_SI,
            compiled_sys.w0 => w0_SI,
            compiled_sys.δ => δ_SI,
            compiled_sys.Mx => 0.12,
            compiled_sys.Mf => 0.05,
            compiled_sys.U => 0.0,
            compiled_sys.tanϕ => 0.0
        )
    )
    sol = solve(prob)

    R_SI = sol[compiled_sys.R]  # m/s

    # Convert back to US customary (ft/min) for comparison
    R_US = R_SI / ftmin_to_ms

    # For Fuel Model 1 (short grass), no wind, 5% moisture, the spread rate
    # should be approximately 1-5 ft/min based on published Rothermel results
    # This is a sanity check that the unit conversion preserves reasonable values
    @test R_US > 0.5   # Minimum reasonable spread rate
    @test R_US < 20.0  # Maximum reasonable spread rate for short grass without wind
end
