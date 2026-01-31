# NFDRS Tests for WildlandFire.jl
# Tests for the National Fire Danger Rating System implementation
# Based on Cohen & Deeming (1985)

@testsnippet NFDRSSetup begin
    using WildlandFire
    using Test
    using ModelingToolkit
    using ModelingToolkit: mtkcompile
    using OrdinaryDiffEqDefault
    using Symbolics

    # Unit conversion factors for tests (matching nfdrs.jl)
    const FT_TO_M = 0.3048
    const M_TO_FT = 1.0 / FT_TO_M
    const TONS_PER_ACRE_TO_KG_PER_M2 = 0.2241702
    const BTU_PER_LB_TO_J_PER_KG = 2326.0
    const FPM_TO_MS = 0.00508
    const MS_TO_FPM = 196.8504
end

@testitem "NFDRS Fuel Models" setup=[NFDRSSetup] tags=[:nfdrs] begin
    # Test that all fuel models are available
    @test length(NFDRS_FUEL_MODELS) == 20

    # Test get_fuel_model function (uses Symbol keys)
    model_a = get_fuel_model(:A)
    @test model_a.description == "Western grasses (annual)"

    # Test all fuel model codes exist
    for code in [:A, :B, :C, :D, :E, :F, :G, :H, :I, :J,
                 :K, :L, :N, :O, :P, :Q, :R, :S, :T, :U]
        @test haskey(NFDRS_FUEL_MODELS, code)
    end

    # Test that M is not in the fuel models (skipped in NFDRS)
    @test !haskey(NFDRS_FUEL_MODELS, :M)

    # Test that fuel models have expected properties
    @test model_a.SG1 > 0
    @test model_a.DEPTH > 0

    # Verify specific fuel model values from Appendix (Cohen & Deeming 1985, p. 15)
    # Values are now in SI units
    @test model_a.SG1 ≈ 3000.0 * M_TO_FT atol=0.1  # m⁻¹ (3000 ft⁻¹ converted)
    @test model_a.W1 ≈ 0.20 * TONS_PER_ACRE_TO_KG_PER_M2 atol=1e-6  # kg/m²
    @test model_a.DEPTH ≈ 0.80 * FT_TO_M atol=1e-6  # m
    @test model_a.MXD == 15.0     # Dead fuel moisture of extinction (%) - unchanged
    @test model_a.HD ≈ 8000.0 * BTU_PER_LB_TO_J_PER_KG atol=1.0  # J/kg
    @test model_a.SCM ≈ 300.0 * FPM_TO_MS atol=1e-6  # m/s
    @test model_a.WNDFC == 0.6    # Wind reduction factor - dimensionless, unchanged

    # Test California chaparral (Model B)
    model_b = get_fuel_model(:B)
    @test model_b.SG1 ≈ 700.0 * M_TO_FT atol=0.1
    @test model_b.WWOOD ≈ 11.50 * TONS_PER_ACRE_TO_KG_PER_M2 atol=1e-5  # Heavy woody load
    @test model_b.HD ≈ 9500.0 * BTU_PER_LB_TO_J_PER_KG atol=1.0  # Higher heat content for chaparral
end

@testitem "Component Creation" setup=[NFDRSSetup] tags=[:nfdrs] begin
    # Test that components can be created (basic smoke test)
    @test EquilibriumMoistureContent() isa System
    @test OneHourFuelMoisture() isa System
    @test TenHourFuelMoisture() isa System
    @test HundredHourFuelMoisture() isa System
    @test ThousandHourFuelMoisture() isa System
    @test HerbaceousFuelMoisture() isa System
    @test WoodyFuelMoisture() isa System
    @test FuelLoadingTransfer() isa System
    @test SpreadComponent() isa System
    @test EnergyReleaseComponent() isa System
    @test BurningIndex() isa System
    @test IgnitionComponent() isa System
    @test HumanFireOccurrenceIndex() isa System
    @test LightningFireOccurrenceIndex() isa System
    @test FireLoadIndex() isa System
end

@testitem "EquilibriumMoistureContent Implementation" setup=[NFDRSSetup] tags=[:nfdrs] begin
    # Test the actual implementation against expected values from
    # Cohen & Deeming (1985), Eq. 1a, 1b, 1c (page 1)
    emc = EquilibriumMoistureContent()
    emc_compiled = mtkcompile(emc)

    # Test case 1: Low RH (Eq. 1a: RH < 10%)
    # At RH=5%, TEMP=70°F (294.26 K)
    # EMC = 0.03229 + 0.281073*5 - 0.000578*70*5 = 1.236%
    temp_K_low = (70.0 + 459.67) * 5 / 9  # 70°F to K
    prob1 = ODEProblem(emc_compiled,
        Dict(emc_compiled.TEMP => temp_K_low, emc_compiled.RH => 0.05),
        (0.0, 1.0))
    sol1 = solve(prob1)
    expected_emc_low = 0.01236  # 1.236% as fraction
    @test sol1[emc_compiled.EMC][end] ≈ expected_emc_low atol=0.002

    # Test case 2: Mid RH (Eq. 1b: 10% <= RH < 50%)
    # At RH=30%, TEMP=70°F
    # EMC = 2.22749 + 0.160107*30 - 0.014784*70 = 5.997%
    prob2 = ODEProblem(emc_compiled,
        Dict(emc_compiled.TEMP => temp_K_low, emc_compiled.RH => 0.30),
        (0.0, 1.0))
    sol2 = solve(prob2)
    expected_emc_mid = 0.05997  # 5.997% as fraction
    @test sol2[emc_compiled.EMC][end] ≈ expected_emc_mid atol=0.002

    # Test case 3: High RH (Eq. 1c: RH >= 50%)
    # At RH=80%, TEMP=70°F
    # EMC = 21.0606 + 0.005565*6400 - 0.00035*80*70 - 0.483199*80 = 16.06%
    prob3 = ODEProblem(emc_compiled,
        Dict(emc_compiled.TEMP => temp_K_low, emc_compiled.RH => 0.80),
        (0.0, 1.0))
    sol3 = solve(prob3)
    expected_emc_high = 0.1606  # 16.06% as fraction
    @test sol3[emc_compiled.EMC][end] ≈ expected_emc_high atol=0.01
end

@testitem "OneHourFuelMoisture Implementation" setup=[NFDRSSetup] tags=[:nfdrs] begin
    # Test against Cohen & Deeming (1985), page 3-4
    fm1 = OneHourFuelMoisture()
    fm1_compiled = mtkcompile(fm1)

    # Test 1: When raining, MC1 = 35% (page 4)
    prob_rain = ODEProblem(fm1_compiled,
        Dict(fm1_compiled.EMCPRM => 0.10,
             fm1_compiled.MC10 => 0.10,
             fm1_compiled.use_fuel_sticks => 0.0,
             fm1_compiled.is_raining => 1.0),
        (0.0, 1.0))
    sol_rain = solve(prob_rain)
    @test sol_rain[fm1_compiled.MC1][end] ≈ 0.35 atol=0.001

    # Test 2: Without sticks, MC1 = 1.03 * EMCPRM (page 3)
    prob_no_sticks = ODEProblem(fm1_compiled,
        Dict(fm1_compiled.EMCPRM => 0.10,
             fm1_compiled.MC10 => 0.15,
             fm1_compiled.use_fuel_sticks => 0.0,
             fm1_compiled.is_raining => 0.0),
        (0.0, 1.0))
    sol_no_sticks = solve(prob_no_sticks)
    @test sol_no_sticks[fm1_compiled.MC1][end] ≈ 1.03 * 0.10 atol=0.001

    # Test 3: With sticks, MC1 = (4.0 * EMCPRM + MC10) / 5.0 (page 3)
    prob_sticks = ODEProblem(fm1_compiled,
        Dict(fm1_compiled.EMCPRM => 0.10,
             fm1_compiled.MC10 => 0.15,
             fm1_compiled.use_fuel_sticks => 1.0,
             fm1_compiled.is_raining => 0.0),
        (0.0, 1.0))
    sol_sticks = solve(prob_sticks)
    expected_mc1_sticks = (4.0 * 0.10 + 0.15) / 5.0
    @test sol_sticks[fm1_compiled.MC1][end] ≈ expected_mc1_sticks atol=0.001
end

@testitem "TenHourFuelMoisture Implementation" setup=[NFDRSSetup] tags=[:nfdrs] begin
    # Test against Cohen & Deeming (1985), page 4
    fm10 = TenHourFuelMoisture()
    fm10_compiled = mtkcompile(fm10)

    # Test: Without fuel sticks, MC10 = 1.28 * EMCPRM (page 4)
    prob = ODEProblem(fm10_compiled,
        Dict(fm10_compiled.EMCPRM => 0.10,
             fm10_compiled.use_fuel_sticks => 0.0,
             fm10_compiled.WT => 0.1,
             fm10_compiled.AGE => 0.0,
             fm10_compiled.CLIMAT => 2.0),
        (0.0, 1.0))
    sol = solve(prob)
    @test sol[fm10_compiled.MC10][end] ≈ 1.28 * 0.10 atol=0.001
end

@testitem "FuelLoadingTransfer Implementation" setup=[NFDRSSetup] tags=[:nfdrs] begin
    # Test equations 5, 6, 7, 8 from Cohen & Deeming (1985), page 7
    flt = FuelLoadingTransfer()
    flt_compiled = mtkcompile(flt)

    # Test at MCHERB = 30% (0.30 as fraction)
    # FCTCUR = 1.33 - 0.0111 * 30 = 0.997
    W1_test = 0.05  # kg/m²
    WHERB_test = 0.02  # kg/m²
    prob = ODEProblem(flt_compiled,
        Dict(flt_compiled.MCHERB => 0.30,
             flt_compiled.W1 => W1_test,
             flt_compiled.WHERB => WHERB_test),
        (0.0, 1.0))
    sol = solve(prob)

    expected_fctcur = min(1.0, 1.33 - 0.0111 * 30)  # 0.997
    @test sol[flt_compiled.FCTCUR][end] ≈ expected_fctcur atol=0.01

    expected_wherbc = expected_fctcur * WHERB_test
    @test sol[flt_compiled.WHERBC][end] ≈ expected_wherbc atol=0.001

    expected_w1p = W1_test + expected_wherbc
    @test sol[flt_compiled.W1P][end] ≈ expected_w1p atol=0.001

    expected_wherbp = WHERB_test - expected_wherbc
    @test sol[flt_compiled.WHERBP][end] ≈ expected_wherbp atol=0.001

    # Test at MCHERB = 120% (0.0% transfer - fully green)
    prob_green = ODEProblem(flt_compiled,
        Dict(flt_compiled.MCHERB => 1.20,
             flt_compiled.W1 => W1_test,
             flt_compiled.WHERB => WHERB_test),
        (0.0, 1.0))
    sol_green = solve(prob_green)
    @test sol_green[flt_compiled.FCTCUR][end] ≈ 0.0 atol=0.01
end

@testitem "BurningIndex Implementation" setup=[NFDRSSetup] tags=[:nfdrs] begin
    # Test BI = 3.01 * (SC * ERC)^0.46 (page 12)
    # SC is in ft/min for the original equation, so we need to convert
    bi = BurningIndex()
    bi_compiled = mtkcompile(bi)

    # Test case: SC=50 ft/min converted to m/s, ERC=20
    SC_fpm = 50.0
    SC_ms = SC_fpm * FPM_TO_MS
    ERC_val = 20.0

    prob = ODEProblem(bi_compiled,
        Dict(bi_compiled.SC => SC_ms,
             bi_compiled.ERC => ERC_val,
             bi_compiled.fuels_wet => 0.0),
        (0.0, 1.0))
    sol = solve(prob)

    # Expected BI using original ft/min units for SC
    expected_bi = 3.01 * (SC_fpm * ERC_val)^0.46
    @test sol[bi_compiled.BI][end] ≈ expected_bi atol=1.0

    # Test wet conditions - BI should be 0
    prob_wet = ODEProblem(bi_compiled,
        Dict(bi_compiled.SC => SC_ms,
             bi_compiled.ERC => ERC_val,
             bi_compiled.fuels_wet => 1.0),
        (0.0, 1.0))
    sol_wet = solve(prob_wet)
    @test sol_wet[bi_compiled.BI][end] ≈ 0.0 atol=0.001
end

@testitem "FireLoadIndex Implementation" setup=[NFDRSSetup] tags=[:nfdrs] begin
    # Test FLI = 0.71 * sqrt(BI^2 + (LOI + MCOI)^2) (page 14)
    fli = FireLoadIndex()
    fli_compiled = mtkcompile(fli)

    # Test case: BI=50, LOI=30, MCOI=20
    prob = ODEProblem(fli_compiled,
        Dict(fli_compiled.BI => 50.0,
             fli_compiled.LOI => 30.0,
             fli_compiled.MCOI => 20.0),
        (0.0, 1.0))
    sol = solve(prob)

    expected_fli = 0.71 * sqrt(50.0^2 + (30.0 + 20.0)^2)
    @test sol[fli_compiled.FLI][end] ≈ expected_fli atol=0.1

    # Test with BI limited to 100
    prob_limit = ODEProblem(fli_compiled,
        Dict(fli_compiled.BI => 150.0,  # Should be limited to 100
             fli_compiled.LOI => 30.0,
             fli_compiled.MCOI => 30.0),
        (0.0, 1.0))
    sol_limit = solve(prob_limit)

    expected_fli_limit = 0.71 * sqrt(100.0^2 + min(100.0, 60.0)^2)
    @test sol_limit[fli_compiled.FLI][end] ≈ expected_fli_limit atol=0.1
end

@testitem "HumanFireOccurrenceIndex Implementation" setup=[NFDRSSetup] tags=[:nfdrs] begin
    # Test MCOI = 0.01 * MRISK * IC (page 13)
    mcoi = HumanFireOccurrenceIndex()
    mcoi_compiled = mtkcompile(mcoi)

    # Test case: MRISK=50, IC=80
    prob = ODEProblem(mcoi_compiled,
        Dict(mcoi_compiled.MRISK => 50.0,
             mcoi_compiled.IC => 80.0),
        (0.0, 1.0))
    sol = solve(prob)

    expected_mcoi = 0.01 * 50.0 * 80.0
    @test sol[mcoi_compiled.MCOI][end] ≈ expected_mcoi atol=0.1
end

@testitem "HerbaceousFuelMoisture Climate Parameters" setup=[NFDRSSetup] tags=[:nfdrs] begin
    # Verify the implementation produces expected values for different climate classes
    # Based on Cohen & Deeming (1985), page 7-8
    hfm = HerbaceousFuelMoisture()
    hfm_compiled = mtkcompile(hfm)

    # Test green stage at X1000 = 15% (0.15)
    # For climate class 1: MCHERB = (-70 + 12.8 * 15) / 100 = 1.22 (122%)
    prob_c1 = ODEProblem(hfm_compiled,
        Dict(hfm_compiled.X1000 => 0.15,
             hfm_compiled.MC1 => 0.05,
             hfm_compiled.CLIMAT => 1.0,
             hfm_compiled.is_annual => 0.0,
             hfm_compiled.GRNDAY => 0.0,
             hfm_compiled.is_greenup => 0.0,
             hfm_compiled.is_cured => 0.0),
        (0.0, 1.0))
    sol_c1 = solve(prob_c1)

    expected_mcherb_c1 = (-70.0 + 12.8 * 15) / 100.0
    @test sol_c1[hfm_compiled.MCHERB][end] ≈ expected_mcherb_c1 atol=0.02
end

@testitem "WoodyFuelMoisture Climate Parameters" setup=[NFDRSSetup] tags=[:nfdrs] begin
    # Verify the implementation produces expected values for different climate classes
    # Based on Cohen & Deeming (1985), page 8-9
    wfm = WoodyFuelMoisture()
    wfm_compiled = mtkcompile(wfm)

    # Test frozen stage for climate class 2 (PREGRN = 60%)
    prob_frozen = ODEProblem(wfm_compiled,
        Dict(wfm_compiled.MC1000 => 0.20,
             wfm_compiled.CLIMAT => 2.0,
             wfm_compiled.GRNDAY => 0.0,
             wfm_compiled.is_greenup => 0.0,
             wfm_compiled.is_frozen => 1.0),
        (0.0, 1.0))
    sol_frozen = solve(prob_frozen)
    @test sol_frozen[wfm_compiled.MCWOOD][end] ≈ 0.60 atol=0.01

    # Test green stage for climate class 2 at MC1000 = 20%
    # MCWOOD = (-5.0 + 8.2 * 20) / 100 = 1.59 (159%)
    prob_green = ODEProblem(wfm_compiled,
        Dict(wfm_compiled.MC1000 => 0.20,
             wfm_compiled.CLIMAT => 2.0,
             wfm_compiled.GRNDAY => 0.0,
             wfm_compiled.is_greenup => 0.0,
             wfm_compiled.is_frozen => 0.0),
        (0.0, 1.0))
    sol_green = solve(prob_green)

    expected_mcwood = (-5.0 + 8.2 * 20) / 100.0
    @test sol_green[wfm_compiled.MCWOOD][end] ≈ expected_mcwood atol=0.02
end

@testitem "SpreadComponent Structure" setup=[NFDRSSetup] tags=[:nfdrs] begin
    # Test that SpreadComponent has expected structure
    sc = SpreadComponent()
    @test length(equations(sc)) > 0
    @test length(parameters(sc)) > 0

    # Verify outputs include SC, ROS, IR
    vars = unknowns(sc)
    var_names = [string(Symbolics.tosymbol(v, escape=false)) for v in vars]
    @test "SC" in var_names
    @test "ROS" in var_names
    @test "IR" in var_names
end

@testitem "EnergyReleaseComponent Structure" setup=[NFDRSSetup] tags=[:nfdrs] begin
    # Test that EnergyReleaseComponent has expected structure
    erc = EnergyReleaseComponent()
    @test length(equations(erc)) > 0
    @test length(parameters(erc)) > 0

    # Verify outputs include ERC, IRE
    vars = unknowns(erc)
    var_names = [string(Symbolics.tosymbol(v, escape=false)) for v in vars]
    @test "ERC" in var_names
    @test "IRE" in var_names
end

@testitem "Fuel Loading Conversion" setup=[NFDRSSetup] tags=[:nfdrs] begin
    # Test new SI conversion function (tons/acre to kg/m²)
    @test WildlandFire.fuel_loading_to_kg_per_sqm(1.0) ≈ TONS_PER_ACRE_TO_KG_PER_M2 atol=1e-6

    # Test deprecated conversion (tons/acre to lb/ft²) still works for backward compatibility
    @test WildlandFire.fuel_loading_to_lb_per_sqft(1.0) ≈ 0.0459137 atol=1e-6

    # Verify: 0.20 tons/acre (Model A) converts correctly
    @test WildlandFire.fuel_loading_to_kg_per_sqm(0.20) ≈ 0.20 * TONS_PER_ACRE_TO_KG_PER_M2 atol=1e-6
    @test WildlandFire.fuel_loading_to_lb_per_sqft(0.20) ≈ 0.00918274 atol=1e-6
end

@testitem "LightningFireOccurrenceIndex Structure" setup=[NFDRSSetup] tags=[:nfdrs] begin
    # Test component creation
    loi = LightningFireOccurrenceIndex()
    @test length(equations(loi)) > 0

    # Verify outputs include expected variables
    vars = unknowns(loi)
    var_names = [string(Symbolics.tosymbol(v, escape=false)) for v in vars]
    @test "LOI" in var_names
    @test "LRISK" in var_names
    @test "CGRATE" in var_names
end

@testitem "IgnitionComponent Structure" setup=[NFDRSSetup] tags=[:nfdrs] begin
    # Test component creation
    ic = IgnitionComponent()
    @test length(equations(ic)) > 0

    # Verify outputs include expected variables
    vars = unknowns(ic)
    var_names = [string(Symbolics.tosymbol(v, escape=false)) for v in vars]
    @test "IC" in var_names
    @test "PI" in var_names
    @test "QIGN" in var_names
end

@testitem "HundredHourFuelMoisture Structure" setup=[NFDRSSetup] tags=[:nfdrs] begin
    # Test component creation
    fm100 = HundredHourFuelMoisture()
    @test length(equations(fm100)) > 0

    # Verify outputs include expected variables
    vars = unknowns(fm100)
    var_names = [string(Symbolics.tosymbol(v, escape=false)) for v in vars]
    @test "MC100" in var_names
    @test "EMCBAR" in var_names
    @test "BNDRYH" in var_names
end

@testitem "ThousandHourFuelMoisture Structure" setup=[NFDRSSetup] tags=[:nfdrs] begin
    # Test component creation
    fm1000 = ThousandHourFuelMoisture()
    @test length(equations(fm1000)) > 0

    # Verify outputs include expected variables
    vars = unknowns(fm1000)
    var_names = [string(Symbolics.tosymbol(v, escape=false)) for v in vars]
    @test "MC1000" in var_names
    @test "BNDRYT" in var_names
end
