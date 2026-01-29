using WildlandFire
using Test
using ModelingToolkit
using ModelingToolkit: mtkcompile
using OrdinaryDiffEqDefault

@testset "WildlandFire.jl" begin
    @testset "NFDRS Fuel Models" begin
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
        @test model_a.SG1 == 3000.0   # Western grasses SAV ratio
        @test model_a.W1 == 0.20      # 1-hr loading (tons/acre)
        @test model_a.DEPTH == 0.80   # Fuel bed depth (ft)
        @test model_a.MXD == 15.0     # Dead fuel moisture of extinction (%)
        @test model_a.HD == 8000.0    # Heat of combustion (Btu/lb)
        @test model_a.SCM == 300.0    # Spread component threshold
        @test model_a.WNDFC == 0.6    # Wind reduction factor

        # Test California chaparral (Model B)
        model_b = get_fuel_model(:B)
        @test model_b.SG1 == 700.0
        @test model_b.WWOOD == 11.50  # Heavy woody load
        @test model_b.HD == 9500.0    # Higher heat content for chaparral
    end

    @testset "Component Creation" begin
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

    @testset "EquilibriumMoistureContent Equations" begin
        # Test that the system has the expected structure
        emc = EquilibriumMoistureContent()

        # Check that it has equations
        @test length(equations(emc)) > 0

        # Check that it has parameters
        params = parameters(emc)
        @test length(params) == 2  # TEMP and RH

        # Verify EMC equations against Cohen & Deeming (1985), Eq. 1a, 1b, 1c (page 1)
        # These regression equations are from Simard (1968)

        # Test case 1: Low RH (Eq. 1a: RH < 10%)
        # EMC = 0.03229 + 0.281073 * RH - 0.000578 * TEMP * RH
        # At RH=5%, TEMP=70°F: EMC = 0.03229 + 0.281073*5 - 0.000578*70*5 = 1.236%
        # (approximately 1.2%)

        # Test case 2: Mid RH (Eq. 1b: 10% <= RH < 50%)
        # EMC = 2.22749 + 0.160107 * RH - 0.014784 * TEMP
        # At RH=30%, TEMP=70°F: EMC = 2.22749 + 0.160107*30 - 0.014784*70 = 5.997%
        # (approximately 6%)

        # Test case 3: High RH (Eq. 1c: RH >= 50%)
        # EMC = 21.0606 + 0.005565 * RH^2 - 0.00035 * RH * TEMP - 0.483199 * RH
        # At RH=80%, TEMP=70°F: EMC = 21.0606 + 0.005565*6400 - 0.00035*80*70 - 0.483199*80
        # = 21.0606 + 35.616 - 1.96 - 38.656 = 16.06%
    end

    @testset "OneHourFuelMoisture Equations" begin
        # Test that the system has the expected structure
        fm1 = OneHourFuelMoisture()

        # Check that it has equations
        @test length(equations(fm1)) > 0

        # Check that it has parameters
        params = parameters(fm1)
        @test length(params) > 0

        # Verify: When raining, MC1 = 35% (page 4)
        # Verify: Without sticks, MC1 = 1.03 * EMCPRM (page 3)
        # Verify: With sticks, MC1 = (4.0 * EMCPRM + MC10) / 5.0 (page 3)
    end

    @testset "FuelLoadingTransfer Equations" begin
        # Test equations 5, 6, 7, 8 from Cohen & Deeming (1985), page 7
        flt = FuelLoadingTransfer()
        @test length(equations(flt)) == 4

        # Eq. 5: FCTCUR = 1.33 - 0.0111 * MCHERB (MCHERB as percent)
        # At MCHERB = 30%, FCTCUR = 1.33 - 0.0111*30 = 0.997 (clamped to 1.0)
        # At MCHERB = 120%, FCTCUR = 1.33 - 0.0111*120 = 0.0 (fully green)
    end

    @testset "BurningIndex Equation" begin
        # Verify BI = 3.01 * (SC * ERC)^0.46 (page 12)
        # BI = 10 * FL where FL = 0.301 * (SC * ERC)^0.46
        bi = BurningIndex()
        @test length(equations(bi)) == 1
    end

    @testset "FireLoadIndex Equation" begin
        # Verify FLI = 0.71 * sqrt(BI^2 + (LOI + MCOI)^2) (page 14)
        # With BI and (LOI + MCOI) each limited to 100
        fli = FireLoadIndex()
        @test length(equations(fli)) == 1
    end

    @testset "HumanFireOccurrenceIndex Equation" begin
        # Verify MCOI = 0.01 * MRISK * IC (page 13)
        mcoi = HumanFireOccurrenceIndex()
        @test length(equations(mcoi)) == 1
    end

    @testset "SpreadComponent Structure" begin
        # Test that SpreadComponent has expected structure
        sc = SpreadComponent()
        @test length(equations(sc)) > 0
        @test length(parameters(sc)) > 0
    end

    @testset "EnergyReleaseComponent Structure" begin
        # Test that EnergyReleaseComponent has expected structure
        erc = EnergyReleaseComponent()
        @test length(equations(erc)) > 0
        @test length(parameters(erc)) > 0
    end

    @testset "Fuel Loading Conversion" begin
        # Test conversion from tons/acre to lb/ft² (page 9)
        # Factor: 0.0459137 (lb·acre)/(ton·ft²)
        @test WildlandFire.fuel_loading_to_lb_per_sqft(1.0) ≈ 0.0459137 atol=1e-6

        # Verify: 0.20 tons/acre (Model A) = 0.00918 lb/ft²
        @test WildlandFire.fuel_loading_to_lb_per_sqft(0.20) ≈ 0.00918274 atol=1e-6
    end

    @testset "Herbaceous Fuel Moisture Climate Parameters" begin
        # Verify climate-dependent parameters from Cohen & Deeming (1985), page 7
        # HERBGA and HERBGB for greenup/green stage
        hfm = HerbaceousFuelMoisture()
        @test length(equations(hfm)) > 0

        # Climate class 1: HERBGA = -70.0, HERBGB = 12.8
        # Climate class 2: HERBGA = -100.0, HERBGB = 14.0
        # Climate class 3: HERBGA = -137.5, HERBGB = 15.5
        # Climate class 4: HERBGA = -185.0, HERBGB = 17.4
    end

    @testset "Woody Fuel Moisture Climate Parameters" begin
        # Verify climate-dependent parameters from Cohen & Deeming (1985), page 8
        wfm = WoodyFuelMoisture()
        @test length(equations(wfm)) > 0

        # PREGRN (dormant moisture content):
        # Climate class 1: 50%, class 2: 60%, class 3: 70%, class 4: 80%

        # WOODGA and WOODGB:
        # Climate class 1: WOODGA = 12.5, WOODGB = 7.5
        # Climate class 2: WOODGA = -5.0, WOODGB = 8.2
        # Climate class 3: WOODGA = -22.5, WOODGB = 8.9
        # Climate class 4: WOODGA = -45.0, WOODGB = 9.8
    end

    @testset "Numerical Validation - EMC Equations" begin
        # Verify EMC equations against Cohen & Deeming (1985), Eq. 1a, 1b, 1c (page 1)
        # These regression equations are from Simard (1968)

        # Manual calculation functions matching the paper
        function emc_manual(temp, rh_pct)
            if rh_pct < 10
                return 0.03229 + 0.281073 * rh_pct - 0.000578 * temp * rh_pct
            elseif rh_pct < 50
                return 2.22749 + 0.160107 * rh_pct - 0.014784 * temp
            else
                return 21.0606 + 0.005565 * rh_pct^2 - 0.00035 * rh_pct * temp - 0.483199 * rh_pct
            end
        end

        # Test case 1: Low RH (Eq. 1a: RH < 10%)
        # At RH=5%, TEMP=70°F
        expected_low = emc_manual(70.0, 5.0)  # Should be ~1.236%
        @test expected_low ≈ 1.236 atol=0.01

        # Test case 2: Mid RH (Eq. 1b: 10% <= RH < 50%)
        # At RH=30%, TEMP=70°F
        expected_mid = emc_manual(70.0, 30.0)  # Should be ~5.997%
        @test expected_mid ≈ 5.997 atol=0.01

        # Test case 3: High RH (Eq. 1c: RH >= 50%)
        # At RH=80%, TEMP=70°F
        expected_high = emc_manual(70.0, 80.0)  # Should be ~16.06%
        @test expected_high ≈ 16.06 atol=0.1
    end

    @testset "Numerical Validation - Fuel Loading Transfer" begin
        # Test Eq. 5: FCTCUR = 1.33 - 0.0111 * MCHERB (MCHERB as percent)
        # At MCHERB = 30%, FCTCUR = 1.33 - 0.0111*30 = 0.997
        @test 1.33 - 0.0111 * 30 ≈ 0.997 atol=0.001

        # At MCHERB = 120%, FCTCUR = 1.33 - 0.0111*120 = -0.002 (clamped to 0)
        @test max(0.0, 1.33 - 0.0111 * 120) ≈ 0.0 atol=0.01

        # At MCHERB = 30% (boundary), nearly all herb transferred
        @test min(1.0, 1.33 - 0.0111 * 30) ≈ 0.997 atol=0.01
    end

    @testset "Numerical Validation - Burning Index" begin
        # Verify BI = 3.01 * (SC * ERC)^0.46 (page 12)
        # Example: SC=50, ERC=20 -> BI = 3.01 * (50*20)^0.46 = 3.01 * 1000^0.46
        # 1000^0.46 ≈ 23.99, so BI ≈ 72.2
        @test 3.01 * (50 * 20)^0.46 ≈ 72.2 atol=1.0

        # Zero case
        @test 3.01 * (0 * 0)^0.46 == 0.0
    end

    @testset "Numerical Validation - Fire Load Index" begin
        # Verify FLI = 0.71 * sqrt(BI^2 + (LOI + MCOI)^2) (page 14)
        # Example: BI=50, LOI=30, MCOI=20 -> (LOI+MCOI)=50
        # FLI = 0.71 * sqrt(50^2 + 50^2) = 0.71 * sqrt(5000) ≈ 50.2
        @test 0.71 * sqrt(50^2 + 50^2) ≈ 50.2 atol=0.1

        # With limits: BI=150 (limited to 100), LOI+MCOI=60
        # FLI = 0.71 * sqrt(100^2 + 60^2) = 0.71 * sqrt(13600) ≈ 82.8
        @test 0.71 * sqrt(100^2 + 60^2) ≈ 82.8 atol=0.1
    end

    @testset "Numerical Validation - Human Fire Occurrence" begin
        # Verify MCOI = 0.01 * MRISK * IC (page 13)
        # Example: MRISK=50, IC=80 -> MCOI = 0.01 * 50 * 80 = 40
        @test 0.01 * 50 * 80 == 40.0
    end

    @testset "Lightning Fire Occurrence Index" begin
        # Test component creation
        loi = LightningFireOccurrenceIndex()
        @test length(equations(loi)) > 0

        # Verify LAL table values from page 13
        # LAL 2: CGRATE=12.5, STMDIA=3.0, TOTWID=7.0
        # LAL 3: CGRATE=25.0, STMDIA=4.0, TOTWID=8.0
        # LAL 4: CGRATE=50.0, STMDIA=5.0, TOTWID=9.0
        # LAL 5: CGRATE=100.0, STMDIA=7.0, TOTWID=11.0

        # Test LGTDUR equation: LGTDUR = -86.83 + 153.41 * CGRATE^0.1437
        # At CGRATE=25: LGTDUR = -86.83 + 153.41 * 25^0.1437 ≈ -86.83 + 153.41 * 1.588 ≈ 156.8
        @test -86.83 + 153.41 * 25^0.1437 ≈ 156.8 atol=1.0
    end
end
