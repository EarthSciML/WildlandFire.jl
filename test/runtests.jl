using WildlandFire
using Test
using ModelingToolkit
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
    end

    @testset "Component Creation" begin
        # Test that components can be created (basic smoke test)
        @test EquilibriumMoistureContent() isa System
        @test OneHourFuelMoisture() isa System
        @test TenHourFuelMoisture() isa System
        @test HerbaceousFuelMoisture() isa System
        @test WoodyFuelMoisture() isa System
        @test FuelLoadingTransfer() isa System
        @test BurningIndex() isa System
        @test IgnitionComponent() isa System
        @test HumanFireOccurrenceIndex() isa System
        @test FireLoadIndex() isa System
    end

    @testset "OneHourFuelMoisture Equations" begin
        # Test that the system has the expected structure
        fm1 = OneHourFuelMoisture()

        # Check that it has equations
        @test length(equations(fm1)) > 0

        # Check that it has parameters
        params = parameters(fm1)
        @test length(params) > 0
    end

    @testset "EquilibriumMoistureContent Equations" begin
        # Test that the system has the expected structure
        emc = EquilibriumMoistureContent()

        # Check that it has equations
        @test length(equations(emc)) > 0

        # Check that it has parameters
        params = parameters(emc)
        @test length(params) > 0
    end
end
