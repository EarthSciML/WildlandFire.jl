@testsnippet CouplingSetup begin
    using WildlandFire
    using ModelingToolkit
    using ModelingToolkit: t
    using EarthSciMLBase
    using DynamicQuantities
    using Symbolics
    using Test
end

@testitem "FuelModelLookup component" setup = [CouplingSetup] tags = [:coupling] begin
    fm = FuelModelLookup()
    @test fm isa ModelingToolkit.AbstractSystem
    @test length(equations(fm)) == 5
    eq_names = [Symbolics.tosymbol(eq.lhs, escape = false) for eq in equations(fm)]
    @test :σ ∈ eq_names
    @test :w0 ∈ eq_names
    @test :δ ∈ eq_names
    @test :Mx ∈ eq_names
    @test :h ∈ eq_names
end

@testitem "TerrainSlope component" setup = [CouplingSetup] tags = [:coupling] begin
    ts = TerrainSlope()
    @test ts isa ModelingToolkit.AbstractSystem
    @test length(equations(ts)) == 2
    eq_names = [Symbolics.tosymbol(eq.lhs, escape = false) for eq in equations(ts)]
    @test :tanϕ ∈ eq_names
    @test :slope_aspect ∈ eq_names
end

@testitem "MidflameWind component" setup = [CouplingSetup] tags = [:coupling] begin
    mw = MidflameWind()
    @test mw isa ModelingToolkit.AbstractSystem
    @test length(equations(mw)) == 2
    eq_names = [Symbolics.tosymbol(eq.lhs, escape = false) for eq in equations(mw)]
    @test :U ∈ eq_names
    @test :ω ∈ eq_names
end

@testitem "Fuel lookup functions" setup = [CouplingSetup] tags = [:coupling] begin
    # Fuel model 1 (short grass)
    @test WildlandFire.fuel_savr(1.0) ≈ 3500.0 / 0.3048  rtol = 0.01
    @test WildlandFire.fuel_load(1.0) ≈ 0.166  rtol = 0.01
    @test WildlandFire.fuel_depth(1.0) ≈ 0.305  rtol = 0.01
    @test WildlandFire.fuel_mce(1.0) ≈ 0.12  rtol = 0.01
    # Non-burnable code should return default
    @test WildlandFire.fuel_savr(98.0) ≈ 3500.0 / 0.3048  rtol = 0.01
end

@testitem "RothermelFireSpread has CoupleType" setup = [CouplingSetup] tags = [:coupling] begin
    r = RothermelFireSpread()
    ct = ModelingToolkit.getmetadata(r, EarthSciMLBase.CoupleType, nothing)
    @test ct === WildlandFire.RothermelCoupler
end

@testitem "FuelModelLookup-Rothermel coupling" setup = [CouplingSetup] tags = [:coupling] begin
    fm = FuelModelLookup()
    r = RothermelFireSpread()
    cs = couple(fm, r)
    sys = convert(System, cs; compile = false)
    @test sys isa ModelingToolkit.AbstractSystem
    # The coupled system should have equations from both components plus connectors
    all_eqs = equations(sys)
    @test length(all_eqs) > 26 + 5  # Rothermel (26) + FuelModelLookup (5) + connectors
end

@testitem "TerrainSlope-Rothermel coupling" setup = [CouplingSetup] tags = [:coupling] begin
    ts = TerrainSlope()
    r = RothermelFireSpread()
    cs = couple(ts, r)
    sys = convert(System, cs)
    @test sys isa ModelingToolkit.AbstractSystem
end

@testitem "MidflameWind-Rothermel coupling" setup = [CouplingSetup] tags = [:coupling] begin
    mw = MidflameWind()
    r = RothermelFireSpread()
    cs = couple(mw, r)
    sys = convert(System, cs)
    @test sys isa ModelingToolkit.AbstractSystem
end

@testitem "TerrainSlope-MidflameWind coupling" setup = [CouplingSetup] tags = [:coupling] begin
    ts = TerrainSlope()
    mw = MidflameWind()
    cs = couple(ts, mw)
    sys = convert(System, cs)
    @test sys isa ModelingToolkit.AbstractSystem
end

@testitem "EMC has CoupleType" setup = [CouplingSetup] tags = [:coupling] begin
    emc = EquilibriumMoistureContent()
    ct = ModelingToolkit.getmetadata(emc, EarthSciMLBase.CoupleType, nothing)
    @test ct === WildlandFire.EMCCoupler
end

@testitem "OneHourFuelMoisture has CoupleType" setup = [CouplingSetup] tags = [:coupling] begin
    fm1 = OneHourFuelMoisture()
    ct = ModelingToolkit.getmetadata(fm1, EarthSciMLBase.CoupleType, nothing)
    @test ct === WildlandFire.OneHourFuelMoistureCoupler
end

@testitem "EMC-OneHourFuelMoisture coupling" setup = [CouplingSetup] tags = [:coupling] begin
    emc = EquilibriumMoistureContent()
    fm1 = OneHourFuelMoisture()
    cs = couple(emc, fm1)
    sys = convert(System, cs; compile = false)
    @test sys isa ModelingToolkit.AbstractSystem
end

@testitem "OneHourFuelMoisture-Rothermel coupling" setup = [CouplingSetup] tags = [:coupling] begin
    fm1 = OneHourFuelMoisture()
    r = RothermelFireSpread()
    cs = couple(fm1, r)
    sys = convert(System, cs; compile = false)
    @test sys isa ModelingToolkit.AbstractSystem
end

@testitem "Full fire model coupling" setup = [CouplingSetup] tags = [:coupling] begin
    fm = FuelModelLookup()
    ts = TerrainSlope()
    mw = MidflameWind()
    emc = EquilibriumMoistureContent()
    fm1 = OneHourFuelMoisture()
    r = RothermelFireSpread()
    cs = couple(fm, ts, mw, emc, fm1, r)
    sys = convert(System, cs; compile = false)
    @test sys isa ModelingToolkit.AbstractSystem
    # Should have all equations from all components plus connectors
    # Rothermel(26) + FuelModelLookup(5) + TerrainSlope(2) + MidflameWind(2) + EMC(1) + OneHourFM(1) + connectors
    @test length(equations(sys)) > 26 + 5 + 2 + 2 + 1 + 1
end
