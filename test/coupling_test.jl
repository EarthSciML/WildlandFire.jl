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

@testitem "USGS3DEP-TerrainSlope coupling" setup = [CouplingSetup] tags = [:coupling] begin
    using EarthSciData
    using Dates

    domain = DomainInfo(
        DateTime(2018, 11, 8), DateTime(2018, 11, 9);
        lonrange = deg2rad(-121.65):deg2rad(0.01):deg2rad(-121.55),
        latrange = deg2rad(39.73):deg2rad(0.01):deg2rad(39.83),
        levrange = 1:1,
    )
    dep = USGS3DEP(domain; resolution = 10.0)
    ts = TerrainSlope()

    cs = couple(dep, ts)
    sys = convert(System, cs; compile = false)
    @test sys isa ModelingToolkit.AbstractSystem
    # Should have USGS3DEP equations (3) + TerrainSlope equations (2) + connectors (2)
    @test length(equations(sys)) >= 3 + 2 + 2
end

@testitem "Full fire model coupling" setup = [CouplingSetup] tags = [:coupling] begin
    fm = FuelModelLookup()
    ts = TerrainSlope()
    mw = MidflameWind()
    emc = EquilibriumMoistureContent()
    fm1 = OneHourFuelMoisture()
    fc = FuelConsumption()
    r = RothermelFireSpread()
    cs = couple(fm, ts, mw, emc, fm1, fc, r)
    sys = convert(System, cs; compile = false)
    @test sys isa ModelingToolkit.AbstractSystem
    # Should have all equations from all components plus connectors
    # Rothermel(26) + FuelModelLookup(5) + TerrainSlope(2) + MidflameWind(2) + EMC(1) + OneHourFM(1) + FuelConsumption(2) + connectors
    @test length(equations(sys)) > 26 + 5 + 2 + 2 + 1 + 1 + 2
end

@testitem "LevelSetFireSpread has CoupleType" setup = [CouplingSetup] tags = [:coupling] begin
    using DomainSets

    @parameters x [unit = u"m"]
    @parameters y [unit = u"m"]
    domain = DomainInfo(
        constIC(0.0, t ∈ Interval(0.0, 60.0)),
        constBC(0.0, x ∈ Interval(0.0, 500.0), y ∈ Interval(0.0, 500.0)),
    )
    ls = LevelSetFireSpread(
        domain;
        initial_condition = (x, y) -> sqrt((x - 250.0)^2 + (y - 250.0)^2) - 10.0,
    )
    @test ls.metadata[EarthSciMLBase.CoupleType] === WildlandFire.LevelSetCoupler
end

@testitem "Rothermel-LevelSet coupling" setup = [CouplingSetup] tags = [:coupling] begin
    using DomainSets

    r = RothermelFireSpread()

    @parameters x [unit = u"m"]
    @parameters y [unit = u"m"]
    domain = DomainInfo(
        constIC(0.0, t ∈ Interval(0.0, 60.0)),
        constBC(0.0, x ∈ Interval(0.0, 500.0), y ∈ Interval(0.0, 500.0)),
    )
    ls = LevelSetFireSpread(
        domain;
        initial_condition = (x, y) -> sqrt((x - 250.0)^2 + (y - 250.0)^2) - 10.0,
    )

    cs = EarthSciMLBase.couple2(
        WildlandFire.RothermelCoupler(r),
        WildlandFire.LevelSetCoupler(ls),
    )

    # The coupling should produce a ConnectorSystem
    @test cs isa EarthSciMLBase.ConnectorSystem

    # The connector equation should link S to R
    @test length(cs.eqs) == 1
    eq = cs.eqs[1]
    lhs_name = Symbolics.tosymbol(eq.lhs, escape = false)
    @test lhs_name == :S

    # The modified level-set system should no longer have S as a parameter
    @test !any(p -> Symbol(p) == :S, cs.from.ps)

    # The Rothermel system should be unchanged
    @test cs.to isa ModelingToolkit.AbstractSystem
end

@testitem "Rothermel-LevelSet couple and convert" setup = [CouplingSetup] tags = [:coupling] begin
    using DomainSets

    r = RothermelFireSpread()

    @parameters x [unit = u"m"]
    @parameters y [unit = u"m"]
    domain = DomainInfo(
        constIC(0.0, t ∈ Interval(0.0, 10.0)),
        constBC(0.0, x ∈ Interval(0.0, 100.0), y ∈ Interval(0.0, 100.0)),
    )
    ls = LevelSetFireSpread(
        domain;
        initial_condition = (x, y) -> sqrt((x - 50.0)^2 + (y - 50.0)^2) - 10.0,
    )

    # Couple using the EarthSciMLBase.couple function
    cs = couple(r, ls, domain)

    # CoupledSystem should contain both systems
    @test length(cs.systems) == 1       # Rothermel (ODE)
    @test length(cs.pdesystems) == 1    # LevelSet (PDE)

    # Convert to a merged PDESystem
    pde = convert(PDESystem, cs)
    @test pde isa PDESystem

    # The merged PDE should have 27 equations: 1 level-set PDE + 26 Rothermel algebraic
    @test length(equations(pde)) == 27

    # The level-set ψ should be a dependent variable
    dv_names = [Symbolics.tosymbol(dv, escape = false) for dv in pde.dvs]
    @test :ψ ∈ dv_names

    # Rothermel R should be a dependent variable (promoted to spatial)
    @test any(n -> occursin("R", string(n)), string.(dv_names))

    # Discretize and solve on a coarse grid.
    # Parameter defaults from @constants metadata must be copied into
    # initial_conditions for MethodOfLines to find them.
    using MethodOfLines, OrdinaryDiffEqSSPRK
    for p in pde.ps
        if ModelingToolkit.hasdefault(p)
            pde.initial_conditions[Symbolics.unwrap(p)] = ModelingToolkit.getdefault(p)
        end
    end
    # Set Rothermel fuel model 1 (short grass) inputs
    for p in pde.ps
        sym = Symbolics.tosymbol(p, escape = false)
        if sym == Symbol("RothermelFireSpread₊σ")
            pde.initial_conditions[Symbolics.unwrap(p)] = 11483.0
        elseif sym == Symbol("RothermelFireSpread₊w0")
            pde.initial_conditions[Symbolics.unwrap(p)] = 0.166
        elseif sym == Symbol("RothermelFireSpread₊δ")
            pde.initial_conditions[Symbolics.unwrap(p)] = 0.3048
        elseif sym == Symbol("RothermelFireSpread₊Mx")
            pde.initial_conditions[Symbolics.unwrap(p)] = 0.12
        elseif sym == Symbol("RothermelFireSpread₊Mf")
            pde.initial_conditions[Symbolics.unwrap(p)] = 0.05
        elseif sym == Symbol("RothermelFireSpread₊U")
            pde.initial_conditions[Symbolics.unwrap(p)] = 2.235
        elseif sym == Symbol("RothermelFireSpread₊tanϕ")
            pde.initial_conditions[Symbolics.unwrap(p)] = 0.0
        end
    end

    dx = 25.0
    discretization = MOLFiniteDifference(
        [pde.ivs[2] => dx, pde.ivs[3] => dx], pde.ivs[1],
    )
    prob = MethodOfLines.discretize(pde, discretization; checks = false)
    sol = solve(prob, SSPRK33(); dt = 0.5, adaptive = false, saveat = 10.0)
    @test sol.retcode == SciMLBase.ReturnCode.Success
end

@testitem "FuelConsumption has CoupleType" setup = [CouplingSetup] tags = [:coupling] begin
    fc = FuelConsumption()
    ct = ModelingToolkit.getmetadata(fc, EarthSciMLBase.CoupleType, nothing)
    @test ct === WildlandFire.FuelConsumptionCoupler
end

@testitem "LevelSet-FuelConsumption coupling" setup = [CouplingSetup] tags = [:coupling] begin
    using DomainSets

    @parameters x [unit = u"m"]
    @parameters y [unit = u"m"]
    domain = DomainInfo(
        constIC(0.0, t ∈ Interval(0.0, 60.0)),
        constBC(0.0, x ∈ Interval(0.0, 500.0), y ∈ Interval(0.0, 500.0)),
    )
    ls = LevelSetFireSpread(
        domain;
        initial_condition = (x, y) -> sqrt((x - 250.0)^2 + (y - 250.0)^2) - 10.0,
    )
    fc = FuelConsumption()

    cs = EarthSciMLBase.couple2(
        WildlandFire.LevelSetCoupler(ls),
        WildlandFire.FuelConsumptionCoupler(fc),
    )

    # The coupling should produce a ConnectorSystem
    @test cs isa EarthSciMLBase.ConnectorSystem

    # The connector equation should set is_burning via smooth Heaviside
    @test length(cs.eqs) == 1
    eq = cs.eqs[1]
    lhs_name = string(Symbolics.tosymbol(eq.lhs, escape = false))
    @test occursin("is_burning", lhs_name)

    # is_burning should have been promoted from parameter to variable
    param_names = [Symbolics.tosymbol(p, escape = false) for p in parameters(cs.from)]
    @test :is_burning ∉ param_names
end

@testitem "FuelConsumption-Rothermel coupling" setup = [CouplingSetup] tags = [:coupling] begin
    fc = FuelConsumption()
    r = RothermelFireSpread()

    cs = EarthSciMLBase.couple2(
        WildlandFire.FuelConsumptionCoupler(fc),
        WildlandFire.RothermelCoupler(r),
    )

    # The coupling should produce a ConnectorSystem
    @test cs isa EarthSciMLBase.ConnectorSystem

    # The connector equation should link w0 to w0_effective
    @test length(cs.eqs) == 1
    eq = cs.eqs[1]
    lhs_name = string(Symbolics.tosymbol(eq.lhs, escape = false))
    @test occursin("w0", lhs_name)

    # w0 should have been promoted from parameter to variable in Rothermel
    param_names = [Symbolics.tosymbol(p, escape = false) for p in parameters(cs.from)]
    @test :w0 ∉ param_names
end

@testitem "FuelModelLookup-FuelConsumption coupling" setup = [CouplingSetup] tags = [:coupling] begin
    fm = FuelModelLookup()
    fc = FuelConsumption()

    cs = EarthSciMLBase.couple2(
        WildlandFire.FuelModelLookupCoupler(fm),
        WildlandFire.FuelConsumptionCoupler(fc),
    )

    @test cs isa EarthSciMLBase.ConnectorSystem
    @test length(cs.eqs) == 1
    eq = cs.eqs[1]
    lhs_name = string(Symbolics.tosymbol(eq.lhs, escape = false))
    @test occursin("w0_initial", lhs_name)

    # w0_initial should have been promoted from parameter to variable
    param_names = [Symbolics.tosymbol(p, escape = false) for p in parameters(cs.from)]
    @test :w0_initial ∉ param_names
end
