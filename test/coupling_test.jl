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
    r = RothermelFireSpread()
    fsd = FireSpreadDirection()
    cs = couple(fm, ts, mw, emc, fm1, r, fsd)
    sys = convert(System, cs; compile = false)
    @test sys isa ModelingToolkit.AbstractSystem
    # Should have all equations from all components plus connectors
    # Rothermel(26) + FuelModelLookup(5) + TerrainSlope(2) + MidflameWind(2)
    # + EMC(1) + OneHourFM(1) + FireSpreadDirection(6) + connectors
    @test length(equations(sys)) > 26 + 5 + 2 + 2 + 1 + 1 + 6
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

@testitem "FireSpreadDirection has CoupleType" setup = [CouplingSetup] tags = [:coupling] begin
    fsd = FireSpreadDirection()
    ct = ModelingToolkit.getmetadata(fsd, EarthSciMLBase.CoupleType, nothing)
    @test ct === WildlandFire.FireSpreadDirectionCoupler
end

@testitem "Rothermel-FireSpreadDirection coupling" setup = [CouplingSetup] tags = [:coupling] begin
    r = RothermelFireSpread()
    fsd = FireSpreadDirection()
    cs = EarthSciMLBase.couple2(
        WildlandFire.RothermelCoupler(r),
        WildlandFire.FireSpreadDirectionCoupler(fsd),
    )

    @test cs isa EarthSciMLBase.ConnectorSystem

    # Should couple 7 parameters: R0, φw, φs, β_ratio, C_coeff, B_coeff, E_coeff
    @test length(cs.eqs) == 7

    eq_lhs_names = Set(string(Symbolics.tosymbol(eq.lhs, escape = false)) for eq in cs.eqs)
    for name in ["R0", "φw", "φs", "β_ratio", "C_coeff", "B_coeff", "E_coeff"]
        @test any(n -> endswith(n, name), eq_lhs_names)
    end
end

@testitem "MidflameWind-FireSpreadDirection coupling" setup = [CouplingSetup] tags = [:coupling] begin
    mw = MidflameWind()
    fsd = FireSpreadDirection()
    cs = EarthSciMLBase.couple2(
        WildlandFire.MidflameWindCoupler(mw),
        WildlandFire.FireSpreadDirectionCoupler(fsd),
    )

    @test cs isa EarthSciMLBase.ConnectorSystem
    @test length(cs.eqs) == 1
    lhs_name = string(Symbolics.tosymbol(cs.eqs[1].lhs, escape = false))
    @test endswith(lhs_name, "ω")
end

@testitem "FireSpreadDirection-LevelSet coupling" setup = [CouplingSetup] tags = [:coupling] begin
    using DomainSets

    fsd = FireSpreadDirection()

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
        WildlandFire.FireSpreadDirectionCoupler(fsd),
        WildlandFire.LevelSetCoupler(ls),
    )

    @test cs isa EarthSciMLBase.ConnectorSystem

    # Should couple 3 parameters: R_H, Z, α
    @test length(cs.eqs) == 3

    eq_lhs_names = Set(Symbolics.tosymbol(eq.lhs, escape = false) for eq in cs.eqs)
    for name in [:R_H, :Z, :α]
        @test name ∈ eq_lhs_names
    end

    # The level-set system should no longer have R_H, Z, α as parameters
    @test !any(p -> Symbol(p) == :R_H, cs.from.ps)
    @test !any(p -> Symbol(p) == :Z, cs.from.ps)
    @test !any(p -> Symbol(p) == :α, cs.from.ps)
end

@testitem "FireSpreadDirection-LevelSet couple and convert" setup = [CouplingSetup] tags = [:coupling] begin
    using DomainSets

    r = RothermelFireSpread()
    mw = MidflameWind()
    fsd = FireSpreadDirection()

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
    cs = couple(r, mw, fsd, ls, domain)

    # CoupledSystem should contain both ODE systems and PDE systems
    @test length(cs.systems) >= 1
    @test length(cs.pdesystems) == 1    # LevelSet (PDE)

    # Convert to a merged PDESystem
    pde = convert(PDESystem, cs)
    @test pde isa PDESystem

    # The level-set ψ should be a dependent variable
    dv_names = [Symbolics.tosymbol(dv, escape = false) for dv in pde.dvs]
    @test :ψ ∈ dv_names

    # Deduplicate equations: the EarthSciMLBase coupling loop fires couple2 for
    # both (i,j) and (j,i) orderings, creating duplicate connector equations.
    # This is harmless for ODE systems (mtkcompile resolves them) but
    # MethodOfLines requires exact equation/unknown matching.
    seen = Set{String}()
    unique_eqs = Symbolics.Equation[]
    for eq in equations(pde)
        s = string(eq)
        if s ∉ seen
            push!(seen, s)
            push!(unique_eqs, eq)
        end
    end
    pde = PDESystem(
        unique_eqs, pde.bcs, pde.domain, pde.ivs, pde.dvs, pde.ps;
        name = pde.name
    )

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
