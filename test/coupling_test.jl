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
    @test length(equations(fm)) == 6
    eq_names = [Symbolics.tosymbol(eq.lhs, escape = false) for eq in equations(fm)]
    @test :σ ∈ eq_names
    @test :w0 ∈ eq_names
    @test :δ ∈ eq_names
    @test :Mx ∈ eq_names
    @test :h ∈ eq_names
    @test :T_f ∈ eq_names
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
    # Fuel weight (T_f = weight / 0.8514)
    @test WildlandFire.fuel_weight(1.0) ≈ 7.0 / 0.8514 rtol = 0.01
    @test WildlandFire.fuel_weight(8.0) ≈ 900.0 / 0.8514 rtol = 0.01
    # Non-burnable code 98 (water): SAVR stays at FM1 default to avoid division by zero
    @test WildlandFire.fuel_savr(98.0) ≈ 3500.0 / 0.3048  rtol = 0.01
    # Heat content should be zero for non-burnable codes
    @test WildlandFire.fuel_heat(98.0) == 0.0
    # Fuel weight for non-burnable codes should return FM1 default
    @test WildlandFire.fuel_weight(98.0) ≈ 7.0 / 0.8514 rtol = 0.01
end

@testitem "Non-burnable fuel lookup" setup = [CouplingSetup] tags = [:coupling] begin
    # All standard non-burnable LANDFIRE FBFM13 codes
    for code in [91.0, 92.0, 93.0, 98.0, 99.0]
        # Heat content must be zero to ensure IR=0 and therefore R=0
        @test WildlandFire.fuel_heat(code) == 0.0
        # SAVR, depth, load, and Mx must remain positive to avoid division by zero in Rothermel
        @test WildlandFire.fuel_savr(code) > 0.0
        @test WildlandFire.fuel_depth(code) > 0.0
        @test WildlandFire.fuel_load(code) > 0.0
        @test WildlandFire.fuel_mce(code) > 0.0
    end

    # Unrecognized codes should also return zero heat content
    @test WildlandFire.fuel_heat(0.0) == 0.0
    @test WildlandFire.fuel_heat(999.0) == 0.0

    # Valid fuel models (1-13) should still return nonzero heat
    @test WildlandFire.fuel_heat(1.0) > 0.0
    @test WildlandFire.fuel_heat(13.0) > 0.0
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
    @test length(all_eqs) > 26 + 6  # Rothermel (26) + FuelModelLookup (6) + connectors
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
    fsd = FireSpreadDirection()
    cs = couple(fm, ts, mw, emc, fm1, fc, r, fsd)
    sys = convert(System, cs; compile = false)
    @test sys isa ModelingToolkit.AbstractSystem
    # Should have all equations from all components plus connectors
    # Rothermel(26) + FuelModelLookup(6) + TerrainSlope(2) + MidflameWind(2)
    # + EMC(1) + OneHourFM(1) + FuelConsumption(2) + FireSpreadDirection(6) + connectors
    @test length(equations(sys)) > 26 + 6 + 2 + 2 + 1 + 1 + 2 + 6
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

@testitem "FuelModelLookup-FuelConsumption coupling" setup = [CouplingSetup] tags = [:coupling] begin
    fm = FuelModelLookup()
    fc = FuelConsumption()

    cs = EarthSciMLBase.couple2(
        WildlandFire.FuelModelLookupCoupler(fm),
        WildlandFire.FuelConsumptionCoupler(fc),
    )

    @test cs isa EarthSciMLBase.ConnectorSystem
    # Should couple both w0_initial and T_f
    @test length(cs.eqs) == 2
    eq_lhs_names = Set(string(Symbolics.tosymbol(eq.lhs, escape = false)) for eq in cs.eqs)
    @test any(n -> occursin("w0_initial", n), eq_lhs_names)
    @test any(n -> occursin("T_f", n), eq_lhs_names)

    # w0_initial and T_f should have been promoted from parameters to variables
    param_names = [Symbolics.tosymbol(p, escape = false) for p in parameters(cs.from)]
    @test :w0_initial ∉ param_names
    @test :T_f ∉ param_names
end

@testitem "EMC-OneHourFuelMoisture numerical verification" setup = [CouplingSetup] tags = [:coupling] begin
    using ModelingToolkit: mtkcompile
    using OrdinaryDiffEqDefault

    emc = EquilibriumMoistureContent()
    fm1 = OneHourFuelMoisture()
    cs = couple(emc, fm1)
    sys = convert(System, cs; compile = false)
    compiled = mtkcompile(sys)

    # TEMP=294.26K (70°F), RH=30%, no fuel sticks, no rain
    # Expected EMC ≈ 0.060 (from Cohen & Deeming 1985, Eq. 1b)
    # Expected MC1 = 1.03 × EMC ≈ 0.062
    temp_K = (70.0 + 459.67) * 5 / 9  # 70°F to K
    prob = ODEProblem(
        compiled,
        Dict(
            compiled.EquilibriumMoistureContent₊TEMP => temp_K,
            compiled.EquilibriumMoistureContent₊RH => 0.3,
            compiled.OneHourFuelMoisture₊use_fuel_sticks => 0.0,
            compiled.OneHourFuelMoisture₊is_raining => 0.0,
            compiled.OneHourFuelMoisture₊MC10 => 0.0,
        ),
        (0.0, 1.0),
    )
    sol = solve(prob)

    emc_val = sol[compiled.EquilibriumMoistureContent₊EMC][end]
    mc1_val = sol[compiled.OneHourFuelMoisture₊MC1][end]

    @test emc_val ≈ 0.06 atol = 0.005
    @test mc1_val ≈ 1.03 * emc_val atol = 0.001
end

@testitem "OneHourFuelMoisture-Rothermel dry conditions" setup = [CouplingSetup] tags = [:coupling] begin
    using ModelingToolkit: mtkcompile
    using OrdinaryDiffEqDefault

    fm1 = OneHourFuelMoisture()
    r = RothermelFireSpread()
    cs = couple(fm1, r)
    sys = convert(System, cs; compile = false)
    compiled = mtkcompile(sys)

    # Dry conditions: EMCPRM=0.04 → MC1=0.0412, well below Mx=0.12
    # Fuel Model 1 (short grass) parameters
    prob = ODEProblem(
        compiled,
        Dict(
            compiled.OneHourFuelMoisture₊EMCPRM => 0.04,
            compiled.OneHourFuelMoisture₊use_fuel_sticks => 0.0,
            compiled.OneHourFuelMoisture₊is_raining => 0.0,
            compiled.OneHourFuelMoisture₊MC10 => 0.0,
            compiled.RothermelFireSpread₊σ => 11483.5,
            compiled.RothermelFireSpread₊w0 => 0.166,
            compiled.RothermelFireSpread₊δ => 0.3048,
            compiled.RothermelFireSpread₊Mx => 0.12,
            compiled.RothermelFireSpread₊h => 18608000.0,
            compiled.RothermelFireSpread₊U => 0.0,
            compiled.RothermelFireSpread₊tanϕ => 0.0,
        ),
        (0.0, 1.0),
    )
    sol = solve(prob)

    mc1_val = sol[compiled.OneHourFuelMoisture₊MC1][end]
    R_val = sol[compiled.RothermelFireSpread₊R][end]
    IR_val = sol[compiled.RothermelFireSpread₊IR][end]

    # MC1 should be 1.03 * 0.04 = 0.0412
    @test mc1_val ≈ 0.0412 atol = 0.001
    # Fire should spread (R > 0, IR > 0) in dry conditions
    @test R_val > 0
    @test IR_val > 0
end

@testitem "OneHourFuelMoisture-Rothermel moisture extinction" setup = [CouplingSetup] tags = [:coupling] begin
    using ModelingToolkit: mtkcompile
    using OrdinaryDiffEqDefault

    fm1 = OneHourFuelMoisture()
    r = RothermelFireSpread()
    cs = couple(fm1, r)
    sys = convert(System, cs; compile = false)
    compiled = mtkcompile(sys)

    # Wet conditions: EMCPRM=0.15 → MC1=0.1545, exceeds Mx=0.12
    # At moisture extinction: rM=1.0, η_M=0, IR=0, R=0
    prob = ODEProblem(
        compiled,
        Dict(
            compiled.OneHourFuelMoisture₊EMCPRM => 0.15,
            compiled.OneHourFuelMoisture₊use_fuel_sticks => 0.0,
            compiled.OneHourFuelMoisture₊is_raining => 0.0,
            compiled.OneHourFuelMoisture₊MC10 => 0.0,
            compiled.RothermelFireSpread₊σ => 11483.5,
            compiled.RothermelFireSpread₊w0 => 0.166,
            compiled.RothermelFireSpread₊δ => 0.3048,
            compiled.RothermelFireSpread₊Mx => 0.12,
            compiled.RothermelFireSpread₊h => 18608000.0,
            compiled.RothermelFireSpread₊U => 0.0,
            compiled.RothermelFireSpread₊tanϕ => 0.0,
        ),
        (0.0, 1.0),
    )
    sol = solve(prob)

    mc1_val = sol[compiled.OneHourFuelMoisture₊MC1][end]
    R_val = sol[compiled.RothermelFireSpread₊R][end]
    IR_val = sol[compiled.RothermelFireSpread₊IR][end]

    # MC1 = 1.03 * 0.15 = 0.1545 exceeds Mx = 0.12
    @test mc1_val > 0.12
    # Fire should NOT spread at moisture extinction
    @test R_val ≈ 0.0 atol = 1.0e-8
    @test IR_val ≈ 0.0 atol = 1.0e-8
end

@testitem "Full EMC chain -- moisture extinction from high humidity" setup = [CouplingSetup] tags = [:coupling] begin
    using ModelingToolkit: mtkcompile
    using OrdinaryDiffEqDefault

    emc = EquilibriumMoistureContent()
    fm1 = OneHourFuelMoisture()
    r = RothermelFireSpread()
    cs = couple(emc, fm1, r)
    sys = convert(System, cs; compile = false)
    compiled = mtkcompile(sys)

    # Fuel Model 1 (short grass) base parameters
    base_params = Dict(
        compiled.OneHourFuelMoisture₊use_fuel_sticks => 0.0,
        compiled.OneHourFuelMoisture₊is_raining => 0.0,
        compiled.OneHourFuelMoisture₊MC10 => 0.0,
        compiled.RothermelFireSpread₊σ => 11483.5,
        compiled.RothermelFireSpread₊w0 => 0.166,
        compiled.RothermelFireSpread₊δ => 0.3048,
        compiled.RothermelFireSpread₊Mx => 0.12,
        compiled.RothermelFireSpread₊h => 18608000.0,
        compiled.RothermelFireSpread₊U => 2.235,
        compiled.RothermelFireSpread₊tanϕ => 0.0,
    )

    # Dry case: TEMP=313.15K (104°F), RH=5% → EMC very low → R > 0
    dry_params = merge(
        base_params, Dict(
            compiled.EquilibriumMoistureContent₊TEMP => 313.15,
            compiled.EquilibriumMoistureContent₊RH => 0.05,
        )
    )
    prob_dry = ODEProblem(compiled, dry_params, (0.0, 1.0))
    sol_dry = solve(prob_dry)
    R_dry = sol_dry[compiled.RothermelFireSpread₊R][end]
    @test R_dry > 0

    # Wet case: TEMP=280.15K (45°F), RH=90% → EMC >> Mx → R ≈ 0
    wet_params = merge(
        base_params, Dict(
            compiled.EquilibriumMoistureContent₊TEMP => 280.15,
            compiled.EquilibriumMoistureContent₊RH => 0.9,
        )
    )
    prob_wet = ODEProblem(compiled, wet_params, (0.0, 1.0))
    sol_wet = solve(prob_wet)
    R_wet = sol_wet[compiled.RothermelFireSpread₊R][end]
    @test R_wet ≈ 0.0 atol = 1.0e-10

    # Sweep: R should decrease monotonically as RH increases
    RH_vals = 0.05:0.1:0.95
    R_vals = Float64[]
    for rh in RH_vals
        params = merge(
            base_params, Dict(
                compiled.EquilibriumMoistureContent₊TEMP => 294.26,
                compiled.EquilibriumMoistureContent₊RH => rh,
            )
        )
        prob = ODEProblem(compiled, params, (0.0, 1.0))
        sol = solve(prob)
        push!(R_vals, sol[compiled.RothermelFireSpread₊R][end])
    end
    # Monotonically non-increasing
    for i in 2:length(R_vals)
        @test R_vals[i] <= R_vals[i - 1] + 1.0e-12
    end
    # At high RH, fire should be extinguished
    @test R_vals[end] ≈ 0.0 atol = 1.0e-10
end
