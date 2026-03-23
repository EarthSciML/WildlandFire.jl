@testsnippet LevelSetSetup begin
    using Test
    using ModelingToolkit
    using ModelingToolkit: t, D
    using DynamicQuantities
    using WildlandFire
    using EarthSciMLBase
    using MethodOfLines
    using OrdinaryDiffEqDefault
    using DomainSets
end

@testitem "LevelSetFireSpread - Structural Verification" setup = [LevelSetSetup] tags = [:levelset] begin
    @parameters x [unit = u"m"]
    @parameters y [unit = u"m"]
    domain = DomainInfo(
        constIC(0.0, t ∈ Interval(0.0, 10.0)),
        constBC(0.0, x ∈ Interval(0.0, 100.0), y ∈ Interval(0.0, 100.0)),
    )
    sys = LevelSetFireSpread(
        domain;
        initial_condition = (x, y) -> sqrt((x - 50.0)^2 + (y - 50.0)^2) - 10.0,
        spread_rate = 1.0,
    )

    # PDESystem structure
    @test sys isa PDESystem
    @test length(equations(sys)) == 1   # One PDE: level-set equation
    @test length(sys.bcs) == 5          # 1 IC + 4 BCs (Neumann on all sides)
    @test length(sys.ivs) == 3          # t, x, y
    @test length(sys.dvs) == 1          # ψ(t, x, y)

    # Parameters include S and psi_ref
    @test length(sys.ps) == 2
end

@testitem "LevelSetFireSpread - Circular Spread (Isotropic)" setup = [LevelSetSetup] tags = [:levelset] begin
    # Analytical test: with constant S and no wind/slope, fire spreads as a circle
    # Initial: circle of radius r₀ = 10m centered at (50, 50)
    # At time t, the radius should be approximately r₀ + S*t
    # Eq. 9, Mandel et al. (2011)

    r0 = 10.0  # initial radius (m)
    S_val = 1.0  # spread rate (m/s)
    domain_size = 100.0
    center = domain_size / 2.0
    t_end = 5.0

    @parameters x [unit = u"m"]
    @parameters y [unit = u"m"]
    domain = DomainInfo(
        constIC(0.0, t ∈ Interval(0.0, t_end)),
        constBC(0.0, x ∈ Interval(0.0, domain_size), y ∈ Interval(0.0, domain_size)),
    )
    sys = LevelSetFireSpread(
        domain;
        initial_condition = (x, y) -> sqrt((x - center)^2 + (y - center)^2) - r0,
        spread_rate = S_val,
    )

    dx = 5.0
    discretization = MOLFiniteDifference(
        [sys.ivs[2] => dx, sys.ivs[3] => dx], sys.ivs[1];
        advection_scheme = WENOScheme()
    )
    prob = MethodOfLines.discretize(sys, discretization; checks = false)
    sol = solve(prob; saveat = t_end)

    @test sol.retcode == SciMLBase.ReturnCode.Success

    # Get solution
    psi = sol[sys.dvs[1]]

    # At t=0, fire area should correspond to a circle of radius r₀
    # Burning cells: ψ ≤ 0
    burning_t0 = count(x -> x <= 0, psi[1, :, :])
    burning_tend = count(x -> x <= 0, psi[end, :, :])

    # Fire should have spread: more burning cells at t_end
    @test burning_tend > burning_t0

    # Check that center value stays negative (well inside burning region)
    mid = div(size(psi, 2), 2) + 1
    @test psi[end, mid, mid] < 0

    # Estimate fire radius from burning area:
    # Area ≈ n_burning * dx^2, radius ≈ sqrt(Area / π)
    area_t0 = burning_t0 * dx^2
    area_tend = burning_tend * dx^2
    r_t0 = sqrt(area_t0 / π)
    r_tend = sqrt(area_tend / π)

    # Expected radius at t_end: r₀ + S*t_end = 10 + 1*5 = 15m
    r_expected = r0 + S_val * t_end

    # Tolerance allows for discretization error on a coarse grid (dx=5m)
    # WENO reduces errors vs standard FD but coarser grid increases spatial error
    @test abs(r_tend - r_expected) / r_expected < 0.35
end

@testitem "LevelSetFireSpread - Custom Boundary Conditions" setup = [LevelSetSetup] tags = [:levelset] begin
    # Test that custom boundary conditions can be provided
    @parameters x [description = "x", unit = u"m"]
    @parameters y [description = "y", unit = u"m"]
    @variables ψ(..) [description = "ψ", unit = u"m"]

    Dx = Differential(x)
    Dy = Differential(y)

    # Dirichlet boundaries: ψ = 100m (far from fire)
    @constants ψ_bc = 100.0 [description = "Boundary level-set value", unit = u"m"]
    custom_bcs = [
        ψ(t, 0.0, y) ~ ψ_bc,
        ψ(t, 100.0, y) ~ ψ_bc,
        ψ(t, x, 0.0) ~ ψ_bc,
        ψ(t, x, 100.0) ~ ψ_bc,
    ]

    domain = DomainInfo(
        constIC(0.0, t ∈ Interval(0.0, 5.0)),
        constBC(0.0, x ∈ Interval(0.0, 100.0), y ∈ Interval(0.0, 100.0)),
    )
    sys = LevelSetFireSpread(
        domain;
        initial_condition = (x, y) -> sqrt((x - 50.0)^2 + (y - 50.0)^2) - 10.0,
        boundary_conditions = custom_bcs,
    )

    @test length(sys.bcs) == 5  # 1 IC + 4 custom BCs
end

@testitem "FuelConsumption - Structural Verification" setup = [LevelSetSetup] tags = [:levelset] begin
    fc = FuelConsumption()
    @test fc !== nothing
    @test length(equations(fc)) == 1   # dF/dt = -is_burning * F / T_f
    @test length(unknowns(fc)) == 2    # F, is_burning
end

@testitem "FuelConsumption - Exponential Decay" setup = [LevelSetSetup] tags = [:levelset] begin
    # Eq. 3, Mandel et al. (2011): F(t) = exp(-t/T_f)
    # Test that fuel fraction decays exponentially when burning
    fc = FuelConsumption()
    fc_nns = ModelingToolkit.toggle_namespacing(fc, false)
    compiled_sys = mtkcompile(fc; inputs = [fc_nns.is_burning])

    T_f_val = 10.0  # 10 second burn time
    prob = ODEProblem(
        compiled_sys,
        merge(Dict(compiled_sys.F => 1.0), Dict(compiled_sys.T_f => T_f_val, compiled_sys.is_burning => 1.0)),
        (0.0, 30.0),
    )
    sol = solve(prob)

    @test sol.retcode == SciMLBase.ReturnCode.Success

    # Check exponential decay at several time points
    for t_check in [T_f_val, 2 * T_f_val, 3 * T_f_val]
        F_expected = exp(-t_check / T_f_val)
        F_computed = sol(t_check)[1]
        @test isapprox(F_computed, F_expected, rtol = 1.0e-4)
    end

    # At t = T_f, F should be approximately 1/e ≈ 0.3679
    @test isapprox(sol(T_f_val)[1], 1 / exp(1), rtol = 1.0e-4)
end

@testitem "FuelConsumption - Not Burning" setup = [LevelSetSetup] tags = [:levelset] begin
    # When is_burning = 0, fuel should not decay
    fc = FuelConsumption()
    fc_nns = ModelingToolkit.toggle_namespacing(fc, false)
    compiled_sys = mtkcompile(fc; inputs = [fc_nns.is_burning])

    prob = ODEProblem(
        compiled_sys,
        merge(Dict(compiled_sys.F => 1.0), Dict(compiled_sys.T_f => 10.0, compiled_sys.is_burning => 0.0)),
        (0.0, 100.0),
    )
    sol = solve(prob)

    @test sol.retcode == SciMLBase.ReturnCode.Success
    @test isapprox(sol(100.0)[1], 1.0, atol = 1.0e-10)
end

@testitem "FireHeatFlux - Structural Verification" setup = [LevelSetSetup] tags = [:levelset] begin
    fhf = FireHeatFlux()
    @test fhf !== nothing
    @test length(equations(fhf)) == 2    # phi_h and phi_q
    @test length(unknowns(fhf)) == 3     # fuel_burn_rate, phi_h, phi_q
end

@testitem "FireHeatFlux - Equation Verification" setup = [LevelSetSetup] tags = [:levelset] begin
    using NonlinearSolve

    # Verify heat flux values against known inputs
    # Eq. 4: phi_h = burn_rate * w_l * h / (1 + M_f)
    # Eq. 5: phi_q = burn_rate * (M_f + 0.56) / (1 + M_f) * L * w_l
    fhf = FireHeatFlux()
    fhf_nns = ModelingToolkit.toggle_namespacing(fhf, false)
    compiled_sys = mtkcompile(fhf; inputs = [fhf_nns.fuel_burn_rate])

    w_l_val = 0.166     # kg/m^2 (fuel model 1)
    h_val = 8000.0 * 1055.06 / 0.453592  # J/kg (8000 BTU/lb in SI)
    M_f_val = 0.08      # 8% moisture
    burn_rate_val = 0.1  # 1/s

    prob = NonlinearProblem(
        compiled_sys,
        merge(
            Dict(compiled_sys.phi_h => 0.0, compiled_sys.phi_q => 0.0),
            Dict(
                compiled_sys.fuel_burn_rate => burn_rate_val,
                compiled_sys.w_l => w_l_val,
                compiled_sys.h_fuel => h_val,
                compiled_sys.M_f => M_f_val,
            ),
        ),
    )
    sol = solve(prob)

    # Expected sensible heat flux — Eq. 4, Mandel et al. (2011)
    phi_h_expected = burn_rate_val * w_l_val * h_val / (1.0 + M_f_val)
    @test isapprox(sol[compiled_sys.phi_h], phi_h_expected, rtol = 1.0e-6)

    # Expected latent heat flux — Eq. 5, Mandel et al. (2011)
    L_water = 2.5e6
    phi_q_expected = burn_rate_val * (M_f_val + 0.56) / (1.0 + M_f_val) * L_water * w_l_val
    @test isapprox(sol[compiled_sys.phi_q], phi_q_expected, rtol = 1.0e-6)
end

@testitem "anderson_fuel_coefficients - All Models" setup = [LevelSetSetup] tags = [:levelset] begin
    # Test that all 13 Anderson fuel models can be computed
    for model in 1:13
        coeffs = anderson_fuel_coefficients(model)

        # All spread rates should be non-negative
        @test coeffs.S_0 >= 0.0
        @test coeffs.R_0 >= 0.0
        @test coeffs.c >= 0.0
        @test coeffs.d >= 0.0
        @test coeffs.e > 0.0

        # Wind exponent b should be positive
        @test coeffs.b > 0.0

        # Fuel burn time should be positive
        @test coeffs.T_f > 0.0

        # Fuel load should be positive
        @test coeffs.w_l > 0.0
    end

    # Invalid model should error
    @test_throws ErrorException anderson_fuel_coefficients(0)
    @test_throws ErrorException anderson_fuel_coefficients(14)
end

@testitem "anderson_fuel_coefficients - Fuel Model 1 Values" setup = [LevelSetSetup] tags = [:levelset] begin
    # Fuel Model 1 (short grass) sanity checks
    coeffs = anderson_fuel_coefficients(1)

    # R_0 should be a reasonable spread rate for short grass
    # Typical R_0 for FM1 at 8% moisture: ~1-3 m/min = 0.017-0.05 m/s
    @test 0.005 < coeffs.R_0 < 0.1

    # Fuel load: 0.166 kg/m^2 (from Anderson)
    @test coeffs.w_l == 0.166

    # T_f for FM1: weight=7, T_f = 7/0.8514 ≈ 8.22 s
    @test isapprox(coeffs.T_f, 7.0 / 0.8514, rtol = 1.0e-4)

    # Heat content: 8000 BTU/lb converted to SI (J/kg)
    h_expected = 8000.0 * 1055.06 / 0.453592  # ≈ 18,608,000 J/kg
    @test isapprox(coeffs.h, h_expected, rtol = 1.0e-4)
end

@testitem "anderson_fuel_coefficients - Moisture Sensitivity" setup = [LevelSetSetup] tags = [:levelset] begin
    # Higher moisture should reduce spread rate — Table 2, Mandel et al. (2011)
    coeffs_dry = anderson_fuel_coefficients(1; M_f = 0.04)
    coeffs_wet = anderson_fuel_coefficients(1; M_f = 0.1)

    @test coeffs_dry.R_0 > coeffs_wet.R_0
end

@testitem "ANDERSON_FUEL_DATA - Table 1 Properties" setup = [LevelSetSetup] tags = [:levelset] begin
    # Verify the raw fuel data matches Table 1 of Mandel et al. (2011)
    @test length(ANDERSON_FUEL_DATA) == 13

    # All models should have positive fuel properties
    for i in 1:13
        fd = ANDERSON_FUEL_DATA[i]
        @test fd.fgi > 0.0       # Total fuel load (kg/m²)
        @test fd.depth > 0.0     # Fuel bed depth (m)
        @test fd.savr > 0        # Surface-area-to-volume ratio (1/ft)
        @test 0.0 < fd.mce < 1.0 # Moisture of extinction (fraction)
        @test fd.dens > 0.0      # Fuel particle density (lb/ft³)
        @test fd.st > 0.0        # Total mineral content
        @test fd.se > 0.0        # Effective mineral content
        @test fd.weight > 0      # Fuel weight parameter (s)
        @test fd.h > 0.0         # Heat content (BTU/lb)
    end

    # Spot-check Fuel Model 1 (short grass) values
    fm1 = ANDERSON_FUEL_DATA[1]
    @test fm1.fgi ≈ 0.166       # 0.74 tons/acre ≈ 0.166 kg/m²
    @test fm1.depth ≈ 0.305     # 1.0 ft ≈ 0.305 m
    @test fm1.savr == 3500      # 1/ft
    @test fm1.mce ≈ 0.12
    @test fm1.dens ≈ 32.0       # lb/ft³
end
