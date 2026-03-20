@testsnippet FSimSetup begin
    using Test
    using ModelingToolkit
    using NonlinearSolve
    using Symbolics
    using WildlandFire
end

# ============================================================================
# FireOccurrenceLogistic Tests
# ============================================================================

@testitem "FireOccurrenceLogistic - Structural Verification" setup = [FSimSetup] tags = [:fsim] begin
    sys = FireOccurrenceLogistic()

    # Verify system structure
    @test sys !== nothing
    @test length(equations(sys)) == 1
    @test length(unknowns(sys)) == 1

    # Check variable names
    var_names = [string(Symbolics.tosymbol(v, escape = false)) for v in unknowns(sys)]
    @test "P_fire" in var_names

    # Check parameter names
    param_names = [string(Symbolics.tosymbol(p, escape = false)) for p in parameters(sys)]
    @test "β₀" in param_names
    @test "β₁" in param_names
    @test "ERC" in param_names
end

@testitem "FireOccurrenceLogistic - Equation Verification" setup = [FSimSetup] tags = [:fsim] begin
    sys = FireOccurrenceLogistic()
    cs = mtkcompile(sys)

    # Test 1: Known logistic values
    # At ERC=0 with β₀=-5, β₁=0.05: P = 1/(1+exp(5)) ≈ 0.00669
    prob = NonlinearProblem(cs, Dict(cs.β₀ => -5.0, cs.β₁ => 0.05, cs.ERC => 0.0))
    sol = solve(prob)
    @test sol[cs.P_fire] ≈ 1 / (1 + exp(5.0)) rtol = 1.0e-10

    # Test 2: At ERC=100 with β₀=-5, β₁=0.05: P = 1/(1+exp(0)) = 0.5
    prob2 = NonlinearProblem(cs, Dict(cs.β₀ => -5.0, cs.β₁ => 0.05, cs.ERC => 100.0))
    sol2 = solve(prob2)
    @test sol2[cs.P_fire] ≈ 0.5 rtol = 1.0e-10

    # Test 3: At ERC=60 with β₀=-5, β₁=0.05: P = 1/(1+exp(2)) ≈ 0.1192
    prob3 = NonlinearProblem(cs, Dict(cs.β₀ => -5.0, cs.β₁ => 0.05, cs.ERC => 60.0))
    sol3 = solve(prob3)
    @test sol3[cs.P_fire] ≈ 1 / (1 + exp(2.0)) rtol = 1.0e-10
end

@testitem "FireOccurrenceLogistic - Qualitative Properties" setup = [FSimSetup] tags = [:fsim] begin
    sys = FireOccurrenceLogistic()
    cs = mtkcompile(sys)

    β₀ = -5.0
    β₁ = 0.05

    # Test monotonicity: P(fire) increases with ERC
    ERC_values = [0.0, 20.0, 40.0, 60.0, 80.0, 100.0]
    P_values = Float64[]
    for erc in ERC_values
        prob = NonlinearProblem(cs, Dict(cs.β₀ => β₀, cs.β₁ => β₁, cs.ERC => erc))
        sol = solve(prob)
        push!(P_values, sol[cs.P_fire])
    end
    for i in 2:length(P_values)
        @test P_values[i] > P_values[i - 1]
    end

    # Test boundedness: 0 < P < 1 for all ERC
    for p in P_values
        @test 0.0 < p < 1.0
    end

    # Test: higher ERC thresholds (larger β₀ magnitude) give lower probabilities
    # Consistent with Fig. 4a: larger fires are less likely
    prob_small = NonlinearProblem(cs, Dict(cs.β₀ => -3.0, cs.β₁ => β₁, cs.ERC => 60.0))
    prob_large = NonlinearProblem(cs, Dict(cs.β₀ => -7.0, cs.β₁ => β₁, cs.ERC => 60.0))
    sol_small = solve(prob_small)
    sol_large = solve(prob_large)
    @test sol_small[cs.P_fire] > sol_large[cs.P_fire]
end

@testitem "FireOccurrenceLogistic - Limiting Behavior" setup = [FSimSetup] tags = [:fsim] begin
    sys = FireOccurrenceLogistic()
    cs = mtkcompile(sys)

    # At very high ERC, probability should approach 1
    prob_high = NonlinearProblem(cs, Dict(cs.β₀ => -5.0, cs.β₁ => 0.05, cs.ERC => 500.0))
    sol_high = solve(prob_high)
    @test sol_high[cs.P_fire] > 0.999

    # At very low ERC, probability should approach 0
    prob_low = NonlinearProblem(cs, Dict(cs.β₀ => -5.0, cs.β₁ => 0.05, cs.ERC => -200.0))
    sol_low = solve(prob_low)
    @test sol_low[cs.P_fire] < 0.001

    # At β₀ + β₁*ERC = 0 (inflection point), probability = 0.5
    prob_half = NonlinearProblem(cs, Dict(cs.β₀ => -5.0, cs.β₁ => 0.05, cs.ERC => 100.0))
    sol_half = solve(prob_half)
    @test sol_half[cs.P_fire] ≈ 0.5 rtol = 1.0e-10
end

# ============================================================================
# FireContainment Tests
# ============================================================================

@testitem "FireContainment - Structural Verification" setup = [FSimSetup] tags = [:fsim] begin
    sys = FireContainment()

    @test sys !== nothing
    @test length(equations(sys)) == 1
    @test length(unknowns(sys)) == 1

    var_names = [string(Symbolics.tosymbol(v, escape = false)) for v in unknowns(sys)]
    @test "P_contain" in var_names

    param_names = [string(Symbolics.tosymbol(p, escape = false)) for p in parameters(sys)]
    @test "α₀" in param_names
    @test "α₁" in param_names
    @test "α₂" in param_names
    @test "α₃" in param_names
    @test "NPI" in param_names
    @test "is_high_spread" in param_names
    @test "is_timber" in param_names
end

@testitem "FireContainment - Equation Verification" setup = [FSimSetup] tags = [:fsim] begin
    sys = FireContainment()
    cs = mtkcompile(sys)

    # Test with known parameters
    α₀, α₁, α₂, α₃ = 0.5, 0.3, -1.5, -0.8

    # Low spread, no timber, NPI=5
    prob = NonlinearProblem(
        cs, Dict(
            cs.α₀ => α₀, cs.α₁ => α₁, cs.α₂ => α₂, cs.α₃ => α₃,
            cs.NPI => 5.0, cs.is_high_spread => 0.0, cs.is_timber => 0.0
        )
    )
    sol = solve(prob)
    expected = 1 / (1 + exp(-(α₀ + α₁ * 5.0)))
    @test sol[cs.P_contain] ≈ expected rtol = 1.0e-10

    # High spread, timber, NPI=5
    prob2 = NonlinearProblem(
        cs, Dict(
            cs.α₀ => α₀, cs.α₁ => α₁, cs.α₂ => α₂, cs.α₃ => α₃,
            cs.NPI => 5.0, cs.is_high_spread => 1.0, cs.is_timber => 1.0
        )
    )
    sol2 = solve(prob2)
    expected2 = 1 / (1 + exp(-(α₀ + α₁ * 5.0 + α₂ * 1.0 + α₃ * 1.0)))
    @test sol2[cs.P_contain] ≈ expected2 rtol = 1.0e-10
end

@testitem "FireContainment - Qualitative Properties" setup = [FSimSetup] tags = [:fsim] begin
    sys = FireContainment()
    cs = mtkcompile(sys)

    α₀, α₁, α₂, α₃ = 0.5, 0.3, -1.5, -0.8

    # Test: containment probability increases with NPI (fire duration)
    # Consistent with Fig. 6: longer fires are more likely to be contained
    P_values = Float64[]
    for npi in [1.0, 5.0, 10.0, 15.0, 20.0]
        prob = NonlinearProblem(
            cs, Dict(
                cs.α₀ => α₀, cs.α₁ => α₁, cs.α₂ => α₂, cs.α₃ => α₃,
                cs.NPI => npi, cs.is_high_spread => 0.0, cs.is_timber => 0.0
            )
        )
        sol = solve(prob)
        push!(P_values, sol[cs.P_contain])
    end
    for i in 2:length(P_values)
        @test P_values[i] > P_values[i - 1]
    end

    # Test: containment lower during high spread intervals
    # Consistent with Fig. 6
    prob_low = NonlinearProblem(
        cs, Dict(
            cs.α₀ => α₀, cs.α₁ => α₁, cs.α₂ => α₂, cs.α₃ => α₃,
            cs.NPI => 5.0, cs.is_high_spread => 0.0, cs.is_timber => 0.0
        )
    )
    prob_high = NonlinearProblem(
        cs, Dict(
            cs.α₀ => α₀, cs.α₁ => α₁, cs.α₂ => α₂, cs.α₃ => α₃,
            cs.NPI => 5.0, cs.is_high_spread => 1.0, cs.is_timber => 0.0
        )
    )
    sol_low = solve(prob_low)
    sol_high = solve(prob_high)
    @test sol_low[cs.P_contain] > sol_high[cs.P_contain]

    # Test: containment lower in timber fuels
    # Consistent with Fig. 6
    prob_no_timber = NonlinearProblem(
        cs, Dict(
            cs.α₀ => α₀, cs.α₁ => α₁, cs.α₂ => α₂, cs.α₃ => α₃,
            cs.NPI => 5.0, cs.is_high_spread => 0.0, cs.is_timber => 0.0
        )
    )
    prob_timber = NonlinearProblem(
        cs, Dict(
            cs.α₀ => α₀, cs.α₁ => α₁, cs.α₂ => α₂, cs.α₃ => α₃,
            cs.NPI => 5.0, cs.is_high_spread => 0.0, cs.is_timber => 1.0
        )
    )
    sol_no_timber = solve(prob_no_timber)
    sol_timber = solve(prob_timber)
    @test sol_no_timber[cs.P_contain] > sol_timber[cs.P_contain]

    # Test: all probabilities bounded in [0, 1]
    for p in P_values
        @test 0.0 < p < 1.0
    end
end

# ============================================================================
# BurnProbability Tests
# ============================================================================

@testitem "BurnProbability - Structural Verification" setup = [FSimSetup] tags = [:fsim] begin
    sys = BurnProbability()

    @test sys !== nothing
    @test length(equations(sys)) == 3
    @test length(unknowns(sys)) == 3

    var_names = [string(Symbolics.tosymbol(v, escape = false)) for v in unknowns(sys)]
    @test "BP" in var_names
    @test "BP_FPU" in var_names
    @test "BP_conditional" in var_names
end

@testitem "BurnProbability - Equation Verification" setup = [FSimSetup] tags = [:fsim] begin
    sys = BurnProbability()
    cs = mtkcompile(sys)

    # Test cell-level burn probability
    prob = NonlinearProblem(
        cs, Dict(
            cs.N_burned => 50.0,
            cs.N_years => 10000.0,
            cs.A_burned => 1.0e9,
            cs.A_FPU => 1.0e10,
            cs.P_FL_given_burn => 0.3
        )
    )
    sol = solve(prob)

    @test sol[cs.BP] ≈ 50.0 / 10000.0 rtol = 1.0e-10
    @test sol[cs.BP_FPU] ≈ 1.0e9 / (1.0e10 * 10000.0) rtol = 1.0e-10
    @test sol[cs.BP_conditional] ≈ (50.0 / 10000.0) * 0.3 rtol = 1.0e-10

    # Test with typical historical values from Table 1
    # Northwest California: BP ≈ 0.00127, FPU area ≈ 3,044,575 ha
    prob2 = NonlinearProblem(
        cs, Dict(
            cs.N_burned => 127.0,
            cs.N_years => 100000.0,
            cs.A_burned => 3.866e9,
            cs.A_FPU => 3.044575e10,
            cs.P_FL_given_burn => 1.0
        )
    )
    sol2 = solve(prob2)
    @test sol2[cs.BP] ≈ 0.00127 rtol = 1.0e-10
    @test sol2[cs.BP_FPU] ≈ 3.866e9 / (3.044575e10 * 100000.0) rtol = 1.0e-6
end

@testitem "BurnProbability - Conservation" setup = [FSimSetup] tags = [:fsim] begin
    sys = BurnProbability()
    cs = mtkcompile(sys)

    # Test: conditional burn probability with P(FL|burn)=1 equals total burn probability
    prob = NonlinearProblem(
        cs, Dict(
            cs.N_burned => 100.0,
            cs.N_years => 10000.0,
            cs.A_burned => 1.0e9,
            cs.A_FPU => 1.0e10,
            cs.P_FL_given_burn => 1.0
        )
    )
    sol = solve(prob)
    @test sol[cs.BP_conditional] ≈ sol[cs.BP] rtol = 1.0e-10

    # Test: conditional burn probability with P(FL|burn)=0 equals 0
    prob2 = NonlinearProblem(
        cs, Dict(
            cs.N_burned => 100.0,
            cs.N_years => 10000.0,
            cs.A_burned => 1.0e9,
            cs.A_FPU => 1.0e10,
            cs.P_FL_given_burn => 0.0
        )
    )
    sol2 = solve(prob2)
    @test sol2[cs.BP_conditional] ≈ 0.0 atol = 1.0e-15
end

# ============================================================================
# ERCTimeSeries Tests
# ============================================================================

@testitem "ERCTimeSeries - Structural Verification" setup = [FSimSetup] tags = [:fsim] begin
    sys = ERCTimeSeries()

    @test sys !== nothing
    @test length(equations(sys)) == 3
    @test length(unknowns(sys)) == 3

    var_names = [string(Symbolics.tosymbol(v, escape = false)) for v in unknowns(sys)]
    @test "ERC_hat" in var_names
    @test "var_a" in var_names
    @test "σ²_daily" in var_names

    param_names = [string(Symbolics.tosymbol(p, escape = false)) for p in parameters(sys)]
    @test "f_t" in param_names
    @test "a_prev" in param_names
    @test "a_current" in param_names
    @test "φ₁" in param_names
    @test "ρ₁" in param_names
    @test "s²" in param_names
end

@testitem "ERCTimeSeries - Equation Verification" setup = [FSimSetup] tags = [:fsim] begin
    sys = ERCTimeSeries()
    cs = mtkcompile(sys)

    # Test Eq. 2: ẑ(t) = f(t) + φ₁·a(t-1) + a(t)
    f_t, a_prev, a_curr, φ₁, ρ₁, s2 = 50.0, 3.0, 1.5, 0.7, 0.75, 100.0
    prob = NonlinearProblem(
        cs, Dict(
            cs.f_t => f_t, cs.a_prev => a_prev, cs.a_current => a_curr,
            cs.φ₁ => φ₁, cs.ρ₁ => ρ₁, cs.s² => s2
        )
    )
    sol = solve(prob)

    @test sol[cs.ERC_hat] ≈ f_t + φ₁ * a_prev + a_curr rtol = 1.0e-10

    # Test Eq. 4: var(a(t)) = s²(t) × (1 - ρ₁φ₁)
    @test sol[cs.var_a] ≈ s2 * (1 - ρ₁ * φ₁) rtol = 1.0e-10

    # Test Eq. 3: σ²(t) = var(a(t)) / (1 - ρ₁φ₁)
    @test sol[cs.σ²_daily] ≈ sol[cs.var_a] / (1 - ρ₁ * φ₁) rtol = 1.0e-10

    # Eqs. 3 and 4 combined: σ²_daily should equal s²
    @test sol[cs.σ²_daily] ≈ s2 rtol = 1.0e-10
end

@testitem "ERCTimeSeries - Qualitative Properties" setup = [FSimSetup] tags = [:fsim] begin
    sys = ERCTimeSeries()
    cs = mtkcompile(sys)

    # Test: ERC follows seasonal trend when noise is zero
    prob = NonlinearProblem(
        cs, Dict(
            cs.f_t => 60.0, cs.a_prev => 0.0, cs.a_current => 0.0,
            cs.φ₁ => 0.7, cs.ρ₁ => 0.75, cs.s² => 100.0
        )
    )
    sol = solve(prob)
    @test sol[cs.ERC_hat] ≈ 60.0 rtol = 1.0e-10

    # Test: variance is positive for valid parameters
    @test sol[cs.var_a] > 0.0
    @test sol[cs.σ²_daily] > 0.0

    # Test: higher autocorrelation (closer to 1) reduces innovation variance
    for ρ in [0.3, 0.5, 0.7, 0.9]
        φ = ρ  # For AR(1), optimal φ₁ = ρ₁
        prob = NonlinearProblem(
            cs, Dict(
                cs.f_t => 50.0, cs.a_prev => 0.0, cs.a_current => 0.0,
                cs.φ₁ => φ, cs.ρ₁ => ρ, cs.s² => 100.0
            )
        )
        sol = solve(prob)
        # Innovation variance should decrease as autocorrelation increases
        @test sol[cs.var_a] ≈ 100.0 * (1 - ρ^2) rtol = 1.0e-10
    end
end

# ============================================================================
# FlameLengthCategory Tests
# ============================================================================

@testitem "FlameLengthCategory - Structural Verification" setup = [FSimSetup] tags = [:fsim] begin
    sys = FlameLengthCategory()

    @test sys !== nothing
    @test length(equations(sys)) == 6
    @test length(unknowns(sys)) == 6

    var_names = [string(Symbolics.tosymbol(v, escape = false)) for v in unknowns(sys)]
    for i in 1:6
        @test "cat_$i" in var_names
    end
end

@testitem "FlameLengthCategory - Equation Verification" setup = [FSimSetup] tags = [:fsim] begin
    sys = FlameLengthCategory()
    cs = mtkcompile(sys)

    # Test each category boundary — Finney et al. (2011), Fig. 10
    test_cases = [
        (0.3, [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]),  # FL ≤ 0.6 m → cat 1
        (0.6, [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]),  # FL = 0.6 m → cat 1 (boundary)
        (0.9, [0.0, 1.0, 0.0, 0.0, 0.0, 0.0]),  # 0.6 < FL ≤ 1.2 → cat 2
        (1.5, [0.0, 0.0, 1.0, 0.0, 0.0, 0.0]),  # 1.2 < FL ≤ 1.8 → cat 3
        (2.0, [0.0, 0.0, 0.0, 1.0, 0.0, 0.0]),  # 1.8 < FL ≤ 2.4 → cat 4
        (3.0, [0.0, 0.0, 0.0, 0.0, 1.0, 0.0]),  # 2.4 < FL ≤ 3.7 → cat 5
        (5.0, [0.0, 0.0, 0.0, 0.0, 0.0, 1.0]),  # FL > 3.7 → cat 6
    ]

    for (fl, expected) in test_cases
        prob = NonlinearProblem(cs, Dict(cs.F_L => fl))
        sol = solve(prob)
        cats = [
            sol[cs.cat_1], sol[cs.cat_2], sol[cs.cat_3],
            sol[cs.cat_4], sol[cs.cat_5], sol[cs.cat_6],
        ]
        @test cats ≈ expected atol = 1.0e-10
    end
end

@testitem "FlameLengthCategory - Mutual Exclusivity" setup = [FSimSetup] tags = [:fsim] begin
    sys = FlameLengthCategory()
    cs = mtkcompile(sys)

    # Test: exactly one category is 1 for any flame length
    for fl in [0.1, 0.5, 0.7, 1.0, 1.3, 1.6, 2.0, 2.5, 3.5, 4.0, 10.0]
        prob = NonlinearProblem(cs, Dict(cs.F_L => fl))
        sol = solve(prob)
        cats = [
            sol[cs.cat_1], sol[cs.cat_2], sol[cs.cat_3],
            sol[cs.cat_4], sol[cs.cat_5], sol[cs.cat_6],
        ]
        @test sum(cats) ≈ 1.0 atol = 1.0e-10
    end
end
