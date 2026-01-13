"""
Example: Heterogeneous fuel model with multiple size classes
"""

using WildlandFire
using OrdinaryDiffEqDefault

println("\n" * "="^80)
println("ROTHERMEL HETEROGENEOUS FUEL MODEL")
println("="^80 * "\n")

println("Creating heterogeneous fuel model with 1 dead and 1 live class...")

sys = create_heterogeneous_rothermel_model(1, 1)

# Parameters for a grass/shrub mixture
params = [
    # Dead fuel (1-hr timelag)
    h_d[1] => 8000.0,
    S_T_d[1] => 0.0555,
    S_e_d[1] => 0.010,
    ρ_p_d[1] => 32.0,
    σ_d[1] => 3500.0,
    w_o_d[1] => 0.10,
    M_f_d[1] => 0.06,

    # Live fuel (herbaceous)
    h_l[1] => 8000.0,
    S_T_l[1] => 0.0555,
    S_e_l[1] => 0.010,
    ρ_p_l[1] => 32.0,
    σ_l[1] => 1800.0,
    w_o_l[1] => 0.15,
    M_f_l[1] => 1.00,  # 100% moisture for live fuel

    # Common parameters
    δ => 1.5,
    U => 352.0,  # 4 mph
    tan_ϕ => 0.0,
]

# Solve
prob = ODEProblem(sys, [], (0.0, 1.0), params)
sol = solve(prob)

R_val = sol[R][end]
println("Rate of spread: $(round(R_val, digits=2)) ft/min")
println("Rate of spread: $(round(R_val * 1.829, digits=2)) chains/hr")

println("\n" * "="^80)
