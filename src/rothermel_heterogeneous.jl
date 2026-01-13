"""
Rothermel Surface Fire Spread Model - Heterogeneous Fuel Implementation
Based on Rothermel (1972) and Albini (1976a)

This version handles multiple fuel size classes (dead and live fuel categories)
as described in Tables 5 and 6a-6c of the reference document.
"""

using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

"""
Create a Rothermel model for heterogeneous fuel with multiple size classes.

Parameters:
- n_dead: number of dead fuel size classes
- n_live: number of live fuel size classes

The model uses indexed parameters for each size class:
- For dead fuel: categories indexed by i = 1 to n_dead
- For live fuel: categories indexed by i = 1 to n_live
"""
function create_heterogeneous_rothermel_model(n_dead::Int, n_live::Int)

    @parameters begin
        # Environmental parameters
        U                    # Midflame wind speed (ft/min)
        tan_ϕ                # Slope steepness (fraction)

        # Fuel particle parameters for each size class (dead fuel)
        h_d[1:n_dead]        # Heat content (Btu/lb) for each dead fuel class
        S_T_d[1:n_dead]    # Total mineral content (fraction) for each dead fuel class
        S_e_d[1:n_dead]    # Effective mineral content (fraction) for each dead fuel class
        ρ_p_d[1:n_dead]    # Particle density (lb/ft³) for each dead fuel class

        # Fuel particle parameters for each size class (live fuel)
        h_l[1:n_live]        # Heat content (Btu/lb) for each live fuel class
        S_T_l[1:n_live]    # Total mineral content (fraction) for each live fuel class
        S_e_l[1:n_live]    # Effective mineral content (fraction) for each live fuel class
        ρ_p_l[1:n_live]    # Particle density (lb/ft³) for each live fuel class

        # Fuel array parameters for each size class (dead)
        σ_d[1:n_dead]        # Surface-area-to-volume ratio (ft²/ft³) for each dead class
        w_o_d[1:n_dead]    # Oven-dry fuel load (lb/ft²) for each dead class
        M_f_d[1:n_dead]    # Fuel moisture (fraction) for each dead class

        # Fuel array parameters for each size class (live)
        σ_l[1:n_live]        # Surface-area-to-volume ratio (ft²/ft³) for each live class
        w_o_l[1:n_live]    # Oven-dry fuel load (lb/ft²) for each live class
        M_f_l[1:n_live]    # Fuel moisture (fraction) for each live class

        # Common fuel bed parameters
        δ                    # Fuel bed depth (ft)
    end

    @variables begin
        # ========================================================================
        # WEIGHTING FACTORS
        # ========================================================================

        # Surface area per unit fuel cell for each size class
        A_d[1:n_dead](t)     # Dead fuel surface areas
        A_l[1:n_live](t)     # Live fuel surface areas

        # Total surface areas
        A_T_d(t)             # Total surface area of dead categories
        A_T_l(t)             # Total surface area of live categories
        A_T(t)               # Total surface area of all fuel

        # Weighting factors
        f_d[1:n_dead](t)     # Weighting factors for dead fuel classes
        f_l[1:n_live](t)     # Weighting factors for live fuel classes
        f_dead(t)            # Weighting factor for all dead fuel
        f_live(t)            # Weighting factor for all live fuel

        # ========================================================================
        # CHARACTERISTIC VALUES (WEIGHTED AVERAGES)
        # ========================================================================

        # Net fuel load
        w_n_d[1:n_dead](t) # Net fuel load for each dead class
        w_n_l[1:n_live](t) # Net fuel load for each live class
        w_n_d_total(t)     # Total net dead fuel load
        w_n_l_total(t)     # Total net live fuel load
        w_n_total(t)         # Total net fuel load

        # Weighted characteristic values
        h_d_char(t)          # Characteristic heat content for dead fuel
        h_l_char(t)          # Characteristic heat content for live fuel
        S_e_d_char(t)      # Characteristic effective mineral content for dead
        S_e_l_char(t)      # Characteristic effective mineral content for live
        M_f_d_char(t)      # Characteristic moisture for dead fuel
        M_f_l_char(t)      # Characteristic moisture for live fuel

        # ========================================================================
        # MOISTURE OF EXTINCTION
        # ========================================================================

        # Live fuel moisture of extinction
        M_x_l(t)           # Live fuel moisture of extinction

        # Dead-to-live fuel ratio
        W(t)                 # Dead-to-live load ratio

        # Fine dead fuel moisture (weighted average of small size classes)
        M_f_d_fine(t)      # "Fine" dead fuel moisture

        # ========================================================================
        # DAMPING COEFFICIENTS
        # ========================================================================

        # Moisture damping coefficients
        η_M_d(t)           # Moisture damping for dead fuel
        η_M_l(t)           # Moisture damping for live fuel
        η_M_total(t)       # Combined moisture damping

        # Mineral damping coefficients
        η_s_d(t)           # Mineral damping for dead fuel
        η_s_l(t)           # Mineral damping for live fuel
        η_s_total(t)       # Combined mineral damping

        # ========================================================================
        # FUEL BED CHARACTERISTICS
        # ========================================================================

        # Surface-area-to-volume ratios
        σ_d_char(t)          # Characteristic σ for dead fuel
        σ_l_char(t)          # Characteristic σ for live fuel
        σ_total(t)           # Overall characteristic σ

        # Mean bulk density
        ρ_b(t)               # Mean bulk density

        # Mean packing ratio
        β(t)                 # Mean packing ratio

        # Optimum packing ratio
        β_opt(t)             # Optimum packing ratio
        β_ratio(t)           # Relative packing ratio

        # ========================================================================
        # REACTION VELOCITY AND INTENSITY
        # ========================================================================

        Γ_max(t)             # Maximum reaction velocity
        A(t)                 # Constant for reaction velocity
        Γ(t)                 # Optimum reaction velocity
        I_R(t)               # Reaction intensity

        # ========================================================================
        # PROPAGATING FLUX RATIO
        # ========================================================================

        ξ(t)                 # Propagating flux ratio
        I_R_ξ(t)             # No-wind, no-slope propagating flux

        # ========================================================================
        # WIND AND SLOPE FACTORS
        # ========================================================================

        ϕ_w(t)               # Wind factor
        ϕ_s(t)               # Slope factor

        # ========================================================================
        # HEAT SINK
        # ========================================================================

        ε(t)                 # Effective heating number
        Q_ig_d[1:n_dead](t) # Heat of preignition for each dead class
        Q_ig_l[1:n_live](t) # Heat of preignition for each live class
        heat_sink(t)         # Total heat sink

        # ========================================================================
        # RATE OF SPREAD
        # ========================================================================

        R_0(t)               # No-wind, no-slope rate of spread
        R(t)                 # Rate of spread
    end

    # ============================================================================
    # BUILD EQUATIONS
    # ============================================================================

    eqs = Equation[]

    # Surface areas for each fuel class
    for i in 1:n_dead
        push!(eqs, A_d[i] ~ σ_d[i] * w_o_d[i] / ρ_p_d[i])
    end
    for i in 1:n_live
        push!(eqs, A_l[i] ~ σ_l[i] * w_o_l[i] / ρ_p_l[i])
    end

    # Total surface areas
    push!(eqs, A_T_d ~ sum(A_d[i] for i in 1:n_dead))
    push!(eqs, A_T_l ~ sum(A_l[i] for i in 1:n_live))
    push!(eqs, A_T ~ A_T_d + A_T_l)

    # Weighting factors for each size class
    for i in 1:n_dead
        push!(eqs, f_d[i] ~ A_d[i] / A_T)
    end
    for i in 1:n_live
        push!(eqs, f_l[i] ~ A_l[i] / A_T)
    end

    # Weighting factors for dead and live categories
    push!(eqs, f_dead ~ A_T_d / A_T)
    push!(eqs, f_live ~ A_T_l / A_T)

    # Net fuel loads (accounting for mineral content)
    for i in 1:n_dead
        push!(eqs, w_n_d[i] ~ w_o_d[i] * (1 - S_T_d[i]))
    end
    for i in 1:n_live
        push!(eqs, w_n_l[i] ~ w_o_l[i] * (1 - S_T_l[i]))
    end

    # Total net fuel loads
    push!(eqs, w_n_d_total ~ sum(w_n_d[i] for i in 1:n_dead))
    push!(eqs, w_n_l_total ~ sum(w_n_l[i] for i in 1:n_live))
    push!(eqs, w_n_total ~ w_n_d_total + w_n_l_total)

    # Characteristic heat content (weighted by size class weighting factors)
    push!(eqs, h_d_char ~ sum(f_d[i] * h_d[i] for i in 1:n_dead))
    push!(eqs, h_l_char ~ sum(f_l[i] * h_l[i] for i in 1:n_live))

    # Characteristic effective mineral content
    push!(eqs, S_e_d_char ~ sum(f_d[i] * S_e_d[i] for i in 1:n_dead))
    push!(eqs, S_e_l_char ~ sum(f_l[i] * S_e_l[i] for i in 1:n_live))

    # Characteristic moisture content
    push!(eqs, M_f_d_char ~ sum(f_d[i] * M_f_d[i] for i in 1:n_dead))
    push!(eqs, M_f_l_char ~ sum(f_l[i] * M_f_l[i] for i in 1:n_live))

    # Characteristic surface-area-to-volume ratios
    push!(eqs, σ_d_char ~ sum(f_d[i] * σ_d[i] for i in 1:n_dead))
    push!(eqs, σ_l_char ~ sum(f_l[i] * σ_l[i] for i in 1:n_live))
    push!(eqs, σ_total ~ sum(f_d[i] * σ_d[i] for i in 1:n_dead) +
                          sum(f_l[i] * σ_l[i] for i in 1:n_live))

    # Dead-to-live fuel ratio
    push!(eqs, W ~ w_n_d_total / w_n_l_total)

    # Fine dead fuel moisture (assuming first two classes are "fine")
    # This is a simplified version; adjust based on actual size class definitions
    if n_dead >= 2
        push!(eqs, M_f_d_fine ~ (sum(w_o_d[i] * exp(-138/σ_d[i]) for i in 1:min(2,n_dead)) /
                                   sum(w_o_d[i] for i in 1:min(2,n_dead))) *
                                  sum(M_f_d[i] * w_o_d[i] * exp(-138/σ_d[i]) for i in 1:min(2,n_dead)) /
                                  sum(w_o_d[i] * exp(-138/σ_d[i]) for i in 1:min(2,n_dead)))
    else
        push!(eqs, M_f_d_fine ~ M_f_d_char)
    end

    # Live fuel moisture of extinction
    # M_x_l = 2.9W(1 - M_f_dead/(M_x)_d) - 0.226 (min = M_f_l)
    # Assuming (M_x)_d = 0.25 for dead fuel (typical value)
    M_x_d_value = 0.25
    push!(eqs, M_x_l ~ max(2.9 * W * (1 - M_f_d_fine / M_x_d_value) - 0.226, M_f_l_char))

    # Mineral damping coefficients
    push!(eqs, η_s_d ~ min(0.174 * S_e_d_char^(-0.19), 1.0))
    push!(eqs, η_s_l ~ min(0.174 * S_e_l_char^(-0.19), 1.0))

    # Moisture damping coefficients
    # Dead fuel: using standard M_x_d = 0.25
    push!(eqs, η_M_d ~ ifelse(M_f_d_char >= M_x_d_value, 0.0,
                                max(1 - 2.59*(M_f_d_char/M_x_d_value) +
                                    5.11*(M_f_d_char/M_x_d_value)^2 -
                                    3.52*(M_f_d_char/M_x_d_value)^3, 0.0)))

    # Live fuel: using calculated M_x_l
    push!(eqs, η_M_l ~ ifelse(M_f_l_char >= M_x_l, 0.0,
                                max(1 - 2.59*(M_f_l_char/M_x_l) +
                                    5.11*(M_f_l_char/M_x_l)^2 -
                                    3.52*(M_f_l_char/M_x_l)^3, 0.0)))

    # Mean bulk density
    push!(eqs, ρ_b ~ sum(w_o_d[i] for i in 1:n_dead) +
                      sum(w_o_l[i] for i in 1:n_live) / δ)

    # Mean packing ratio
    mean_particle_density = (sum(w_o_d[i] for i in 1:n_dead) +
                            sum(w_o_l[i] for i in 1:n_live)) /
                           (sum(w_o_d[i]/ρ_p_d[i] for i in 1:n_dead) +
                            sum(w_o_l[i]/ρ_p_l[i] for i in 1:n_live))
    push!(eqs, β ~ ρ_b / mean_particle_density)

    # Optimum packing ratio (based on total characteristic σ)
    push!(eqs, β_opt ~ 3.348 * σ_total^(-0.8189))
    push!(eqs, β_ratio ~ β / β_opt)

    # Maximum reaction velocity
    push!(eqs, Γ_max ~ σ_total^1.5 / (495 + 0.0594 * σ_total^1.5))

    # Constant A
    push!(eqs, A ~ 133 * σ_total^(-0.7913))

    # Optimum reaction velocity
    push!(eqs, Γ ~ Γ_max * β_ratio^A * exp(A * (1 - β_ratio)))

    # Combined damping coefficients
    push!(eqs, η_s_total ~ η_s_d * η_s_l)
    push!(eqs, η_M_total ~ η_M_d * η_M_l)

    # Reaction intensity
    push!(eqs, I_R ~ Γ * (w_n_d_total * h_d_char + w_n_l_total * h_l_char) *
                     η_M_total * η_s_total)

    # Propagating flux ratio
    push!(eqs, ξ ~ exp((0.792 + 0.681 * σ_total^0.5) * (β + 0.1)) / (192 + 0.2595 * σ_total))

    # No-wind, no-slope propagating flux
    push!(eqs, I_R_ξ ~ I_R * ξ)

    # Wind factor
    # ϕ_w = C * U^B * (β/β_opt)^(-E)
    C_wind = 7.47
    B_wind = 0.02526
    E_wind = 0.715
    push!(eqs, ϕ_w ~ C_wind * exp(-0.133 * σ_total^0.55) *
                     U^(B_wind * σ_total^0.54) *
                     β_ratio^(-E_wind * exp(-3.59e-4 * σ_total)))

    # Slope factor
    push!(eqs, ϕ_s ~ 5.275 * β^(-0.3) * tan_ϕ^2)

    # Effective heating number (use characteristic σ)
    push!(eqs, ε ~ exp(-138 / σ_total))

    # Heat of preignition for each size class
    for i in 1:n_dead
        push!(eqs, Q_ig_d[i] ~ 250 + 1116 * M_f_d[i])
    end
    for i in 1:n_live
        push!(eqs, Q_ig_l[i] ~ 250 + 1116 * M_f_l[i])
    end

    # Heat sink (sum over all size classes, weighted by effective heating and fuel load)
    push!(eqs, heat_sink ~ ρ_b * ε * sum(f_d[i] * Q_ig_d[i] * exp(-138/σ_d[i]) for i in 1:n_dead) +
                           ρ_b * ε * sum(f_l[i] * Q_ig_l[i] * exp(-138/σ_l[i]) for i in 1:n_live))

    # No-wind, no-slope rate of spread
    push!(eqs, R_0 ~ I_R_ξ / heat_sink)

    # Final rate of spread
    push!(eqs, R ~ I_R_ξ * (1 + ϕ_w + ϕ_s) / heat_sink)

    # Create and return the system
    @named rothermel_hetero = ODESystem(eqs, t)
    return structural_simplify(rothermel_hetero)
end

# Example usage is provided in examples/heterogeneous_example.jl
