"""
Example usage of the Rothermel Surface Fire Spread Model

This file demonstrates how to use the WildlandFire.jl package
for various fuel types and conditions.
"""

using WildlandFire
using OrdinaryDiffEqDefault
using Plots

"""
Example 1: Basic fire spread calculation for a single fuel type
"""
function example_basic_fire_spread()
    println("=" ^ 80)
    println("Example 1: Basic Fire Spread Calculation")
    println("=" ^ 80)

    # Parameter values for a typical grass fuel (similar to NFFL Fuel Model 1)
    params = [
        h => 8000.0,          # Low heat content (Btu/lb)
        S_T => 0.0555,        # Total mineral content (5.55%)
        S_e => 0.010,         # Effective mineral content (1%)
        ρ_p => 32.0,          # Particle density (lb/ft³)
        σ => 3500.0,          # Surface-area-to-volume ratio (ft²/ft³) - fine grass
        w_o => 0.138,         # Oven-dry fuel load (lb/ft²)
        δ => 1.0,             # Fuel bed depth (ft)
        M_x => 0.12,          # Moisture of extinction (12%)
        M_f => 0.05,          # Fuel moisture content (5%)
        U => 352.0,           # Wind velocity (ft/min) ≈ 4 mph
        tan_ϕ => 0.0,         # Slope (0 = flat ground)
    ]

    # Create problem (steady-state calculation)
    u0 = []
    tspan = (0.0, 1.0)
    prob = ODEProblem(rothermel_simplified, u0, tspan, params)
    sol = solve(prob)

    # Extract results
    R_ftmin = sol[R][end]
    R_mmin = R_ftmin * 0.3048
    R_mhr = R_mmin * 60
    R_chhr = R_ftmin * 1.829  # chains per hour (1 chain = 66 ft)

    println("\nResults:")
    println("  Rate of spread: $(round(R_ftmin, digits=2)) ft/min")
    println("  Rate of spread: $(round(R_mmin, digits=3)) m/min")
    println("  Rate of spread: $(round(R_mhr, digits=2)) m/hr")
    println("  Rate of spread: $(round(R_chhr, digits=2)) chains/hr")
    println("\nIntermediate values:")
    println("  Reaction intensity: $(round(sol[I_R][end], digits=2)) Btu/ft²·min")
    println("  Packing ratio: $(round(sol[β][end], digits=4))")
    println("  Wind factor ϕ_w: $(round(sol[ϕ_w][end], digits=3))")
    println("  Propagating flux ratio ξ: $(round(sol[ξ][end], digits=4))")

    return sol
end

"""
Example 2: Sensitivity analysis - effect of wind speed on rate of spread
"""
function example_wind_sensitivity()
    println("\n" * "=" ^ 80)
    println("Example 2: Wind Speed Sensitivity Analysis")
    println("=" ^ 80)

    # Base parameters
    base_params = Dict(
        h => 8000.0,
        S_T => 0.0555,
        S_e => 0.010,
        ρ_p => 32.0,
        σ => 3500.0,
        w_o => 0.138,
        δ => 1.0,
        M_x => 0.12,
        M_f => 0.05,
        tan_ϕ => 0.0,
    )

    # Wind speeds in mph and ft/min
    wind_speeds_mph = 0:1:15
    wind_speeds_ftmin = wind_speeds_mph .* 88.0  # 1 mph ≈ 88 ft/min

    ros_results = Float64[]

    for U_val in wind_speeds_ftmin
        params = [base_params..., U => U_val]
        prob = ODEProblem(rothermel_simplified, [], (0.0, 1.0), params)
        sol = solve(prob)
        push!(ros_results, sol[R][end] * 1.829)  # Convert to chains/hr
    end

    println("\nWind Speed (mph) | Rate of Spread (chains/hr)")
    println("-" ^ 45)
    for (wind_mph, ros) in zip(wind_speeds_mph, ros_results)
        println("  $(lpad(wind_mph, 3))           | $(round(ros, digits=2))")
    end

    # Plot results
    plot(wind_speeds_mph, ros_results,
         xlabel="Wind Speed (mph)",
         ylabel="Rate of Spread (chains/hr)",
         title="Wind Speed vs Rate of Spread",
         marker=:circle,
         linewidth=2,
         legend=false,
         grid=true)

    savefig("wind_sensitivity.png")
    println("\nPlot saved as wind_sensitivity.png")

    return wind_speeds_mph, ros_results
end

"""
Example 3: Sensitivity analysis - effect of fuel moisture on rate of spread
"""
function example_moisture_sensitivity()
    println("\n" * "=" ^ 80)
    println("Example 3: Fuel Moisture Sensitivity Analysis")
    println("=" ^ 80)

    base_params = Dict(
        h => 8000.0,
        S_T => 0.0555,
        S_e => 0.010,
        ρ_p => 32.0,
        σ => 3500.0,
        w_o => 0.138,
        δ => 1.0,
        M_x => 0.12,
        U => 352.0,  # 4 mph wind
        tan_ϕ => 0.0,
    )

    # Moisture content from 1% to 11% (below extinction)
    moisture_values = 0.01:0.01:0.11

    ros_results = Float64[]

    for M_f_val in moisture_values
        params = [base_params..., M_f => M_f_val]
        prob = ODEProblem(rothermel_simplified, [], (0.0, 1.0), params)
        sol = solve(prob)
        push!(ros_results, sol[R][end] * 1.829)  # Convert to chains/hr
    end

    println("\nMoisture (%) | Rate of Spread (chains/hr)")
    println("-" ^ 45)
    for (moisture, ros) in zip(moisture_values .* 100, ros_results)
        println("  $(lpad(round(moisture, digits=1), 4))       | $(round(ros, digits=2))")
    end

    # Plot results
    plot(moisture_values .* 100, ros_results,
         xlabel="Fuel Moisture Content (%)",
         ylabel="Rate of Spread (chains/hr)",
         title="Fuel Moisture vs Rate of Spread",
         marker=:circle,
         linewidth=2,
         legend=false,
         grid=true)

    savefig("moisture_sensitivity.png")
    println("\nPlot saved as moisture_sensitivity.png")

    return moisture_values, ros_results
end

"""
Example 4: Effect of slope on rate of spread
"""
function example_slope_effect()
    println("\n" * "=" ^ 80)
    println("Example 4: Slope Effect on Rate of Spread")
    println("=" ^ 80)

    base_params = Dict(
        h => 8000.0,
        S_T => 0.0555,
        S_e => 0.010,
        ρ_p => 32.0,
        σ => 3500.0,
        w_o => 0.138,
        δ => 1.0,
        M_x => 0.12,
        M_f => 0.05,
        U => 0.0,  # No wind to isolate slope effect
    )

    # Slope from 0% to 50%
    slope_percent = 0:5:50
    slope_tan = slope_percent ./ 100

    ros_results = Float64[]

    for tan_ϕ_val in slope_tan
        params = [base_params..., tan_ϕ => tan_ϕ_val]
        prob = ODEProblem(rothermel_simplified, [], (0.0, 1.0), params)
        sol = solve(prob)
        push!(ros_results, sol[R][end] * 1.829)  # Convert to chains/hr
    end

    println("\nSlope (%) | Rate of Spread (chains/hr)")
    println("-" ^ 45)
    for (slope, ros) in zip(slope_percent, ros_results)
        println("  $(lpad(slope, 3))     | $(round(ros, digits=2))")
    end

    # Plot results
    plot(slope_percent, ros_results,
         xlabel="Slope (%)",
         ylabel="Rate of Spread (chains/hr)",
         title="Slope vs Rate of Spread (No Wind)",
         marker=:circle,
         linewidth=2,
         legend=false,
         grid=true)

    savefig("slope_effect.png")
    println("\nPlot saved as slope_effect.png")

    return slope_percent, ros_results
end

"""
Example 5: Different fuel types (NFFL fuel models)
"""
function example_fuel_types()
    println("\n" * "=" ^ 80)
    println("Example 5: Comparison of Different Fuel Types")
    println("=" ^ 80)

    # Standard environmental conditions
    standard_conditions = Dict(
        M_f => 0.08,      # 8% moisture
        U => 352.0,       # 4 mph wind
        tan_ϕ => 0.0,     # Flat ground
    )

    # Define several fuel types (simplified NFFL models)
    fuel_types = Dict(
        "Grass (FM 1)" => Dict(
            h => 8000.0, S_T => 0.0555, S_e => 0.010, ρ_p => 32.0,
            σ => 3500.0, w_o => 0.138, δ => 1.0, M_x => 0.12
        ),
        "Timber Grass (FM 2)" => Dict(
            h => 8000.0, S_T => 0.0555, S_e => 0.010, ρ_p => 32.0,
            σ => 3000.0, w_o => 0.092, δ => 1.0, M_x => 0.15
        ),
        "Tall Grass (FM 3)" => Dict(
            h => 8000.0, S_T => 0.0555, S_e => 0.010, ρ_p => 32.0,
            σ => 1500.0, w_o => 0.138, δ => 2.5, M_x => 0.25
        ),
        "Chaparral (FM 4)" => Dict(
            h => 8000.0, S_T => 0.0555, S_e => 0.010, ρ_p => 32.0,
            σ => 2000.0, w_o => 0.230, δ => 6.0, M_x => 0.20
        ),
        "Brush (FM 5)" => Dict(
            h => 8000.0, S_T => 0.0555, S_e => 0.010, ρ_p => 32.0,
            σ => 1750.0, w_o => 0.046, δ => 2.0, M_x => 0.20
        ),
    )

    println("\nFuel Type            | Rate of Spread (chains/hr)")
    println("-" ^ 55)

    for (fuel_name, fuel_params) in fuel_types
        params = [fuel_params..., standard_conditions...]
        prob = ODEProblem(rothermel_simplified, [], (0.0, 1.0), params)
        sol = solve(prob)
        ros = sol[R][end] * 1.829
        println("  $(rpad(fuel_name, 20)) | $(round(ros, digits=2))")
    end
end

# Run all examples
println("\n")
println("╔" * "═"^78 * "╗")
println("║" * " "^20 * "ROTHERMEL FIRE SPREAD MODEL EXAMPLES" * " "^22 * "║")
println("╚" * "═"^78 * "╝")

example_basic_fire_spread()
example_wind_sensitivity()
example_moisture_sensitivity()
example_slope_effect()
example_fuel_types()

println("\n" * "=" ^ 80)
println("All examples completed successfully!")
println("=" ^ 80 * "\n")
