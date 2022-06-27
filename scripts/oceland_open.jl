# OCELAND MODEL - OPEN VERSION

using DrWatson
@quickactivate "Oceland Model"
using BenchmarkTools
using CairoMakie
using CSV
using DataFrames
using Distributions
using DynamicalSystems
using StaticArrays
using Random
include(srcdir("model_versions.jl"))
include(srcdir("parametrisations.jl"))
include(srcdir("utils.jl"))
include(srcdir("create_model_output.jl"))


"""
    om_compute_EQ_data(nb_runs::Int, model_version::String)

Computes the equilibrium solution to one of the open model versions 
(model_verison = paper, closed or w_tracked) for a number (nb_runs) of different randomly chosen
combinations of parameter values.
"""
function om_compute_EQ_data(nb_runs::Int, model_version::String = "paper")

    col_names = [string(el) for el in keys(om_rand_params())]
    col_names = append!(col_names, ["s", "w1", "w2", "w3", "status"])
    sol = Array{Float64}(undef, 0, length(col_names))
    x0 = @SVector [0.5, 50.0, 50.0, 50.0]

    for n = 1:nb_runs
        sol = [sol; om_equilibrium_solution(x0, model_version)]
    end

    sol_df = DataFrame(sol, col_names)
    CSV.write(datadir("sims", "open model", "test_$(nb_runs)runs.csv"), sol_df)
    return sol_df

end



"""
    om_likelihood_moistening_scenarios(d::DataFrame)

Computes the percentage of open model runs (different parameter combinations) that exhibit a 
certain atmospheric moistening behavior, represented by different branches of a tree diagram. 
Branching in the diagram corresponds to whether net advection into the next atmospheric
box is positive (subdomain dries the air) or negative (subdomain moistens the air). 
E.g. branch one is the case in which the atmosphere dries monotically in wind direction.
"""
function om_likelihood_moistening_scenarios(d::DataFrame)

    lh1 = nrow(d[(d.Δw1 .> 0) .& (d.Δw2 .> 0) .& (d.Δw3 .> 0),:]) * 100/ nrow(d)
    lh11 = nrow(d[(d.Δw1 .> 0) .& (d.Δw2 .> 0) .& (d.Δw3 .> 0) .& (d.Δwtot .> 0),:]) * 100/ nrow(d[(d.Δw1 .> 0) .& (d.Δw2 .> 0) .& (d.Δw3 .> 0),:])
    lh12 = nrow(d[(d.Δw1 .> 0) .& (d.Δw2 .> 0) .& (d.Δw3 .> 0) .& (d.Δwtot .< 0),:]) * 100/ nrow(d[(d.Δw1 .> 0) .& (d.Δw2 .> 0) .& (d.Δw3 .> 0),:])
    lh2 = nrow(d[(d.Δw1 .> 0) .& (d.Δw2 .> 0) .& (d.Δw3 .< 0),:]) * 100/ nrow(d)
    lh21 = nrow(d[(d.Δw1 .> 0) .& (d.Δw2 .> 0) .& (d.Δw3 .< 0) .& (d.Δwtot .> 0),:]) * 100/ nrow(d[(d.Δw1 .> 0) .& (d.Δw2 .> 0) .& (d.Δw3 .< 0),:])
    lh22 = nrow(d[(d.Δw1 .> 0) .& (d.Δw2 .> 0) .& (d.Δw3 .< 0) .& (d.Δwtot .< 0),:]) * 100/ nrow(d[(d.Δw1 .> 0) .& (d.Δw2 .> 0) .& (d.Δw3 .< 0),:])
    lh3 = nrow(d[(d.Δw1 .> 0) .& (d.Δw2 .< 0) .& (d.Δw3 .> 0),:]) * 100/ nrow(d)
    lh4 = nrow(d[(d.Δw1 .> 0) .& (d.Δw2 .< 0) .& (d.Δw3 .< 0),:]) * 100/ nrow(d)
    lh5 = nrow(d[(d.Δw1 .< 0) .& (d.Δw2 .> 0) .& (d.Δw3 .> 0),:]) * 100/ nrow(d)
    lh6 = nrow(d[(d.Δw1 .< 0) .& (d.Δw2 .> 0) .& (d.Δw3 .< 0),:]) * 100/ nrow(d)
    lh61 = nrow(d[(d.Δw1 .< 0) .& (d.Δw2 .> 0) .& (d.Δw3 .< 0) .& (d.Δwtot .> 0),:]) * 100/ nrow(d[(d.Δw1 .< 0) .& (d.Δw2 .> 0) .& (d.Δw3 .< 0),:])
    lh62 = nrow(d[(d.Δw1 .< 0) .& (d.Δw2 .> 0) .& (d.Δw3 .< 0) .& (d.Δwtot .< 0),:]) * 100/ nrow(d[(d.Δw1 .< 0) .& (d.Δw2 .> 0) .& (d.Δw3 .< 0),:])
    lh7 = nrow(d[(d.Δw1 .< 0) .& (d.Δw2 .< 0) .& (d.Δw3 .> 0),:]) * 100/ nrow(d)
    lh8 = nrow(d[(d.Δw1 .< 0) .& (d.Δw2 .< 0) .& (d.Δw3 .< 0),:]) * 100/ nrow(d)
    lh9 = nrow(d[d.Δwtot .> 0,:]) * 100 / nrow(d)
    lh10 = nrow(d[d.Δwtot .< 0,:]) * 100 / nrow(d)

    println("The likelihood for branch 1 is: $(lh1) %,")
    println("thereof $(lh11) % where Δw_tot > 0 and $(lh12) % where Δw_tot < 0.")
    println("The likelihood for branch 2 is: $(lh2) %.")
    println("thereof $(lh21) % where Δw_tot > 0 and $(lh22) % where Δw_tot < 0.")
    println("The likelihood for branch 3 is: $(lh3) %.")
    println("The likelihood for branch 4 is: $(lh4) %.")
    println("The likelihood for branch 5 is: $(lh5) %.")
    println("The likelihood for branch 6 is: $(lh6) %.")
    println("thereof $(lh61) % where Δw_tot > 0 and $(lh62) % where Δw_tot < 0.")
    println("The likelihood for branch 7 is: $(lh7) %.")
    println("The likelihood for branch 8 is: $(lh8) %.")
    println("The likelihood for Δw_tot > 0 is: $(lh9) %.")
    println("The likelihood for Δw_tot > 0 is: $(lh10) %.")

end