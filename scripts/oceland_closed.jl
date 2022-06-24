# OCELAND MODEL - CLOSED VERSION

# written by Luca Schmidt
# started 29th September 2020

using DrWatson
@quickactivate "Oceland Model"
using CSV
using DataFrames
using Distributions
using DynamicalSystems
using Random
using StaticArrays
include(srcdir("parametrisations.jl"))
include(srcdir("utils.jl"))
include(srcdir("model.jl"))
include(srcdir("model_versions.jl"))
include(srcdir("sensitivity_analysis.jl"))


"""
    cm_compute_EQ_data(nb_runs::Int, tau::Bool=true))

Compute the fixed points of the closed Oceland model for a number (nb_runs) of 
combinations of parameter values that are randomly chosen from ranges defined in 
cm_rand_params().
The input parameter tau determines whether the equations are formulated with τ (tau = true) 
or with u and L individually (tau = false).
The function returns a DataFrame with parameter values and corresponding equilibrium
values for soil moisture s and land and ocean water vapor passes.
"""
function cm_compute_EQ_data(nb_runs::Int, tau::Bool=true)

    col_names = [string(el) for el in keys(cm_rand_params(tau))]
    col_names = append!(col_names, ["s", "wl", "wo"])
    sol = Array{Float64}(undef, 0, length(col_names))
    x0 = @SVector [0.5, 50.0, 50.0]

    for n = 1:nb_runs
        sol = [sol; cm_equilibrium_solution(x0, tau)]
    end

    sol_df = DataFrame(sol, col_names)

    if tau == true
        ps = "tau"
    else
        ps = "uL"
    end

    #CSV.write(datadir("sims", "closed model data", "Final runs", "cm_$(ps)_fixedpoints_$(nb_runs)runs.csv"), sol_df)
    
    return sol_df

end