# OCELAND MODEL - OPEN VERSION

# written by Luca Schmidt
# started 16th November 2021

# The open oceland model comprises three domains, two ocean parts and a land part in the middle, 
# representing a tropical island.

# The aim of the model is to constrain the ratio of precipitation over 
# land and ocean, P_l/P_o, based on simple equations of water vapor 
# and soil moisture conservation. 

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
include(srcdir("open_model.jl"))


function om_compute_EQ_data(nb_runs, model_version, tau::Bool=false)

    col_names = [string(el) for el in keys(om_rand_params())]
    col_names = append!(col_names, ["s", "w1", "w2", "w3", "status"])
    sol = Array{Float64}(undef, 0, length(col_names))
    x0 = @SVector [0.5, 50.0, 50.0, 50.0]

    for n = 1:nb_runs
        sol = [sol; om_equilibrium_solution(x0, model_version)]
    end

    sol_df = DataFrame(sol, col_names)
    #d = Int(round(mean(sol_df.L) .* mm2km(1.0), digits = 1))
    return sol_df
    #CSV.write(datadir("sims", "open model pmscan", "om_$(system)_MC_fixedpoints_runs$(nb_runs)_sym_updatedparams-5.csv"), sol_df)
end

#om_compute_EQ_data(10000, "v2") 