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
using DifferentialEquations
using Distributions
using DynamicalSystems
using StaticArrays
using Random
include(srcdir("model_versions.jl"))
include(srcdir("parametrisations.jl"))
include(srcdir("utils.jl"))
include(srcdir("om_analysis.jl"))


function om_MC_fixedpoints(nb_runs, system, tau::Bool=false)

    col_names = [string(el) for el in keys(om_rand_params())]
    col_names = append!(col_names, ["s", "w1", "w2", "w3", "status"])
    sol = Array{Float64}(undef, 0, length(col_names))
    x0 = @SVector [0.5, 50.0, 50.0, 50.0]

    for n = 1:nb_runs
        sol = [sol; om_fixedpoints(x0, system)]
    end

    sol_df = DataFrame(sol, col_names)
    #d = Int(round(mean(sol_df.L) .* mm2km(1.0), digits = 1))
    CSV.write(datadir("sims", "open model pmscan", "om_$(system)_MC_fixedpoints_runs$(nb_runs)_sym_updatedparams.csv"), sol_df)
end

om_MC_fixedpoints(10000, "v2") 

# p = cm_rand_params()
# x0 = @SVector [0.6, 40.0, 40.0]
# dynsys = ContinuousDynamicalSystem(closed_model_smooth, x0, p)
# diffeq = (alg = Vern9(), adaptive = false, dt = 0.001, reltol = 1e-8, abstol = 1e-8)
# sg  = range(0.0, 1.0; length = 100)
# wlg = wog = range(0.0, p[:wsat]; length = 100)
# basins, attractors = basins_of_attraction((sg, wlg, wog), dynsys)
#fig = cm_basins_plot(basins[1,:,:], wlg, wog)

#@btime closed_model_smooth($(x0), $params, 0.0)  