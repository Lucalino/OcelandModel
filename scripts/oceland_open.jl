# OCELAND MODEL - OPEN VERSION

# written by Luca Schmidt
# started 11th May 2021

# The open oceland model comprises three domains, two ocean parts and a land part in the middle, 
# representing a tropical island.

# The aim of the model is to constrain the ratio of precipitation over 
# land and ocean, P_l/P_o, based on simple equations of water vapor 
# and soil moisture conservation. 

using DrWatson
@quickactivate "Oceland Model"
#using PyPlot
#pygui(true)
using DifferentialEquations
using IntervalRootFinding
using IntervalArithmetic
using StaticArrays
using Random
using Distributions
using DataFrames
using BenchmarkTools
using CSV
include(srcdir("parametrisations.jl"))
include(srcdir("utils.jl"))
include(srcdir("analysis.jl"))


##### EQUILIBRIUM SYSTEM ######

function f_one(variables, parameters, t)
    w1, w2, w3, s = variables
    @unpack spwp, sfc, Ep, eo, ϵ, r, nZr, w0, L, Li, Lo1, Lo2, u, a, b, w_sat = parameters
    return SVector(2 * u * (w0 - w1) / Lo1 + eo - exp(a * (w1/w_sat - b)), 
                   2 * u * (2 * w1 - w2 - w0) / Li - exp(a * (w2/w_sat - b)),
                   2 * u * (2 * w2 - 2 * w1 - w3 + w0) / Lo2 + eo - exp(a * (w3/w_sat - b)),
                   exp(a * (w2/w_sat - b)) * infiltration(s, ϵ, r)
                   )
end


function f_two(variables, parameters, t)
    w1, w2, w3, s = variables
    @unpack spwp, sfc, Ep, eo, ϵ, r, nZr, w0, L, Li, Lo1, Lo2, u, a, b, w_sat = parameters
    return SVector(2 * u * (w0 - w1) / Lo1 + eo - exp(a * (w1/w_sat - b)), 
                   2 * u * (2 * w1 - w2 - w0) / Li + Ep/(sfc - spwp) * (s - spwp) - exp(a * (w2/w_sat - b)),
                   2 * u * (2 * w2 - 2 * w1 - w3 + w0) / Lo2 + eo - exp(a * (w3/w_sat - b)),
                   exp(a * (w2/w_sat - b)) * infiltration(s, ϵ, r) - Ep/(sfc - spwp) * (s - spwp)
                   )
end


function f_three(variables, parameters, t)
    w1, w2, w3, s = variables #this line just gives arguments in the right order, "variables" can be any iterable object
    @unpack spwp, sfc, Ep, eo, ϵ, r, nZr, w0, L, Li, Lo1, Lo2, u, a, b, w_sat = parameters
    return SVector(2 * u * (w0 - w1) / Lo1 + eo - exp(a * (w1/w_sat - b)), 
                   2 * u * (2 * w1 - w2 - w0) / Li + Ep - exp(a * (w2/w_sat - b)),
                   2 * u * (2 * w2 - 2 * w1 - w3 + w0) / Lo2 + eo - exp(a * (w3/w_sat - b)),
                   exp(a * (w2/w_sat - b)) * infiltration(s, ϵ, r) - Ep
                   )
end


##### RANDOM SAMPLING AND RESPECTIVE SOLUTIONS ######

function monte_carlo_run(nb_runs, w0_fixed = Nothing)

    output_columns = 20
    col_names = ["spwp", "sfc", "Ep", "eo", "ϵ", "r", "nZr", 
                 "w0", "L", "Li", "Lo1", "Lo2", "u", "a", "b",
                 "w_sat", "w1", "w2", "w3", "s"]

    #initialising output matrix
    sol = Array{Float64}(undef, 0, output_columns)

    for i=1:nb_runs
        params = rand_params(w0_fixed)
        sol = [sol; om_eq_solution(f_one, f_two, f_three, params, output_columns)]
    end

    sol_df = DataFrame(sol, col_names)
    CSV.write(datadir("sims", "om_eq_MonteCarlo_scan_$(nb_runs)_runs.csv"), sol_df)
end

#@btime monte_carlo_run(100000)
monte_carlo_run(100000)

