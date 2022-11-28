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
using DifferentialEquations
using Statistics
include(srcdir("parametrizations.jl"))
include(srcdir("utils.jl"))
include(srcdir("create_model_output.jl"))
include(srcdir("model_versions.jl"))
include(srcdir("sensitivity_analysis.jl"))

#USED FOR PAPER
"""
    cm_compute_EQ_data(nb_runs::Int, tau::Bool=true))

Compute the fixed points of the closed Oceland model, which is purely based on water balance equations without
a diurnal cycle, for a number (nb_runs) of combinations of parameter values that are randomly chosen from ranges 
defined in cm_rand_params().
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

    #CSV.write(datadir("sims", "closed model", "test_$(nb_runs)runs.csv"), sol_df)
    
    return sol_df

end


#Computes analytical solutions to the model equations with linear fluxes 
#for nb_runs randomly chosen combinations of parameter values. 
function cm_compute_EQ_data_lin(nb_runs::Int)
    col_names = [string(el) for el in keys(cm_lin_rand_params())]
    col_names = append!(col_names, ["s", "wl", "wo"])
    sol = Array{Float64}(undef, 0, length(col_names))
    for n = 1:nb_runs
        p = cm_lin_rand_params()
        s, wl, wo = cm_lin(p)
        solrow = [transpose(collect(values(p))) s wl wo]
        sol = [sol; solrow]
    end
    return DataFrame(sol, col_names)     

end


#USED IN SUBSEQUENT PROJECT
#solution computed directly with DifferentialEquations.jl
function cm_DC_diffeq(t_length::Float64)
    u0 = [0.0; 0.0; 0.0]
    tspan = (0.0, t_length)
    p = cm_fixed_params(false)
    prob = ODEProblem(closed_model_DC_wind!, u0, tspan, p)
    sol = solve(prob, Tsit5(), reltol = 1e-8, abstol = 1e-8, dtmax = 0.0001)
    d = cm_DC_derived_quantities(sol, p)
    return d
end

#solution computed with DynamicalSystems.jl (adds more functionality)
function cm_DC_dynsys(t_length::Float64=1000.0)
    x0 = [0.0, 0.0, 0.0]
    p = cm_fixed_params(false)
    ds = ContinuousDynamicalSystem(closed_model_DC, x0, p)
    tr = trajectory(ds, t_length)
    t = collect(0.0:0.01:t_length)
    raw_data = hcat(t, tr)
    data = cm_DC_derived_quantities(raw_data, p)
    eq_mean_state = cm_DC_EQmeanvalues(ds, x0, p)
    return data, eq_mean_state
end


function one_at_a_time_sensitivity(param::Symbol, llim::Float64, ulim::Float64, t_length = 1000.0)
    x0 = [0.5, 50.0, 50.0]
    p = cm_fixed_params(false)
    ds = ContinuousDynamicalSystem(closed_model_DC, x0, p)
    EQstates = Array{Float64}(undef, 0, 3)
    param_range = collect(llim:(ulim-llim)/100:ulim)
    for n in param_range
        p[param] = n
        EQstates = [EQstates; transpose(cm_DC_EQmeanvalues(ds, x0, p))]
    end
    data = hcat(param_range, EQstates)
    return data
end