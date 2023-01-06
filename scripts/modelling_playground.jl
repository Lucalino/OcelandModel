using DrWatson
@quickactivate "Oceland Model"
using CSV
using DataFrames
using Distributions
using DynamicalSystems
using Random
using StaticArrays
using DifferentialEquations
#using ModelingToolkit
include(srcdir("parametrizations.jl"))
include(srcdir("utils.jl"))
include(srcdir("create_model_output.jl"))
include(srcdir("model_versions.jl"))
include(srcdir("sensitivity_analysis.jl"))


function cm_DC_wind_xresolved!(du, u, p, t)
    @unpack nZr, eo, α, L = p
    du[1] = (precip(u[2], p) * infiltration(u[1], p) - evap_scaling(t, p) * El_tanh(u[1], p)) / nZr
    du[2] = evap_scaling(t, p) * El_tanh(u[1], p) - precip(u[2], p) + 2 * advected_moisture(u[2], u[3], t, p) * wind_DC(t, p) / (α * L)
    du[3] = eo - precip(u[3], p) - 2 * advected_moisture(u[2], u[3], t, p) * wind_DC(t, p) / ((1-α) * L)
    #a = eo - precip(u[3], p)
    #b = 2 * advected_moisture(u[2], u[3], t, p) * wind_DC(t, p) / ((1-α) * L)
    #@show a, b
end



function cm_DC_diffeq_xt(t_length::Float64)
    tspan = (0.0, t_length)
    pdict = cm_fixed_params(false)
    
    # Need to give parameters as vector in newest version of ODEProblem.
    # Hopefully this gets fixed, so that I can use the dictionary directly again.
    parr = [pdict[:dx], pdict[:x]]
    w0 = fill(40.0, pdict[:nb_boxes])
    #w0 = [40.0 + 10*cos(2*pi * x / pdict[:L]) for x in pdict[:x]]
    prob = ODEProblem(cm_DC_w_xt!, w0, tspan, parr)
    sol = solve(prob, KenCarp4(), reltol = 1e-8, abstol = 1e-8, dtmax = 0.00001) #Vern9(),AutoTsit5(Rosenbrock23()),Tsit5(), KenCarp4(), lsoda()
    return sol
end


   



