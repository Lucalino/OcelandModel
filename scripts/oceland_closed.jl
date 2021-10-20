# OCELAND MODEL - CLOSED VERSION

# written by Luca Schmidt
# started 29th September 2020

# The oceland model comprises two domains, an ocean and a land part, 
# representing a tropical island.

# The aim of the model is to constrain the ratio of precipitation over 
# land and ocean, P_l/P_o, based on simple equations of water vapor 
# and soil moisture conservation. 

using DrWatson
@quickactivate "Oceland Model"
using BenchmarkTools
using CairoMakie
using Colors
using CSV
using DataFrames
using DifferentialEquations
using Distributions
using DynamicalSystems
using IntervalRootFinding
using NLsolve
#using OrdinaryDiffEq
using Random
using StaticArrays
include(srcdir("parametrisations.jl"))
include(srcdir("utils.jl"))
include(srcdir("cm_analysis.jl"))
include(srcdir("model_versions.jl"))


fs=14.0
lw=2.0

# Model versions:
# 1  : soil moisture ODE without parametrised precipitation but given advection Pa
# 2  : system solved with NLsolve.jl as Monte-Carlo runs. Precipitaiton is parametrised, E_l is piecewise defined
# 3  : set of ODEs for the three state variables s, wl, wo. Precipitation is parametrised, E_l is piecewise defined
# 4  : system solved with DynamicalSystems.jl package and piecewise defined E_l
# 5  : system solved with DynamicalSystems.jl package and tanh-version of E_l
# 6  : system solved with IntervalRootFinding.jl. Precipitation parametrised, E_l piecewise defined

calc_mode = 6

params = cm_rand_params()


if calc_mode == 1

    fig = Figure()
    ax  = Axis(fig[1,1], xlabel = L"time $t$ [days]", ylabel = L"soil moisture saturation $s$")

    for i = 1:3

        #Solve the ODE for s(t)
        p = cm_fixed_params(i)
        s0 = 0.3 #initial condition
        tspan = (0.0,150.0) #in days
        prob = ODEProblem(soilmoisture,s0,tspan,p)
        sol = solve(prob, reltol=1e-6,saveat=1.0)
        @unpack pa = p
        #Plot s(t) over t
        lines!(sol.t,sol.u,label= "Pa = $pa")

    end

    ax.xlabelsize = fs
    ax.ylabelsize = fs
    hidespines!(ax, :t, :r)
    axislegend(position = :lt, labelsize = 16, framevisible = false)
    fig
    #savefig(plotsdir(string("Closed model/Dynamics/sm_Pa0.5-0.9.pdf")))


    #Compute Pl(t), Po(t) and PR(t) using the computed s(t)
    # Pl = zeros(length(sol.t))
    # Po = zeros(length(sol.t))
    # PR = zeros(length(sol.t))

    # for t = 1:length(sol.u)
    #     Pl[t] = land_evapo(spwp,sfc,ep,sol.u[t]) + pa
    #     Po[t] = (ptot - α * Pl[t])/(1-α)
    #     PR[t] = Pl[t]/Po[t]
    # end

    # #Fluxes
    # Pl_fl = Pl * α
    # Po_fl = Po * (1-α)
    # PR_fl = Pl_fl ./ Po_fl

elseif calc_mode == 2

    ### VERSION 2 ###
    #with parametrised precipitation and s, wl, wo as state variables of the system

    # Defining the system of equations

    function f!(F, x, p)

        @unpack eo, α, nZr, u, L = params

        #Variable naming: x[1] = s, x[2] = wl, x[3] = wo
        F[1] = 1/(nZr) * (precip(x[2], params) * infiltration(x[1], params) - land_evap(x[1], params))
        F[2] = land_evap(x[1], params) - precip(x[2], params) + (x[3] - x[2]) * u / (α*L)
        F[3] = eo - precip(x[3], params) - (x[3] - x[2]) * u / ((1-α) * L)
    end


    function cm_monte_carlo_run(nb_runs)

        param_names = keys(params)
        col_names = [string(el) for el in param_names]
        col_names = append!(col_names, ["s", "wl", "wo", "convergence"])
        output_columns = length(col_names)

        #initialising output matrix
        sol = Array{Float64}(undef, 0, output_columns)
    
        for i=1:nb_runs
            init_cond = [0.5, 50.0, 50.0]
            sol = [sol; cm_eq_solution(f!, init_cond)]
        end
    
        sol_df = DataFrame(sol, col_names)
        #CSV.write(datadir("sims", "cm_eq_MonteCarlo_scan_$(nb_runs)_runs.csv"), sol_df)
        println(sol_df)
    end

    cm_monte_carlo_run(5)
    



elseif calc_mode == 3

    tspan = (0.0, 100.0)
    x0 = [0.6, 30.0, 60.0]
    prob = ODEProblem(closed_model_pw, x0, tspan, params)
    sol = solve(prob)

    st  = sol.u[1][1]
    wlt = sol.u[1][2]
    wot = sol.u[1][3]

    for n = 2:length(sol.t)
        st  = [st; sol.u[n][1]]
        wlt = [wlt; sol.u[n][2]]
        wot = [wot; sol.u[n][3]]
    end

    cm_t_evolution_plot(st, wlt, wot, sol.t)


elseif calc_mode == 4

    x0 = @SVector [0.3, 40.0, 60.0] 
    dynsys = ContinuousDynamicalSystem(closed_model_piecewise, x0, params)
    box = cm_state_space(params)
    fp, eigs, stable = fixedpoints(dynsys, box)

elseif calc_mode == 5

    x0 = @SVector [0.3, 60.0, 50.0] 
    dynsys = ContinuousDynamicalSystem(closed_model_smooth, x0, params)
    box = cm_state_space(params)
    fp, eigs, stable = fixedpoints(dynsys, box)
    # tspan = (0.0,100.0)
    # prob = ODEProblem(closed_model_smooth, x0, tspan, params)
    # solve(prob, reltol=1e-8)
    #@benchmark solve(prob,Vern9(), save_every_step = false)    
    #@btime closed_model_smooth($(x0), $params, 0.0)  
    
elseif calc_mode == 6

    closure_pw = x -> closed_model_piecewise(x, params, 0.0)
    box = cm_state_space(params)
    rts = IntervalRootFinding.roots(closure_pw, box)
    println(rts)

else
    println("Specify which version of model formulation you want to work with")
end