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
include(srcdir("cm_plotting.jl"))


# Model versions:
# 1  : soil moisture ODE without parametrised precipitation but given advection Pa
# 2  : system solved with NLsolve.jl as Monte-Carlo runs. Precipitaiton is parametrised, E_l is piecewise defined
# 3  : set of ODEs for the three state variables s, wl, wo. Precipitation is parametrised, E_l is piecewise defined
# 4  : system solved with DynamicalSystems.jl package with smooth OR piecewise defined E_l
# 5  : system solved with IntervalRootFinding.jl. Precipitation parametrised, E_l piecewise defined

calc_mode = 40


if calc_mode == 1

    fs= 14.0 #fontsize
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

    p = cm_rand_params()

    function f!(F, x, p)

        @unpack eo, α, nZr, u, L = p

        #Variable naming: x[1] = s, x[2] = wl, x[3] = wo
        F[1] = 1/(nZr) * (precip(x[2], p) * infiltration(x[1], p) - land_evap(x[1], p))
        F[2] = land_evap(x[1], p) - precip(x[2], p) + (x[3] - x[2]) * u / (α*L)
        F[3] = eo - precip(x[3], p) - (x[3] - x[2]) * u / ((1-α) * L)
    end


    function cm_MC_nlsolve(nb_runs)

        p_names = keys(p)
        col_names = [string(el) for el in p_names]
        col_names = append!(col_names, ["s", "wl", "wo", "convergence"])
        output_columns = length(col_names)

        #initialising output matrix
        sol = Array{Float64}(undef, 0, output_columns)
        x0 = [0.5, 50.0, 50.0]
    
        for i=1:nb_runs            
            sol = [sol; cm_eq_nlsolve(f!, x0)]
        end
    
        sol_df = DataFrame(sol, col_names)
        #CSV.write(datadir("sims", "closed model pmscan", "cm_piecewise_eq_MC_nlsolve_$(nb_runs)_runs.csv"), sol_df)
        println(sol_df)
    end

    #cm_MC_nlsolve(5)
    

elseif calc_mode == 3

    p = cm_fixed_params(1)

    tspan = (0.0, 100.0)
    x0 = [0.3, 60.0, 40.0]
    prob = ODEProblem(closed_model_pw, x0, tspan, p)
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

    function cm_MC_fixedpoints(nb_runs, system)

        col_names = [string(el) for el in keys(cm_rand_params())]
        col_names = append!(col_names, ["s", "wl", "wo"])
        sol = Array{Float64}(undef, 0, length(col_names))
        x0 = @SVector [0.5, 50.0, 50.0]

        for n = 1:nb_runs
            sol = [sol; cm_eq_fixedpoints(system, x0)]
        end

        sol_df = DataFrame(sol, col_names)
        CSV.write(datadir("sims", "closed model pmscan", "cm_$(system)_eq_MC_fixedpoints_$(nb_runs)_runs_domain1000.csv"), sol_df)
        #println(sol_df)

    end

    cm_MC_fixedpoints(10000, "smooth")

    # p = cm_rand_params()
    # x0 = @SVector [0.6, 40.0, 40.0]
    # dynsys = ContinuousDynamicalSystem(closed_model_smooth, x0, p)
    # diffeq = (alg = Vern9(), adaptive = false, dt = 0.001, reltol = 1e-8, abstol = 1e-8)
    # sg  = range(0.0, 1.0; length = 100)
    # wlg = wog = range(0.0, p[:wsat]; length = 100)
    # basins, attractors = basins_of_attraction((sg, wlg, wog), dynsys)
    #fig = cm_basins_plot(basins[1,:,:], wlg, wog)

    #@btime closed_model_smooth($(x0), $params, 0.0)  

    
elseif calc_mode == 5

    closure_pw = x -> closed_model_piecewise(x, params, 0.0)
    box = cm_state_space(params)
    rts = IntervalRootFinding.roots(closure_pw, box)
    println(rts)

else
    println("Specify a valid version of the model formulation that you want to work with!")
end