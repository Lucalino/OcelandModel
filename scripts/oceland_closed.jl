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
include(srcdir("parametrisations.jl"))
include(srcdir("utils.jl"))
include(srcdir("cm_analysis.jl"))
using BenchmarkTools
using CairoMakie
using ChaosTools
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


fs=14.0
lw=2.0

# Model versions:
# 1  : soil moisture ODE without parametrised precipitation but given advection Pa
# 2  : system solved with NLsolve.jl as Monte-Carlo runs. Precipitaiton is parametrised, E_l is piecewise defined
# 3  : set of ODEs for the three state variables s, wl, wo. Precipitation is parametrised, E_l is piecewise defined
# 4  : system solved with DynamicalSystems.jl package and piecewise defined E_l
# 5  : system solved with DynamicalSystems.jl package and tanh-version of E_l
# 6  : system solved with IntervalRootFinding.jl. Precipitation parametrised, E_l piecewise defined

model_version = 3

# dictionary of fixed parameters
spwp = 0.2  #permanent wilting point
sfc = spwp + 0.3  #field capacity
ep = 4.38    #[mm/day] potential evaporation over land in mm/day, taken from [1]
pt = 10
ptot = 3.0   #[mm/day] average total precipitation rate over the full tropics
eo = 3.0     #[mm/day] ocean evaporation rate
pa = [0.5,0.7,0.9] #[mm/day] advected precipitation component
ϵ = 1.0      #numerical parameter from Rodriguez-Iturbe et al. (1991)
r = 2.0      #numerical parameter from Rodriguez-Iturbe et al. (1991)
α = 0.1      #land fraction
nZr = 100.0  #[mm] reservoir depth/"field storage capacity of the soil"
a = 15.6     #numerical parameter from Bretherton et al. (2004)
b = 0.603    #numerical parameter from Bretherton et al. (2004)
wsat = 72.0  #[mm] saturation water vapour pass derived from plots in Bretherton et al. ( 2004)
u = 5.0 * m2mm(1)/s2day(1) #[mm/day] wind speed
L = 10000 * km2mm(1) #[mm] domain size
w_crit = 25.0  #[mm] critical wvp-value below which P(w)=0

p_fixed = @dict spwp sfc ep ptot eo pa ϵ r α nZr a b wsat u L w_crit pt


function closed_model_piecewise(x, p, t)
    #This format is required by ContinuousDynamicalSystem.jl algorithm

    @unpack spwp, sfc, ep, ptot, eo, ϵ, r, α, nZr, a, b, wsat, u, L, w_crit, pt = p

    ds = 1/(nZr) * (precip(x[2], wsat, a, b) * infiltration(x[1], ϵ, r) - land_evap(x[1], spwp, sfc, ep))
    dwl = land_evap(x[1], spwp, sfc, ep) - precip(x[2], wsat, a, b) + (x[3] - x[2]) * u / (α*L)
    dwo = eo - precip(x[3], wsat, a, b) - (x[3] - x[2]) * u / ((1-α) * L)
    
    return SVector{3}(ds, dwl, dwo)
end

function closed_model_smooth(x, p, t)
    #This format is required by ContinuousDynamicalSystem.jl algorithm
   
    @unpack spwp, sfc, ep, ptot, eo, ϵ, r, α, nZr, a, b, wsat, u, L, w_crit, pt = p

    ds = 1/(nZr) * (precip(x[2], wsat, a, b) * infiltration(x[1], ϵ, r) - evap_tanh(x[1], spwp, sfc, ep, pt))
    dwl = evap_tanh(x[1], spwp, sfc, ep, pt) - precip(x[2], wsat, a, b) + (x[3] - x[2]) * u / (α*L)
    dwo = eo - precip(x[3], wsat, a, b) - (x[3] - x[2]) * u / ((1-α) * L)
    
    return SVector{3}(ds, dwl, dwo)
end


if model_version == 1

    function soilmoisture(s,p,t)
        spwp, sfc, ep, ptot, eo, pa, ϵ, r, α, nZr = p
        #ds =( ( (1-α) * pa + land_evap(s, spwp, sfc, ep) ) * infiltration(s, ϵ,r) - land_evap(s, spwp, sfc, ep) ) / nZr
        ds = (pa * infiltration(s,ϵ,r) + land_evap(s,spwp,sfc,ep) * (infiltration(s,ϵ,r) - 1)) / nZr
    end

    @unpack spwp, sfc, ep, ptot, eo, pa, ϵ, r, α, nZr = p_fixed
    fig = Figure()
    ax  = Axis(fig[1,1], xlabel = L"time $t$ [days]", ylabel = L"soil moisture saturation $s$")

    for i = 1:3

        #Solve the ODE for s(t)
        p = (spwp,sfc,ep,ptot,eo,pa[i],ϵ,r,α,nZr)
        s0 = 0.3 #initial condition
        tspan = (0.0,150.0) #in days
        prob = ODEProblem(soilmoisture,s0,tspan,p)
        sol = solve(prob, reltol=1e-6,saveat=1.0)

        #Plot s(t) over t
        lines!(sol.t,sol.u,label= "Pa = $(pa[i])")

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

elseif model_version == 2

    ### VERSION 2 ###
    #with parametrised precipitation and s, wl, wo as state variables of the system

    # Defining the system of equations

    function f_trans!(F, x, p)

        @unpack spwp, sfc, ep, eo, ϵ, r, α, nZr, a, b, wsat, u, L = p

        #Variable naming: x[1] = s, x[2] = wl, x[3] = wo
        F[1] = 1/(nZr) * (precip(x[2], wsat, a, b) * infiltration(x[1], ϵ, r) - land_evap_trans(x[1], spwp, sfc, ep))
        F[2] = land_evap_trans(x[1], spwp, sfc, ep) - precip(x[2], wsat, a, b) + (x[3] - x[2]) * u / (α*L)
        F[3] = eo - precip(x[3], wsat, a, b) - (x[3] - x[2]) * u / ((1-α) * L)
    end


    function cm_monte_carlo_run(nb_runs)

        param_names = keys(cm_rand_params())
        col_names = [string(el) for el in param_names]
        typeof(col_names)
        col_names = append!(col_names, ["s", "wl", "wo", "convergence"])
        output_columns = length(col_names)

        #initialising output matrix
        sol = Array{Float64}(undef, 0, output_columns)
    
        for i=1:nb_runs
            init_cond = [0.5, 50.0, 50.0]
            sol = [sol; cm_eq_solution(f_trans!, init_cond)]
        end
    
        sol_df = DataFrame(sol, col_names)
        CSV.write(datadir("sims", "cm_eq_MonteCarlo_scan_$(nb_runs)_runs.csv"), sol_df)
        println(sol_df)
    end

    cm_monte_carlo_run(10)



elseif model_version == 3

    function closed_model(dx, x, p, t)
        @unpack spwp, sfc, ep, ptot, eo, ϵ, r, α, nZr, a, b, wsat, u, L = p
        dx[1] = 1/(nZr) * (precip(x[2], wsat, a, b) * infiltration(x[1], ϵ, r) - land_evap(x[1], spwp, sfc, ep))
        dx[2] = land_evap(x[1], spwp, sfc, ep) - precip(x[2], wsat, a, b) + (x[3] - x[2]) * u / (α*L)
        dx[3] = eo - precip(x[3], wsat, a, b) - (x[3] - x[2]) * u / ((1-α) * L)
    end

    tspan = (0.0, 100.0)
    ic = [0.6, 30, 60]
    prob = ODEProblem(closed_model, ic, tspan, p_fixed)
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


elseif model_version == 4

    dynsys = ContinuousDynamicalSystem(closed_model_piecewise, [0.1, 10, 60], p_fixed)

    s_range = 0.0..1.0
    wl_range = 0.0..wsat 
    wo_range = 0.0..wsat
    box = s_range × wl_range × wo_range
    fp, eigs, stable = fixedpoints(dynsys, box)

elseif model_version == 5

    dynsys = ContinuousDynamicalSystem(closed_model_smooth, [0.1, 10, 60], p_fixed)

    s_range = 0.0..1.0
    wl_range = 0.0..wsat 
    wo_range = 0.0..wsat
    box = s_range × wl_range × wo_range
    fp, eigs, stable = fixedpoints(dynsys, box)

    
elseif model_version == 6

    closure_pw = x -> closed_model_piecewise(x, p_fixed, 0.0)
    s_range = 0.0..1.0
    wl_range = 0.0..wsat 
    wo_range = 0.0..wsat
    box = s_range × wl_range × wo_range
    rts = IntervalRootFinding.roots(closure_pw, box)
    println(rts)

else
    println("Specify which version of model formulation you want to work with")
end