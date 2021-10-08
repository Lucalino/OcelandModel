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
using PyPlot
pygui(true)
using DifferentialEquations
using NLsolve

fs=14.0
lw=2.0

# Model versions:
# 1  : soil moisture ODE without parametrised precipitation but given advection Pa
# 2  : set of equilibrium equations for the three state variables s, wl, wo. Precipitaiton is parametrised
# 3  : set of ODEs for the three state variables s, wl, wo. Precipitation is parametrised

model_version = 2

if model_version == 1

    spwp = 0.1 #permanent wilting point, NEEDS SOME MORE RESEARCH!
    sfc = 0.8  #field capacity
    Ep = 4.38  #potential evaporation over land in mm/day, taken from [1]
    ptot = 3.0 #average total precipitation rate over the full tropics
    eo = 3.0   #ocean evaporation rate
    pa = [0.5,0.7,0.9]   #advection rate into land domain (which entirely falls as rain due to closed system)
    ϵ = 1.0
    r = 2.0
    α = 0.3    #land fraction
    nZr = 100  #reservoir depth/"field storage capacity of the soil" [mm]


    function soilmoisture(s,p,t)
        spwp,sfc,Ep,ptot,eo,pa,ϵ,r,α,nZr = p
        #ds =( ( (1-α) * pa + land_evapo(spwp, sfc, Ep, s) ) * infiltration(ϵ,r,s) - land_evapo(spwp, sfc, Ep, s) ) / nZr
        ds = (pa * infiltration(ϵ,r,s) + land_evapo(spwp,sfc,Ep,s) * (infiltration(ϵ,r,s) - 1)) / nZr
    end


    figure()
    for i = 1:3

        #Solve the ODE for s(t)
        p = (spwp,sfc,Ep,ptot,eo,pa[i],ϵ,r,α,nZr)
        s0 = 0.3 #initial condition
        tspan = (0.0,150.0) #in days
        prob = ODEProblem(soilmoisture,s0,tspan,p)
        sol = solve(prob,reltol=1e-6,saveat=1.0)


        #Plot s(t) over t
        plot(sol.t,sol.u,label=string("Pa =",pa[i]))

    end
    xlabel("time \$t\$ [days]",fontsize=fs)
    ylabel("soil moisture saturation \$s\$",fontsize=fs)
    legend()
    #savefig(plotsdir(string("Closed model/Dynamics/sm_Pa0.5-0.9.pdf")))


    #Compute Pl(t), Po(t) and PR(t) using the computed s(t)
    # Pl = zeros(length(sol.t))
    # Po = zeros(length(sol.t))
    # PR = zeros(length(sol.t))

    # for t = 1:length(sol.u)
    #     Pl[t] = land_evapo(spwp,sfc,Ep,sol.u[t]) + pa
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

    spwp = 0.35 #permanent wilting point, NEEDS SOME MORE RESEARCH!
    sfc = spwp + 0.3  #field capacity
    ep = 4.38  #[mm/day] potential evaporation over land in mm/day, taken from [1]
    ptot = 3.0 #[mm/day] average total precipitation rate over the full tropics
    eo = 3.0   #[mm/day] ocean evaporation rate
    ϵ = 1.0    #numerical parameter from Rodriguez-Iturbe et al. (1991)
    r = 2.0    #numerical parameter from Rodriguez-Iturbe et al. (1991)
    α = 0.2    #land fraction
    nZr = 100.0  #[mm] reservoir depth/"field storage capacity of the soil" [mm]
    a = 15.6   #numerical parameter from Bretherton et al. (2004)
    b = 0.603  #numerical parameter from Bretherton et al. (2004)
    wsat = 72.0  #[mm] saturation water vapour pass derived from plots in Bretherton et al. ( 2004)
    u = 5.0 * m2mm(1)/s2day(1) #[mm/day] wind speed
    L = 1000 * km2mm(1) #[mm] domain size


    # Defining the system of equations

    function f_trans!(F, x)
        #Variable naming: x[1] = s, x[2] = wl, x[3] = wo
        F[1] = 1/(nZr) * (precip(x[2], wsat, a, b) * infiltration(x[1], ϵ, r) - land_evap_trans(x[1], spwp, sfc, ep))
        F[2] = land_evap_trans(x[1], spwp, sfc, ep) - precip(x[2], wsat, a, b) + (x[3] - x[2]) * u / (α*L)
        F[3] = eo - precip(x[3], wsat, a, b) - (x[3] - x[2]) * u / ((1-α) * L)
    end

    sol = nlsolve(f_trans!, [0.9, 100, 100])


    function cm_monte_carlo_run(nb_runs)

        output_columns = 16
        col_names = ["spwp", "sfc", "ep", "eo", "ϵ", "r", "α", "nZr", 
                     "a", "b", "w_sat", "u", "L", "s", "wl", "wo", "converged if 1"]
    
        #initialising output matrix
        sol = Array{Float64}(undef, 0, output_columns)
    
        for i=1:nb_runs
            params = cm_rand_params(w0_fixed)
            sol = [sol; om_eq_solution(f_one, f_two, f_three, params, output_columns)]
        end
    
        sol_df = DataFrame(sol, col_names)
        CSV.write(datadir("sims", "om_eq_MonteCarlo_scan_$(nb_runs)_runs_domain10000.csv"), sol_df)
    end







elseif model_version == 3

    #parameters
    spwp = 0.35 #permanent wilting point
    sfc = spwp + 0.3  #field capacity
    ep = 4.38  #[mm/day] potential evaporation over land in mm/day, taken from [1]
    ptot = 3.0 #[mm/day] average total precipitation rate over the full tropics
    eo = 3.0   #[mm/day] ocean evaporation rate
    ϵ = 1.0    #numerical parameter from Rodriguez-Iturbe et al. (1991)
    r = 2.0    #numerical parameter from Rodriguez-Iturbe et al. (1991)
    α = 0.2    #land fraction
    nZr = 100.0  #[mm] reservoir depth/"field storage capacity of the soil" [mm]
    a = 15.6   #numerical parameter from Bretherton et al. (2004)
    b = 0.603  #numerical parameter from Bretherton et al. (2004)
    wsat = 72.0  #[mm] saturation water vapour pass derived from plots in Bretherton et al. ( 2004)
    u = 5.0 * m2mm(1)/s2day(1) #[mm/day] wind speed
    L = 1000 * km2mm(1) #[mm] domain size
    p = @dict spwp, sfc, ep, ptot, eo, ϵ, r, α, nZr, a, b, wsat, u, L


    function closed_model(x, p, t)
        @unpack spwp, sfc, ep, ptot, eo, ϵ, r, α, nZr, a, b, wsat, u, L = p
        ds = 1/(nZr) * (precip(x[2], wsat, a, b) * infiltration(x[1], ϵ, r) - land_evap_trans(x[1], spwp, sfc, ep))
        dwl = land_evap_trans(x[1], spwp, sfc, ep) - precip(x[2], wsat, a, b) + (x[3] - x[2]) * u / (α*L)
        dwo = eo - precip(x[3], wsat, a, b) - (x[3] - x[2]) * u / ((1-α) * L)
        return SVector{3}[ds, dwl, dwo]
    end
else
    println("Specify which version of model formulation you want to work with")
end