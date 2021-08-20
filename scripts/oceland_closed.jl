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
using PyPlot
pygui(true)
using DifferentialEquations

fs=14.0
lw=2.0

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