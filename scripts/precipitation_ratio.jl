using DrWatson
@quickactivate "Oceland Model"
include(srcdir("parametrisations.jl"))

using PyPlot

#General definitions
lw = 3.0 #linewidth
fs = 18.0 #fontsize

# Called functions
# land_evap(s, spwp, sfc, Ep)
# infiltration(s,ϵ,r)

constraint = "Ptot"


spwp = 0.2 # NEEDS SOME MORE RESEARCH
sfc = 0.75
Ep = 4.38 #potential evaporation over land in mm/day, taken from [1]
eo = 3.0
#ϵ = collect(0.9:0.1:1.1)
ϵ = 1.0
Ptot = 3.0

#r = collect(1.9:0.1:2.1)
r = 2.0
α = 0.4 #land fraction
s = collect(0.0:0.01:1.0) #soil moisture
PR = zeros(Float64,length(s),length(r)*length(ϵ))
#PR = zeros(length(s)) #precipitation ratio
El  = zeros(length(s)) #land evapotranspiration rate
Eo  = zeros(length(s)) #ocean evaporation rate
infilt  = zeros(length(s)) #infiltration rate
Po      = zeros(length(s)) #ocean precipitation rate
Pl      = zeros(length(s)) #land precipitation rate
Pl_wvpb = zeros(length(s)) #land precipitation rate
Aol     = zeros(length(s)) #advection rate into land domain
R       = zeros(length(s)) #runoff rate from land domain
R_test  = zeros(length(s))
fluxes_init = zeros(length(s)) #initial water influx in model domain
fluxes_der = zeros(length(s)) #fluxes derived from initial influxes

pygui(true)



if constraint == "Eo"

    #Precipitation ratio for combination of parameters
    #= leg = Array{String,1}[]
    combi = 1 #keeps track of combinations of parameters (used to access columns of infilt)
    for m = 1:length(ϵ)
        for n = 1:length(r)    
            for i = 1:length(s)
                PR[i,combi] = land_evap(s[i],spwp,sfc,Ep) * (1-α) / (α * land_evap(s[i],spwp,sfc,Ep) * (infiltration(ϵ[m],r[n],s[i]) - 1) + eo * infiltration(ϵ[m],r[n],s[i]) * (1-α))
            end
            local p    = plot(s,PR[:],combi])
            leg_text = string("ϵ=",ϵ[m]," r=",r[n])
            global leg  = [leg;leg_text]
            combi = combi+1
        end
    end

    legend(leg)
    xlabel("Relative soil saturation \$s\$",fontsize=fs)
    ylabel("Precipitation ratio \$P_l/P_o\$",fontsize=fs)
    #text(0.405,0.5,L"$\frac{P_l}{P_o}(s)=\frac{E_l(1-α)}{αE_l(Φ-1)+E_o Φ (1-α)}$",fontsize=fs) #equation
    #savefig(plotsdir("PR_param_combi_0.5-0.7.pdf")) =#

    #Precipitation ratio for one set of parameters

    for i = 1:length(s)

        El[i] = land_evap(s[i],spwp,sfc,Ep)
        Eo[i] = eo         
        Po[i] = eo + α*land_evap(s[i],spwp,sfc,Ep)/(1-α) * (1-1/infiltration(s[i],ϵ,r))
        Pl[i] = land_evap(s[i],spwp,sfc,Ep)/infiltration(s[i],ϵ,r) #from soil moisture balance
        #Pl_wvpb[i]= (1-α)*(eo-Po[i])/α + land_evap(s[i],spwp,sfc,Ep) #from water vapour balance (equivalent to Pl)
        infilt[i] = Pl[i] * infiltration(s[i],ϵ,r) 
        Aol[i]    = eo - Po[i]
        R[i]      = Pl[i] - El[i]   #equivalent to (Pl[i]-infilt[i]) or (-Pl*(Φ-1)=-Pl*ϵs^r)
        PR[i]     = land_evap(s[i],spwp,sfc,Ep) * (1-α) / (α * land_evap(s[i],spwp,sfc,Ep) * (infiltration(s[i],ϵ,r) - 1) + eo * infiltration(s[i],ϵ,r) * (1-α))
        fluxes_init[i]= (1-α) * eo + α * land_evap(s[i],spwp,sfc,Ep)  #water sources in model 
        fluxes_der[i]= α * (Pl[i] - R[i]) + (1-α) * (Po[i] + Aol[i])  #derived fluxes (including transfer fluxes), should never exceed maxflux
    end

    ll = 1 #lower limit of s range (min 1)
    ul = 81 #upper limit of s range (max 101)
    figure(figsize=(16.0,11.0))
    fig4 = plot(s[ll:ul],(1-α) * Po[ll:ul],label="ocean precipitation",linewidth=lw)
    #fig4 = plot(s[ll:ul],infilt[ll:ul], label="infiltration",linewidth=lw)
    fig4 = plot(s[ll:ul],α * Pl[ll:ul];label="land precipitation",lw)
    fig4 = plot(s[ll:ul],α * El[ll:ul]; label="land evapotranspiration",lw)
    fig4 = plot(s[ll:ul],(1-α) * Eo[ll:ul]; label="ocean evaporation",lw)
    fig4 = plot(s[ll:ul],α * R[ll:ul]; label="runoff",lw)
    fig4 = plot(s[ll:ul],(1-α) * Aol[ll:ul]; label="advection",lw=3.5,linestyle="dotted", color="deeppink")
    fig4 = plot(s[ll:ul],fluxes_init[ll:ul]; label="initial fluxes: (1-α)Eo + αEl",lw, color="darkgreen")
    fig4 = plot(s[ll:ul],fluxes_der[ll:ul]; label="derived fluxes: (1-α)[Po + Aol] + α[Pl - R]",lw=3.5, linestyle="dotted", color="deepskyblue")
    axhline(0.0; linestyle="dotted", color= "black")
    xlabel("Relative soil moisture saturation s",fontsize=fs)
    ylabel("Water fluxes [mm/day]",fontsize=fs) 
    title("Total fluxes (scaled by domain fraction)", fontsize=fs)
    legend(fontsize=14)
    #savefig(plotsdir("Fluxes/Water_fluxes_bettercolors.pdf"))
    

elseif constraint == "Ptot"

    #Precipitation ratio and fluxes for one set of parameters

    for i = 1:length(s)
        El[i] =  land_evap(s[i],spwp,sfc,Ep)
        Eo[i] = (Ptot - α * El[i])/(1-α)  
        Pl[i]  = El[i]/infiltration(s[i],ϵ,r) #from soil moisture balance       
        Po[i] = Pl[i] * (Ptot * infiltration(s[i],ϵ,r) - α * El[i])/(El[i]*(1-α))
        infilt[i] = Pl[i] * infiltration(s[i],ϵ,r) 
        Aol[i]    = Eo[i] - Po[i]
        R[i]      = Pl[i] - El[i]   #equivalent to (Pl[i]-infilt[i]) or (-Pl*(Φ-1)=-Pl*ϵs^r)
        PR[i]     = El[i] * (1-α) / (Ptot * infiltration(s[i],ϵ,r) - α*El[i])
        fluxes_init[i]= (1-α) * Eo[i] + α * El[i]  #water sources in model 
        fluxes_der[i]= α * (Pl[i] - R[i]) + (1-α) * (Po[i] + Aol[i])  #derived fluxes (including transfer fluxes), should never exceed maxflux
    end

    ll = 1 #lower limit of s range (min 1)
    ul = 83 #upper limit of s range (max 101)
    # add_text = string("ϵ=",ϵ," r=",r)

    # figure()
    # plot(s,PR,linewidth=lw)
    # xlabel("Soil moisture saturation s",fontsize=fs)
    # ylabel(L"Precipitation Ratio $P_l/P_o$",fontsize=fs) 
    # text(0.0,100,add_text)
    # savefig(plotsdir("PR/Ptot constrained/PR_fixed_param.pdf"))

    #OLD PLOT WITH INITIAL and DERIVED FLUXES
    # figure(figsize=(16.0,11.0))
    # plot(s[ll:ul],(1-α) * Po[ll:ul],label="ocean precipitation",linewidth=lw)
    # #fig4 = plot(s[ll:ul],infilt[ll:ul], label="infiltration",linewidth=lw)
    # plot(s[ll:ul],α * Pl[ll:ul];label="land precipitation",lw)
    # plot(s[ll:ul],α * El[ll:ul]; label="land evapotranspiration",lw)
    # plot(s[ll:ul],(1-α) * Eo[ll:ul]; label="ocean evaporation",lw)
    # plot(s[ll:ul],α * R[ll:ul]; label="runoff",lw)
    # plot(s[ll:ul],(1-α) * Aol[ll:ul]; label="advection",lw=3.5,linestyle="dotted", color="deeppink")
    # plot(s[ll:ul], α * El[ll:ul] + (1-α) * Eo[ll:ul]; label="αEl + (1-α)Eo", lw, color="darkgreen")
    # plot(s[ll:ul], α * Pl[ll:ul] + (1-α) * Po[ll:ul]; label="αPl + (1-α)Po", lw=3.5, linestyle="dotted", color="deepskyblue")
    
    # #plot(s[ll:ul],fluxes_init[ll:ul]; label="initial fluxes: (1-α)Eo + αEl",lw, color="darkgreen")
    # #plot(s[ll:ul],fluxes_der[ll:ul]; label="derived fluxes: (1-α)[Po + Aol] + α[Pl - R]",lw=3.5, linestyle="dotted", color="deepskyblue")
    # axhline(0.0; linestyle="dotted", color= "black")
    # xlabel("Relative soil moisture saturation s",fontsize=fs)
    # ylabel("Water fluxes [mm/day]",fontsize=fs) 
    # #title(string("Total fluxes, α=",α), fontsize=fs)
    # legend(fontsize=12)
    # #f.savefig(plotsdir("Fluxes/Ptot constrained/Water_fluxes_panel_report.pdf"))

    #NEW PLOT WITH SUM OF EVAPORATION AND PRECIPITATION FLUXES AND UPDATED LAYOUT
    fig, ax = plt.subplots(figsize=(12.0,8.0),constrained_layout=true)
    ax.plot(s[ll:ul],(1-α) * Po[ll:ul],label="ocean precipitation",linewidth=lw)
    ax.plot(s[ll:ul],α * Pl[ll:ul];label="land precipitation",lw)
    ax.plot(s[ll:ul],α * El[ll:ul]; label="land evapotranspiration",lw)
    ax.plot(s[ll:ul],(1-α) * Eo[ll:ul]; label="ocean evaporation",lw)
    ax.plot(s[ll:ul],α * R[ll:ul]; label="runoff",lw)
    ax.plot(s[ll:ul],(1-α) * Aol[ll:ul]; label="advection",lw=3.5,linestyle="dotted", color="deeppink")
    ax.plot(s[ll:ul], α * El[ll:ul] + (1-α) * Eo[ll:ul]; label="αEl + (1-α)Eo", lw, color="darkgreen")
    ax.plot(s[ll:ul], α * Pl[ll:ul] + (1-α) * Po[ll:ul]; label="αPl + (1-α)Po", lw=3.5, linestyle="dotted", color="deepskyblue")
    axhline(0.0; linestyle="dotted", color= "black")
    ax.set_xlabel("Relative soil moisture saturation s",fontsize=fs)
    ax.set_ylabel("Water fluxes [mm/day]",fontsize=fs) 
    ax.legend(fontsize=16, frameon=false, loc="upper left")
    ax.spines["right"].set_visible(false)
    ax.spines["top"].set_visible(false)
    ax.set_ylim([-1.5,6.0])
    ax.set_xlim([0.0,0.82])
    ax.tick_params(labelsize=14)
    #savefig(plotsdir("Closed Model/Fluxes/Ptot constrained/Water_fluxes_panel_report.png"))


    # #just rates
    # figure(figsize=(16.0,11.0))
    # plot(s[ll:ul],Po[ll:ul],label="ocean precipitation",linewidth=lw)
    # #fig4 = plot(s[ll:ul],infilt[ll:ul], label="infiltration",linewidth=lw)
    # plot(s[ll:ul],Pl[ll:ul];label="land precipitation",lw)
    # plot(s[ll:ul],El[ll:ul]; label="land evapotranspiration",lw)
    # plot(s[ll:ul],Eo[ll:ul]; label="ocean evaporation",lw)
    # plot(s[ll:ul],R[ll:ul]; label="runoff",lw)
    # plot(s[ll:ul],Aol[ll:ul]; label="advection",lw=3.5,linestyle="dotted", color="deeppink")
    # #plot(s[ll:ul],fluxes_init[ll:ul]; label="initial fluxes: (1-α)Eo + αEl",lw, color="darkgreen")
    # #plot(s[ll:ul],fluxes_der[ll:ul]; label="derived fluxes: (1-α)[Po + Aol] + α[Pl - R]",lw=3.5, linestyle="dotted", color="deepskyblue")
    # axhline(0.0; linestyle="dotted", color= "black")
    # xlabel("Relative soil moisture saturation s",fontsize=fs)
    # ylabel("Water flux rates [mm/day]",fontsize=fs) 
    # title(string("Flux rates, α=",α), fontsize=fs)
    # legend(fontsize=10)

    #Precipitation ratio for range of parameters

    # leg = Array{String,1}[]
    # combi = 1 #keeps track of combinations of parameters (used to access columns of infilt)
    # figure(figsize=(9,6))
    # for m = 1:length(ϵ)
    #     for n = 1:length(r)    
    #         for i = 1:length(s)
    #             PR[i,combi] = land_evap(s[i],spwp,sfc,Ep) * (1-α) / (Ptot * infiltration(s[i],ϵ[m],r[n]) - α * land_evap(s[i], spwp, sfc, Ep))
    #         end
    #         local p    = plot(s[ll:ul],PR[ll:ul,combi])
    #         leg_text = string("ϵ=",ϵ[m]," r=",r[n])
    #         global leg  = [leg;leg_text]
    #         combi = combi+1
    #     end
    # end
    
    # legend(leg)
    # xlabel("Relative soil saturation \$s\$",fontsize=fs)
    # ylabel("Precipitation ratio \$P_l/P_o\$",fontsize=fs)
    # text(0.2,2.5,L"$\frac{P_l}{P_o}(s)=\frac{E_l(1-α)}{P_{tot} \Phi - \alpha E_l}$",fontsize=fs) #equation
    # #savefig(plotsdir("PR/Ptot constrained/PR_param_combi_0.0-0.6.pdf"))
end

#Precipitation ratio
#= add_text = string("ϵ=",ϵ," r=",r)
figure()
fig3 = plot(s,PR,linewidth=lw)
xlabel("Soil moisture saturation s",fontsize=fs)
ylabel(L"Precipitation Ratio $P_l/P_o$",fontsize=fs) 
text(0.0,47,add_text)
#savefig(plotsdir("PR.pdf")) =#


#= #Evapotranspiration
figure()
fig1 = plot(s,evap,linewidth=lw)
xlabel("Soil moisture saturation s",fontsize=fs)
ylabel("Evapotranspiration [mm/day]",fontsize=fs)
grid("on")
savefig(plotsdir("evap.pdf")) =#


#= #Infiltration function
figure()
fig2 = plot(s,infilt,linewidth=lw)
xlabel("soil moisture staturation s", fontsize=fs)
ylabel("infiltration function Φ", fontsize=fs)
text(0.7, 0.9, add_text, fontsize=fs) #PROBLEM

savefig(plotsdir("infilt.pdf")) =#

