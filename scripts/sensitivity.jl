using DrWatson
@quickactivate "Oceland Model"
include(srcdir("parametrisations.jl"))

using PyPlot
pygui(true)

# Called functions
# land_evap(s, spwp, sfc, Ep)
# infiltration(s,ϵ,r)

#General definitions
lw = 3.0 #linewidth
fs = 16.0 #fontsize

#Computation mode
level = 2
ll = 1 #lower limit of s range (min 1)
ul = 56 #upper limit of s range (max 101)
constraint = "Ptot" # or "Eo"
calcmode = "normal" # or "abs"


spwp = 0.2 #permanent wilting point, NEEDS SOME MORE RESEARCH!
sfc = 0.75  #field capacity
Ep = 4.38  #potential evaporation over land in mm/day, taken from [1]
ptot = 3.0 #average total precipitation rate over the full tropics
eo = 3.0   #fixed ocean evaporation
#ϵ = collect(0.9:0.1:1.1)
ϵ = 1.0
#r = collect(1.9:0.1:2.1)
r = 2.0
α = 0.4 #land fraction
s = collect(0.0:0.01:1.0) #soil moisture



#Initialisation of quantities
El      = zeros(length(s)) #land evapotranspiration rate
Eo      = zeros(length(s)) #ocean evaporation rate
Ptot    = zeros(length(s))
infilt  = zeros(length(s)) #infiltration rate
Po      = zeros(length(s))
Pl_sb   = zeros(length(s))
Aol     = zeros(length(s)) #advection rate into land domain
R       = zeros(length(s)) #runoff rate from land domain
PR      = zeros(length(s)) #precipitation ratio
#PR     = zeros(Float64,length(s),length(r)*length(ϵ))
dPR_El  = zeros(length(s))
dPR_Eo  = zeros(length(s))
dPR_Ptot= zeros(length(s))
dPR_Φ   = zeros(length(s))
dPR_α  = zeros(length(s))
dPR_Ep  = zeros(length(s))
dPR_spwp= zeros(length(s)) 
dPR_sfc = zeros(length(s))
dPR_s   = zeros(length(s))
dPR_ϵ   = zeros(length(s))
dPR_r   = zeros(length(s))


#Calculation of relevant quantities for derivatives
for i = 1:length(s)
    El[i]     = land_evap(s[i],spwp,sfc,Ep)
    Eo[i]     = eo         
    Ptot[i]   = ptot
    Po[i]     = eo + α*El[i]/(1-α) * (1-1/infiltration(s[i],ϵ,r))
    Pl_sb[i]  = El[i]/infiltration(s[i],ϵ,r) #from soil moisture balance
    infilt[i] = Pl_sb[i] * infiltration(s[i],ϵ,r) 
    Aol[i]    = eo - Po[i]
    R[i]      = Pl_sb[i] - El[i]   #equivalent to (Pl_sb[i]-infilt[i])*α)
end

if constraint == "Eo"
    #Level 1: Calculation of precipitation ratio and its partial derivatives wrt El, Eo, Φ and α
    denom1 = zeros(length(s)) #denominater for level 1 analysis
    for i = 1:length(s)
        PR[i]     = El[i] * (1-α) / (α * El[i] * (infiltration(s[i],ϵ,r) - 1) + Eo[i] * infiltration(s[i],ϵ,r) * (1-α))
        denom1[i] = (α * El[i] * (infiltration(s[i],ϵ,r)-1) + Eo[i] * (1-α) * infiltration(s[i],ϵ,r))^2
        dPR_El[i] = (α -1)^2 * Eo[i] * infiltration(s[i],ϵ,r) / denom1[i]
        dPR_Eo[i] = -(α -1)^2 * El[i] * infiltration(s[i],ϵ,r) / denom1[i]
        dPR_Φ[i]  = (α -1) * El[i] * (Eo[i] * (1-α) + El[i] * α) / denom1[i]
        dPR_α[i] = El[i]^2 * (1-infiltration(s[i],ϵ,r)) / denom1[i]
    end

    if level==1

        #Level 1 Plot
        figure(figsize=(8.0,5.0))
        #plot(s[ll:ul],PR[ll:ul]; label="precipitation ratio",lw)
        plot(s[ll:ul],dPR_El[ll:ul]; label=L"|\partial PR/\partial E_l|",lw)
        plot(s[ll:ul],dPR_Eo[ll:ul]; label=L"|\partial PR/\partial E_o|",lw)
        plot(s[ll:ul],dPR_Φ[ll:ul]; label=L"|\partial PR/\partial Φ|",lw)
        plot(s[ll:ul],dPR_α[ll:ul]; label=L"|\partial PR/\partial α|",lw)
        axhline(0.0; linestyle="dotted", color= "grey")
        xlabel("relative soil moisture saturation \$s\$",fontsize=fs)
        ylabel("abs. partial derivative of \$|P_l/P_o|\$",fontsize=fs)
        legend(fontsize=fs)
        #savefig(plotsdir("Sensitivity/Sens_lev1_0.0-0.7_abs.pdf"))

    elseif level ==2
        
        #Level 2: Calculation of additional terms based on the derivatives of PR wrt El and Φ from level 1
        for i = 1:length(s)

            if s[i] <= spwp
                dPR_Ep[i]   = 0.0
                dPR_sfc[i]  = 0.0
                dPR_spwp[i] = 0.0
                dPR_s[i]    = 0.0
            elseif s[i] <= sfc
                dPR_Ep[i]   = dPR_El[i]*(s[i]-spwp)/(sfc-spwp)
                dPR_sfc[i]  = dPR_El[i]*Ep*(spwp-s[i])/(sfc-spwp)^2
                dPR_spwp[i] = dPR_El[i]*Ep*(s[i]-sfc)/(sfc-spwp)^2
                dPR_s[i]    = dPR_El[i]*Ep/(sfc-spwp)+dPR_Φ[i]*(-1)ϵ*r*s[i]^(r-1)
            else
                dPR_Ep[i]   = dPR_El[i]
                dPR_sfc[i]  = 0.0 
                dPR_spwp[i] = 0.0
                dPR_s[i]    = dPR_Φ[i]*(-1)ϵ*r*s[i]^(r-1)
            end

            dPR_ϵ[i] = dPR_Φ[i]*(-1)*s[i]^r
            dPR_r[i] = dPR_Φ[i]*(-1)*ϵ*s[i]^r*log(s[i])
        end

        #Level 2 Plot
        f=figure(figsize=(12.0,7.0))
        #plot(s[ll:ul],PR[ll:ul]; label="precipitation ratio",lw)
        plot(s[ll:ul],dPR_spwp[ll:ul]; label=L"\partial PR/\partial s_{pwp}",lw)
        plot(s[ll:ul],dPR_Ep[ll:ul]; label=L"\partial PR/\partial E_p",lw)
        plot(s[ll:ul],dPR_r[ll:ul]; label=L"\partial PR/\partial r",lw)
        plot(s[ll:ul],dPR_ϵ[ll:ul]; label=L"\partial PR/\partial ϵ",lw)
        plot(s[ll:ul],dPR_Eo[ll:ul]; label=L"\partial PR/\partial E_o",lw)
        plot(s[ll:ul],dPR_s[ll:ul]; label=L"\partial PR/\partial s",lw)
        plot(s[ll:ul],dPR_sfc[ll:ul]; label=L"\partial PR/\partial s_{fc}",lw)
        plot(s[ll:ul],dPR_α[ll:ul]; label=L"\partial PR/\partial α",lw)
        axhline(0.0; linestyle="dotted", color= "grey")
        xlabel("relative soil moisture saturation \$s\$",fontsize=fs)
        ylabel("partial derivative of \$|P_l/P_o|\$",fontsize=fs)
        legend(fontsize=fs)
        #ax = f.get_axes()[1]
        #ax.spines["top"].set_visible(false)
        #ax.spines["right"].set_visible(false)
        #savefig(plotsdir("Sensitivity/Sens_lev2_0.0-1.0.pdf"))
    end
elseif constraint == "Ptot"
    
    #Level 1: Calculation of the partial derivatives of PR wrt El, Ptot, Φ and α
    denom1 = zeros(length(s)) #denominater for level 1 analysis
    for i = 1:length(s)
        denom1[i] = (El[i]*α-Ptot[i]*infiltration(s[i],ϵ,r))^2
        dPR_El[i] = Ptot[i] * infiltration(s[i],ϵ,r) * (1-α) / denom1[i]
        dPR_Ptot[i] = -El[i] * infiltration(s[i],ϵ,r) *(1-α) / denom1[i]
        dPR_Φ[i]  = -El[i] * Ptot[i] * (1-α) / denom1[i]
        dPR_α[i] = El[i] * (El[i] - Ptot[i] * infiltration(s[i],ϵ,r)) / denom1[i]
    end

    if level==1

        #Level 1 Plot
        figure(figsize=(8.0,5.0))
        plot(s[ll:ul],dPR_El[ll:ul]; label=L"\partial PR/\partial E_l",lw)
        plot(s[ll:ul],dPR_Ptot[ll:ul]; label=L"\partial PR/\partial P_{tot}",lw)
        plot(s[ll:ul],dPR_Φ[ll:ul]; label=L"\partial PR/\partial Φ",lw)
        plot(s[ll:ul],dPR_α[ll:ul]; label=L"\partial PR/\partial α",lw)
        axhline(0.0; linestyle="dotted", color= "grey")
        xlabel("relative soil moisture saturation \$s\$",fontsize=fs)
        ylabel("abs. partial derivative of \$P_l/P_o\$",fontsize=fs)
        legend(fontsize=fs)
        #savefig(plotsdir("Sensitivity/Sens_lev1_Ptot_0.0-0.7_abs.pdf"))

    elseif level==2

        if calcmode == "normal"

            #Level 2: Calculation of additional terms based on the derivatives of PR wrt El and Φ from level 1
            for i = 1:length(s)

                if s[i] <= spwp
                    dPR_Ep[i]   = 0.0
                    dPR_sfc[i]  = 0.0
                    dPR_spwp[i] = 0.0
                    dPR_s[i]    = 0.0
                elseif s[i] <= sfc
                    dPR_Ep[i]   = dPR_El[i]*(s[i]-spwp)/(sfc-spwp)
                    dPR_sfc[i]  = dPR_El[i]*Ep*(spwp-s[i])/(sfc-spwp)^2
                    dPR_spwp[i] = dPR_El[i]*Ep*(s[i]-sfc)/(sfc-spwp)^2
                    dPR_s[i]    = dPR_El[i]*Ep/(sfc-spwp)+dPR_Φ[i]*(-1)ϵ*r*s[i]^(r-1)
                else
                    dPR_Ep[i]   = dPR_El[i]
                    dPR_sfc[i]  = 0.0 
                    dPR_spwp[i] = 0.0
                    dPR_s[i]    = dPR_Φ[i]*(-1)ϵ*r*s[i]^(r-1)
                end

                dPR_ϵ[i] = dPR_Φ[i]*(-1)*s[i]^r
                dPR_r[i] = dPR_Φ[i]*(-1)*ϵ*s[i]^r*log(s[i])
            end

            #Level 2 plot

            # #OLD PLOT
            # figure(figsize=(12.0,7.0))
            # #plot(s[ll:ul],PR[ll:ul]; label="precipitation ratio",lw)
            # plot(s[ll:ul],dPR_spwp[ll:ul]; label=L"\partial PR/\partial s_{pwp}",lw)
            # plot(s[ll:ul],dPR_Ep[ll:ul]; label=L"\partial PR/\partial E_p",lw)
            # plot(s[ll:ul],dPR_r[ll:ul]; label=L"\partial PR/\partial r",lw)
            # plot(s[ll:ul],dPR_ϵ[ll:ul]; label=L"\partial PR/\partial ϵ",lw)
            # plot(s[ll:ul],dPR_Ptot[ll:ul]; label=L"\partial PR/\partial P_{tot}",lw)
            # plot(s[ll:ul],dPR_s[ll:ul]; label=L"\partial PR/\partial s",lw)
            # plot(s[ll:ul],dPR_sfc[ll:ul]; label=L"\partial PR/\partial s_{fc}",lw)
            # plot(s[ll:ul],dPR_α[ll:ul]; label=L"\partial PR/\partial α",lw)
            # axhline(0.0; linestyle="dotted", color= "grey")
            # xlabel("relative soil moisture saturation \$s\$",fontsize=fs)
            # ylabel("partial derivative of \$P_l/P_o\$",fontsize=fs)
            # legend(fontsize=fs)
            # #ax = f.get_axes()[1]
            # #ax.spines["top"].set_visible(false)
            # #ax.spines["right"].set_visible(false)
            # #savefig(plotsdir("Sensitivity/Ptot fixed/Sens_lev2_Ptot_0.0-0.5.pdf"))


            #NEW PLOT WITH UPDATED LAYOUT
            fig, ax = plt.subplots(figsize=(12.0,8.0),constrained_layout=true)
            ax.plot(s[ll:ul],dPR_spwp[ll:ul]; label=L"\partial PR/\partial s_{pwp}",lw)
            ax.plot(s[ll:ul],dPR_Ep[ll:ul]; label=L"\partial PR/\partial E_p",lw)
            ax.plot(s[ll:ul],dPR_r[ll:ul]; label=L"\partial PR/\partial r",lw)
            ax.plot(s[ll:ul],dPR_ϵ[ll:ul]; label=L"\partial PR/\partial ϵ",lw)
            ax.plot(s[ll:ul],dPR_Ptot[ll:ul]; label=L"\partial PR/\partial P_{tot}",lw)
            ax.plot(s[ll:ul],dPR_s[ll:ul]; label=L"\partial PR/\partial s",lw)
            ax.plot(s[ll:ul],dPR_sfc[ll:ul]; label=L"\partial PR/\partial s_{fc}",lw)
            ax.plot(s[ll:ul],dPR_α[ll:ul]; label=L"\partial PR/\partial α",lw, color="tab:cyan")
            axhline(0.0; linestyle="dotted", color= "grey")
            ax.set_xlabel("relative soil moisture saturation \$s\$",fontsize=fs)
            ax.set_ylabel("partial derivative of \$P_l/P_o\$",fontsize=fs)
            ax.legend(fontsize=fs, frameon=false)
            ax.spines["right"].set_visible(false)
            ax.spines["top"].set_visible(false)
            ax.set_xlim([0.18,0.55])
            ax.set_ylim([-6.0,12.5])
            ax.tick_params(labelsize=14)
            savefig(plotsdir("Closed Model/Sensitivity/Ptot fixed/Sens_lev2_Ptot_0.18-0.55_panel_report.png"))

        elseif calcmode == "abs"
            #Level 2: Calculation of additional terms based on the derivatives of PR wrt El and Φ from level 1
            for i = 1:length(s)

                if s[i] <= spwp
                    dPR_Ep[i]   = 0.0
                    dPR_sfc[i]  = 0.0
                    dPR_spwp[i] = 0.0
                    dPR_s[i]    = 0.0
                elseif s[i] <= sfc
                    dPR_Ep[i]   = abs(dPR_El[i]*(s[i]-spwp)/(sfc-spwp))
                    dPR_sfc[i]  = abs(dPR_El[i]*Ep*(spwp-s[i])/(sfc-spwp)^2)
                    dPR_spwp[i] = abs(dPR_El[i]*Ep*(s[i]-sfc)/(sfc-spwp)^2)
                    dPR_s[i]    = abs(dPR_El[i]*Ep/(sfc-spwp)+dPR_Φ[i]*(-1)ϵ*r*s[i]^(r-1))
                else
                    dPR_Ep[i]   = abs(dPR_El[i])
                    dPR_sfc[i]  = 0.0 
                    dPR_spwp[i] = 0.0
                    dPR_s[i]    = abs(dPR_Φ[i]*(-1)ϵ*r*s[i]^(r-1))
                end

                dPR_ϵ[i] = abs(dPR_Φ[i]*(-1)*s[i]^r)
                dPR_r[i] = abs(dPR_Φ[i]*(-1)*ϵ*s[i]^r*log(s[i]))
                dPR_Ptot[i] = abs(dPR_Ptot[i])
                dPR_α[i] = abs(dPR_α[i])
            end

            #Level 2 plot
            figure(figsize=(12.0,7.0))
            #plot(s[ll:ul],PR[ll:ul]; label="precipitation ratio",lw)
            plot(s[ll:ul],dPR_spwp[ll:ul]; label=L"|\partial PR/\partial s_{pwp}|",lw)
            plot(s[ll:ul],dPR_Ep[ll:ul]; label=L"|\partial PR/\partial E_p|",lw)
            plot(s[ll:ul],dPR_r[ll:ul]; label=L"|\partial PR/\partial r|",lw)
            plot(s[ll:ul],dPR_ϵ[ll:ul]; label=L"|\partial PR/\partial ϵ|",lw)
            plot(s[ll:ul],dPR_Ptot[ll:ul]; label=L"|\partial PR/\partial P_{tot}|",lw)
            plot(s[ll:ul],dPR_s[ll:ul]; label=L"|\partial PR/\partial s|",lw)
            plot(s[ll:ul],dPR_sfc[ll:ul]; label=L"|\partial PR/\partial s_{fc}|",lw)
            plot(s[ll:ul],dPR_α[ll:ul]; label=L"|\partial PR/\partial α|",lw)
            axhline(0.0; linestyle="dotted", color= "grey")
            xlabel("relative soil moisture saturation \$s\$",fontsize=fs)
            ylabel("abs. partial derivative of \$|P_l/P_o|\$",fontsize=fs)
            legend(fontsize=fs)
            #ax = f.get_axes()[1]
            #ax.spines["top"].set_visible(false)
            #ax.spines["right"].set_visible(false)
            #savefig(plotsdir("Sensitivity/Ptot fixed/Sens_lev2_Ptot_0.0-0.5_abs.pdf"))
        end
    end
end




#SUPPLEMENTARY CODE

# #Level 2 direct calculation wrt to all parameters
# denom21  = zeros(length(s)) #denominator for level 2 analysis between spwp and sfc
# denom22 = zeros(length(s))  #denominator for level 2 analysis above sfc
# for i = 1:length(s)
#     denom21[i]  = (Ep * s[i]^r * (s[i]-spwp) * α * ϵ - Eoi] * (sfc-spwp) * (1-α) * infiltration(s[i],ϵ,r))^2
#     denom22[i] = s[i] * (Ep * s[i]^r * α * ϵ + Eo[i] * (α-1) * (1-ϵ*s[i]^r))^2
#     if s[i] <= spwp
#         dPR_Ep[i]   = 0.0
#         dPR_sfc[i]  = 0.0
#         dPR_spwp[i] = 0.0
#         dPR_s[i]    = 0.0
#         dPR_Eo[i]   = 0.0
#         dPR_α[i]   = 0.0
#         dPR_ϵ[i]    = 0.0
#         dPR_r[i]    = 0.0
#     elseif s[i] <= sfc
#         dPR_Ep[i]   = (Eo[i] * (s[i] - spwp) * (sfc-spwp) * (α-1)^2 * infiltration(s[i],ϵ,r))/denom21[i]
#         dPR_sfc[i]  = (Eo[i] * Ep * (spwp - s[i]) * (α-1)^2 * infiltration(s[i],ϵ,r))/denom21[i]
#         dPR_spwp[i] = (Eo[i] * Ep * (s[i]-sfc) * (α-1)^2 * infiltration(s[i],ϵ,r))/denom21[i]
#         dPR_s[i]    = (Ep * (1-α) * (Ep * r * s[i]^r * (s[i]-spwp)^2 * α * ϵ + Eo[i] * (sfc-spwp) * (α-1) * (-s[i]-(r-1) * s[i]^(1+r) * ϵ + r * s[i]^r * spwp * ϵ)))/(s[i]*denom21[i])
#         dPR_Eo[i]   = (Ep * (s[i]-spwp) * (spwp-sfc) * (α-1)^2 * infiltration(s[i],ϵ,r))/denom21[i]
#         dPR_α[i]   = (Ep^2 * s[i]^r * (s[i]-spwp)^2 * ϵ)/denom21[i]
#         dPR_ϵ[i]    = (Ep * s[i]^r * (spwp - s[i]) * (α-1) * (Eo[i] * (spwp - sfc) * (α-1) + Ep * (s[i] - spwp) * α))/denom21[i]
#         dPR_r[i]    = (Ep * s[i]^r * (spwp - s[i]) * (α-1) * (Eo[i] * (spwp - sfc) * (α-1) + Ep * (s[i] - spwp) * α) * ϵ * log(s[i]))/denom21[i]
#     else
#         dPR_Ep[i]   = (Eo[i] * s[i] * (α-1)^2 * (1 - ϵ*s[i]^r))/denom22[i]
#         dPR_sfc[i]  = 0.0
#         dPR_spwp[i] = 0.0
#         dPR_s[i]    = (Ep * s[i]^r * r * ϵ * (α-1) * (Eo[i] *(α-1) - Ep * α))/denom22[i]
#         dPR_Eo[i]   = (Ep * s[i] * (α-1)^2 * (-1) * (1 - ϵ*s[i]^r))/denom22[i]
#         dPR_α[i]   = (Ep * s[i]^r * s[i] * Ep * ϵ)/denom22[i]
#         dPR_ϵ[i]    = (Ep * s[i]^r * s[i] * (α-1) * (Eo[i] * (α-1) - Ep * α))/denom22[i]
#         dPR_r[i]    = (Ep * s[i]^r * s[i] * (α-1) * (Eo[i] * (α-1) - Ep * α) * ϵ * log(s[i]))/denom22[i]
#     end
# end