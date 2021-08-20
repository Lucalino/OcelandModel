using DrWatson
@quickactivate "Oceland Model"
include(srcdir("parametrisations.jl"))

#using PyPlot
#using Plots

#General definitions
lw = 3.0 #linewidth
fs = 18.0 #fontsize
#pygui(true)

# Called functions
# land_evapo(spwp, sfc, Ep, s)
# infiltration(ϵ,r,s)

# In this script I want to investigate the sign change of dPR/d alpha.
# I will do so by plotting Pl and Po below and above the critical soil moisture value
# where the sign chang occurs. 
# Because I don't find the right key on Steffis keyboard, alpha is called a in this script.
# Can be replaced by alpha later on.

 spwp = 0.1 # NEEDS SOME MORE RESEARCH
 sfc = 0.8
 Ep = 4.38 #potential evaporation over land in mm/day, taken from [1]
 eo = 3.0
 ϵ = 1.0
 Ptot = 3.0
 r = 2.0
#a = collect(Float64,0.0:0.1:1.0)
# ss = 0.3
# sl = 0.6
#s = collect(Float64,0.0:0.1:1.0) #soil moisture
# #PR = zeros(length(s)) #precipitation ratio
# #Els  = zeros(length(a)) #land evapotranspiration rate
# #Ell  = zeros(length(a))
# #Eo  = zeros(length(a)) #ocean evaporation rate
# #inflts  = zeros(length(a))
# #infltl  = zeros(length(a))
# #infilt  = zeros(length(s)) #infiltration rate
# #Pos      = zeros(length(a)) #ocean precipitation rate
# #Pol      = zeros(length(a)) 
# #Pls      = zeros(length(a)) #land precipitation rate
# #Pll      = zeros(length(a)) 
# El       = zeros(Float64, length(s))
# infilt   = zeros(Float64, length(s))
# Pl       = zeros(Float64, length(s), length(a))
# Po       = zeros(Float64, length(s), length(a))

# figure()
# for i = 1:length(s)
#     infilt[i] = infiltration(ϵ,r,s[i])
#     El[i] = land_evapo(spwp, sfc, Ep, s[i])
#     for n = 1:length(a)
#     #Els[i] = land_evapo(spwp, sfc, Ep, ss)
#     #Ell[i] = land_evapo(spwp, sfc, Ep, sl)
#     #inflts[i] = infiltration(ϵ,r,ss)
#     #infltl[i] = infiltration(ϵ,r,sl)
#         Pl[i,n] = El[i]/infilt[i]
#     #Pls[i] = Els[i]/inflts[i]
#         Po[i,n] = Pl[i,n] * (Ptot * infilt[i] - a[n] * El[i])/(El[i]*(1-a[n]))
#     #Pos[i] = Pls[i] * (Ptot * inflts[i] - a[i] * Els[i])/(Els[i]*(1-a[i]))

#     #Pll[i] = Ell[i]/infltl[i]
#     #Pol[i] = Pll[i] * (Ptot * infltl[i] - a[i] * Ell[i])/(Ell[i]*(1-a[i]))
#     end
#     plot(a,Po[i,:], label=string("s=",s[i]))
# end
# axhline(0.0; linestyle="dotted", color= "black")
# #legend()
# xlabel("Land fraction α")
# ylabel("Ocean precipitation [mm/day]")
# legend()

# figure()
# plot(a,Pls,linewidth=lw, label="Pl_s")
# plot(a,Pll,linewidth=lw, label="Pl_l")
# plot(a,Pos,linewidth=lw, label="Po_s")
# plot(a,Pol,linewidth=lw, label="Po_l")
# xlabel("Land fraction α",fontsize=fs)
# ylabel("Precipitation [mm/day]",fontsize=fs)
# legend()

# for n = 1:length(s)
#     for m=1:length(a)
#         Pl[n,m] = land_evapo(spwp, sfc, Ep, s[n])/infiltration(ϵ,r,s[n])
#         Po[n,m] = Pl[n,m] * (Ptot * infiltration(ϵ,r,s[n]) - a[m] * land_evapo(spwp, sfc, Ep, s[n]))/(land_evapo(spwp, sfc, Ep, s[n])*(1-a[m]))
#     end
# end

using Plots; pyplot()
a = range(0.2,stop=0.5,length=50)
s = range(0.1,stop=0.6, length=50)
f(s,a)=3.0
Pl_func(s,a) = land_evapo(spwp, sfc, Ep, s)/infiltration(ϵ,r,s)
Po_func(s,a) = (land_evapo(spwp, sfc, Ep, s)/infiltration(ϵ,r,s)) * (Ptot * infiltration(ϵ,r,s) - a * land_evapo(spwp, sfc, Ep, s))/(land_evapo(spwp, sfc, Ep, s)*(1-a))
plot(s,a,Po_func, st=:surface,c= cgrad(:greens),title = "Po",camera=(10,10))
plot!(s,a,f, st=:surface)
plot!(size(400,800))
#surface(s,a,Po_func,α=0.9,camera=(160,10))
#surface(s,a,Pl_func,α=0.9,camera=(160,10))
xlabel!("s")
ylabel!("α")
#avefig(plotsdir("Sensitivity/Ptot fixed/3D_alpha_0.2-0.5view10_10.pdf"))
# x=[1,2,3,4,5]
# y=[1,2,3,4,5]
# x,y = np.meshgrid(x,y)
# z(x,y)=x+y
# figure()
# plot_surface(x,y,z)