using DrWatson
@quickactivate "Oceland Model"
include(srcdir("parametrisations.jl"))
using PyPlot

fs = 14

s   = collect(0.0:0.01:1.0)     # relative soil moisture saturation, unitless
r   = [1.5,2.0,2.5] # value in [1]: c=1 - linear model
ϵ   = [0.5,1.0,1.5]

infilt = zeros(Float64,length(s),length(r)*length(ϵ))

pygui(true)
leg = Array{String,1}[]
combi = 1 #keeps track of combinations of parameters (used to access columns of infilt)
for m = 1:length(ϵ)
    for i = 1:length(r)    
        for n = 1:length(s)
            infilt[n,combi] = infiltration(ϵ[m],r[i],s[n])
        end
        local p    = plot(s,infilt[:,combi])
        leg_text = string("ϵ=",ϵ[m]," r=",r[i])
        global leg  = [leg;leg_text]
        combi = combi+1
    end
end

legend(leg)
xlabel("relative soil saturation \$s\$",fontsize=fs)
ylabel("infiltration \$infiltration(s)\$",fontsize=fs)
#title("Evapotranspiration Parametrization",fontsize=fs)
text(0.7,0.95,"\$infiltration(s)=1 - ϵ*s^r\$",fontsize=fs)
savefig(plotsdir("runoff_parametrisation.pdf"))