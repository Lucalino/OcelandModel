using DrWatson
@quickactivate "Oceland Model"
using PyPlot

fs = 14

s   = collect(0.0:0.01:1.0)     # relative soil moisture saturation, unitless
E_p = 4.38                      # value from [1] in mm/day
c   = [0.1, 0.5,1.0]                     # value in [1]: c=1
r   = [2.0]                     # value in [1]: r=2
ϵ   = [1.0]                     # value in [1]: ϵ=1
leg = Array{String,1}[]

for i = 1:length(c)
    for n = 1:length(r)
        for m = 1:length(ϵ)
            local P    = E_p * s.^c[i] ./ (1 .- ϵ[m].*s.^r[n]) # Eqn. (1) from [1] for equilibrium
            local p    = plot(s,P)
            text = string("c=",c[i]," r=",r[n]," ϵ=",ϵ[m])
            global leg  = [leg;text]
        end
    end
end

legend(leg)
xlabel("relative soil moisture saturation \$s\$",fontsize=fs)
ylabel("rain fall rate \$P\$ [mm/day]",fontsize=fs)
title("\$E_p=4.38\$ mm/day")
#grid("on")

# [1] Rodriguez‐Iturbe (1991), DOI: 10.1029/91WR01035

# subplot(211)
# p1 = plot(a,b, linewidth=3)
# subplot(212)
# p2 = plot(a,c,linewidth=1)
# errorbar(a,b,yerr=r)
# xlabel("x")
# ylabel("y")
# title("Title")

#Comments: pygui(true) will display plots in new window.