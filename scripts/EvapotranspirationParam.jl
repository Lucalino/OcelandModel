using DrWatson
@quickactivate "Oceland Model"
using PyPlot

fs = 14

s   = collect(0.0:0.01:1.0)     # relative soil moisture saturation, unitless
E_p = 4.38                      # value from [1] in mm/day
c   = [0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0] # value in [1]: c=1 - linear model
leg = Array{String,1}[]

for m = 1:length(c)
    EE_p= E_p .* s.^c[m] ./ E_p
    local p    = plot(s,EE_p)
    text = string("c=",c[m])
    global leg  = [leg;text]
end

legend(leg)
xlabel("relative soil saturation \$s\$",fontsize=fs)
ylabel("\$E(s)\$/\$E_p\$",fontsize=fs)
#title("Evapotranspiration Parametrization",fontsize=fs)
text(0.0,0.95,"\$E(s)=E_p * s^c\$",fontsize=fs)