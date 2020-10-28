# OCELAND TOY MODEL 

# written by Luca Schmidt
# started 29th September 2020

# The oceland model comprises two domains, an ocean and a land part, 
# representing a tropical island.

# The aim of the model is to constrain the ratio of precipitation over 
# land and ocean, P_l/P_o, based on simple equations of water vapor 
# and soil moisture conservation. 

using DrWatson
@quickactivate "Oceland Model"
using PyPlot

const D   = 100            #spatial extent of entire domain
const L   = 10             #spatial extent of island
const O   = D-L            #spatial extent of ocean 
const E_p = 100            #potential evapotranspiration
const c   = 1              #positive const
const r   = 1              #positive const
const ϵ   = 1              #positive const
const D   = 1              #linear scale of horizontal air flow
s   = collect(0.0:0.2:1.0) #soil moisture
# P_a = collect(0.0:20.0:E_p) #advected water vapor, ranges between zero and E_p

# for n = 1:length(P_a)
#     for i = 1:length(s)
#         P_o = D/2*E_p*ϵ*s[i]^(c+r) - (D*E_p-2E_p)*s[i]^c + P_a[n]*ϵ*s[i]^r + E_p - 2P_a[n];
#         P_o_exp = E_p - P_a[n];
#         println("s=",s[i])
#         println("P_o=",P_o)
#         println("Expected P_o=",P_o_exp)
#     end
# end





# E_l(s) = E_p * s^c         #evapotranspiration from land
# Φ(s)   = 1 - ϵ*s^r         #infiltration function (runoff)
# P_l(s) = P_a * (1+s^c/Ω)

# for n in s
#     P_o = E_o + 2E_l(s[n]) - (Φ(s[n])+1)P_l(s[n])
# end

