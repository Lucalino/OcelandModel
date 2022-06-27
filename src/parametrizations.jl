"""
    El_tanh(s, p::Dict{Symbol, Float64})

Compute the mean evaporation rate from the land surface as a function of soil moisture, s, and the model parameters 
permanent wilting point, spwp, field capacity, sfc, and potential evaporation, Ep. 
The shape of the function is inspired by a sketch in Seneviratne et al. (2010), DOI:10.1016/j.earscirev.2010.02.004.
"""
function El_tanh(s, p::Dict{Symbol, Float64})
    @unpack spwp, sfc, ep, pt = p
    return ep/2 * tanh( pt * (s - (spwp+sfc)/2 ) ) + ep/2
end




"""
    infiltration(s, p::Dict{Symbol, Float64})

Compute the fraction of land precipitation that infiltrates the soil as a function of soil moisture, s, and the model 
parameters ϵ and r. Defined e.g. in Rodriguez‐Iturbe et al. (1991), DOI: 10.1029/91WR01035.
"""
function infiltration(s, p::Dict{Symbol, Float64})
    @unpack ϵ, r = p
    return 1 - ϵ*s^r
end




"""
    precip(w,p::Dict{Symbol, Float64})

Computes the mean precipitation rate as a function of water vapour pass, w, and the model parameters
a, b., and saturation water vapour pass, w_sat. 
First introduced by Bretherton et al. (2004), DOI:10.1175/1520-0442(2004)017<1517:RBWVPA>2.0.CO;2

"""
function precip(w,p::Dict{Symbol, Float64})
    @unpack wsat, a, b = p
    return exp(a*(w / wsat - b))
end







