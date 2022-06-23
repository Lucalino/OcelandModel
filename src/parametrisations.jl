using Base: Number
# Parametrisations of physical processes

# written by Luca Schmidt on 14.01.2021
# LAST EDITED ON 1. July 2021

# References:
# [1] Rodriguez‐Iturbe (1991), DOI: 10.1029/91WR01035
# [2] Bretherton et al. (2004), DOI:10.1175/1520-0442(2004)017<1517:RBWVPA>2.0.CO;2



# 1. Evapotranspiration E as function of soil moisture saturation s 
# defined by e.g.[INSERT SOURCE HERE!!!]

# Between the pwp and fc limits, E increases linearly with s

# spwp: permanent wilting point
# sfc: field capacity
# Ep: potential evapotranspiration

"""
    land_evapo(s::Float64, p::Dict{::Symbol, ::Float64})

Compute the evaporation rate of a land domain.

Input parameters are relative soil moisture saturation s and its values for the permanent wilting point, spwp, and field capacity, sfc,
and potential evaporation, Ep. 
"""
function El_piecewise(s, p::Dict{Symbol, Float64})     
    @unpack spwp, sfc, ep = p

    if s < spwp
        return 0.0

    elseif s > sfc
        return ep

    else
        return ep/(sfc-spwp)*(s-spwp)
    end

end


function land_evap_dry()
    return 0.0
end

function land_evap_trans(s::Number, spwp::Number,sfc::Number,Ep::Number)
    return Ep/(sfc-spwp)*(s-spwp)
end

function land_evap_wet(Ep::Number)
    return Ep
end


"""
    function evap_tanh(x, p::Dict{Symbol, Real})

Smooth parametrisation of land evapotranspiration that takes three parameters
A = Ep/2, B = tuning_parameter and C = (spwp + sfc)/2 in the form of a dictionary.

"""
function El_tanh(s, p::Dict{Symbol, Float64})
    @unpack spwp, sfc, ep, pt = p
    return ep/2 * tanh( pt * (s - (spwp+sfc)/2 ) ) + ep/2
end
# function evap_tanh(s, spwp, sfc, ep, pt)
#     return ep/2 * tanh( pt * (s - (spwp+sfc)/2 ) ) + ep/2
# end



"""
    land_evap_diff(s::Number, spwp::Number,sfc::Number,Ep::Number)

Compute the derivative of the evaporation rate of a land domain with respect to s.

Input parameters are relative soil moisture saturation s and its values for the permanent wilting point, spwp, and field capacity, sfc,
and potential evaporation, Ep. 
"""
function land_evap_diff(s::Number, spwp::Number,sfc::Number,Ep::Number)

    if s < spwp
        return 0.
    elseif s > sfc
        return 0.
    else
        return Ep / (sfc - spwp)
    end 
end


# 2. Infiltration function Φ as a function of soil moisture saturation s 
# defined by e.g.[1]

# Exponential function with free parameters ϵ and r.
# Φ gets multiplied with the total precipitation over land P_l.
# The subtracted fraction of P_l described by ϵ, r and s constitutes the runoff.

"""
    infiltration(ϵ::Float64, ϵ::Float64, r::Float64)

Compute the fraction of precipitation that infiltrates the soil. 

Take numerical parameters ϵ, r and relative soil moisture saturation s as inputs. 

To obtain the absolute infiltration amount, the result needs to be multiplied by the precipitation amount/rate.

Defined e.g. in Rodriguez‐Iturbe et al. (1991), DOI: 10.1029/91WR01035
"""
function infiltration(s, p::Dict{Symbol, Float64})
    @unpack ϵ, r = p
    return 1 - ϵ*s^r
end
# function infiltration(s::Float64, ϵ::Float64, r::Float64)
#     return 1 - ϵ*s^r
# end




"""
    infiltration_diff(s::Float64, ϵ::Float64, r::Float64)

Compute the derivative of the fraction of precipitation that infiltrates the soil with respect to s. 

Take numerical parameters ϵ, r and relative soil moisture saturation s as inputs. 

Infiltration defined e.g. in Rodriguez‐Iturbe et al. (1991), DOI: 10.1029/91WR01035
"""
function infiltration_diff(s::Float64, ϵ::Float64, r::Float64)
    return - r * ϵ * s^(r-1)
end

# 3. Precipitation as a function of water vapour pass
# defined by [2]

# Exponential function with 

"""
    precip(w::Float64,w_sat::Float64,a::Float64,b::Float64)

Computes the precipitation rate from input arguments water vapour pass, w, 
saturation water vapour pass, w_sat, and two numerical parameters, a and b.

Defined by Bretherton et al. (2004), DOI:10.1175/1520-0442(2004)017<1517:RBWVPA>2.0.CO;2

"""
function precip(w,p::Dict{Symbol, Float64})
    @unpack wsat, a, b = p
    return exp(a*(w / wsat - b))
end

# function precip(w::Float64,w_sat::Float64,a::Float64,b::Float64)
#     return exp(a*(w / w_sat - b))
# end


"""
    precip_pw(w::Float64, p::Dict)    

Computes precipitation as a function of water vapour pass but defined piece-wise (pw) such that 
the precipitation is zero below some critical value w_crit.

"""
function precip_pw(w::Float64, p::Dict)

    @unpack w_crit, a, b, w_sat = p

    if w < w_crit
        return 0.0

    else
        return exp(a*((w - w_crit) / w_sat - b))
    end
end

#the precip_trans function needs more research
function precip_trans(w::Number,w_crit::Number,w_sat::Number,a::Number,b::Number)
    return exp(a*((w - w_crit) / w_sat - b))
end

"""
    precip_from_water_content(w::Number,L::Number,w_sat::Number=72.,a::Number=15.6,b::Number=0.603)

Computes the precipitation rate from input arguments atmospheric water content, w, domain length, L,
saturation water vapour pass, w_sat, and two numerical parameters, a and b.
w_sat is an optional parameter and increases the validity of the function in different regions and across seasons.

Defined by Bretherton et al. (2004), DOI:10.1175/1520-0442(2004)017<1517:RBWVPA>2.0.CO;2

"""
function precip_from_water_content(W::Number,w_sat::Number,a::Number,b::Number)
    return exp(a*(W / (w_sat * L) - b))
end

"""
    precip_diff_from_water_content(w::Number,L::Number,w_sat::Number=72.,a::Number=15.6,b::Number=0.603)

Computes the derivative of the precipitation rate with respect to w from input arguments atmospheric water content, w, domain length, L, 
saturation water vapour pass, w_sat, and two numerical parameters, a and b.
w_sat is an optional parameter and increases the validity of the function in different regions and across seasons.

Defined by Bretherton et al. (2004), DOI:10.1175/1520-0442(2004)017<1517:RBWVPA>2.0.CO;2

"""
function precip_diff(w::Number,w_sat::Number,a::Number,b::Number)
    return a * exp(a*(w / w_sat - b)) / w_sat
end


# Fake linear parametrisations of El(s) and P(w)

function lin_evap(s, p::Dict{Symbol, Float64})
    @unpack ep = p
    return ep * s
end

function lin_precip(w, p::Dict{Symbol, Float64})
    @unpack wsat = p
    return 40.0 / wsat * w
end


