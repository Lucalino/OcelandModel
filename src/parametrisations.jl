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
    land_evapo(spwp::Number,sfc::Number,Ep::Number,s::Number)

Compute the evaporation rate of a land domain.

Input parameters are relative soil moisture saturation s and its values for the permanent wilting point, spwp, and field capacity, sfc,
and potential evaporation, Ep. 
"""
function land_evap(s::Number, spwp::Number,sfc::Number,Ep::Number)     

    if s < spwp
        return 0.0

    elseif s > sfc
        return Ep

    else
        return Ep/(sfc-spwp)*(s-spwp)
    end

end


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
    infiltration(ϵ::Number,r::Number,s::Number)

Compute the fraction of precipitation that infiltrates the soil. 

Take numerical parameters ϵ, r and relative soil moisture saturation s as inputs. 

To obtain the absolute infiltration amount, the result needs to be multiplied by the precipitation amount/rate.

Defined e.g. in Rodriguez‐Iturbe et al. (1991), DOI: 10.1029/91WR01035
"""
function infiltration(s::Number, ϵ::Number, r::Number)
    return 1 - ϵ*s^r
end

"""
    infiltration_diff(ϵ::Number,r::Number,s::Number)

Compute the derivative of the fraction of precipitation that infiltrates the soil with respect to s. 

Take numerical parameters ϵ, r and relative soil moisture saturation s as inputs. 

Infiltration defined e.g. in Rodriguez‐Iturbe et al. (1991), DOI: 10.1029/91WR01035
"""
function infiltration_diff(s::Number, ϵ::Number, r::Number)
    return - r * ϵ * s^(r-1)
end

# 3. Precipitation as a function of water vapour pass
# defined by [2]

# Exponential function with 

"""
    precip(w::Number,w_sat::Number=72.,a::Number=15.6,b::Number=0.603)

Computes the precipitation rate from input arguments water vapour pass, w, 
saturation water vapour pass, w_sat, and two numerical parameters, a and b.

Defined by Bretherton et al. (2004), DOI:10.1175/1520-0442(2004)017<1517:RBWVPA>2.0.CO;2

"""
function precip(w::Number,w_sat::Number,a::Number,b::Number)
    return exp(a*(w / w_sat - b))
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

