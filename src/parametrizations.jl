function advected_moisture(wl, wo, t, p)
    if wind_DC(t,p) > 0
        return wo
    else
        return wl
    end
end



"""
    density_moist_air(T::Float64, q::Float64, p::Dict{Symbol, Float64})

Compute the density of moist air with specific humidity q in kg/kg and
temperature T in K using the ideal gas law.
"""
function density_moist_air(T::Float64, q::Float64, p::Dict{Symbol, Float64})
    @unpack P0, Rd = p
    return P0 / (Rd * virtual_temperature(T, q))
end


function evap_scaling(t, p)
    @unpack f_a = p
    if f_a == 0.0
        return 1.0
    else
        return f_a * (1 + cos(2π * (t - 0.5)))
    end
end

function wind_DC(t, p::Dict{Symbol, Float64})
    @unpack u_max, t_shift = p
    return u_max * cos(2π * (t - 0.5 - t_shift))
end


"""
    El_tanh(s, p::Dict{Symbol, Float64})

Compute the mean evaporation rate from the land surface as a function of soil moisture, s, and the model parameters 
permanent wilting point, spwp, field capacity, sfc, and potential evaporation, Ep. 
The shape of the function is inspired by a sketch in Seneviratne et al. (2010), DOI:10.1016/j.earscirev.2010.02.004.
"""
function El_tanh(s, p::Dict{Symbol, Any})
    @unpack spwp, sfc, ep, pt = p
    return ep/2 * (1 + tanh( pt * (s - (spwp+sfc)/2 ) ) ) - (ep/2 * (1 + tanh( pt * (-(spwp+sfc)/2) )))
end


#CAUTION: NOT FINISHED
function El_from_LH(T::Float64, s, p::Dict{Symbol, Float64})
    @unpack P0, λ, Rd = p
    q0 = saturation_specific_humidity(T, p)
    q_2m = specific_humidity(2, q0, p)
    ρ0 = density_moist_air(T, q0)
    r_s = 
    
    return density_moist_air()
end



"""
    infiltration(s, p::Dict{Symbol, Float64})

Compute the fraction of land precipitation that infiltrates the soil as a function of soil moisture, s, and the model 
parameters ϵ and r. Defined e.g. in Rodriguez‐Iturbe et al. (1991), DOI: 10.1029/91WR01035.
"""
function infiltration(s, p::Dict{Symbol, Any})
    @unpack ϵ, r = p
    return 1 - ϵ*s^r
end




"""
    precip(w,p::Dict{Symbol, Float64})

Compute the mean precipitation rate as a function of water vapour pass, w, and the model parameters
a, b., and saturation water vapour pass, w_sat. 
First introduced by Bretherton et al. (2004), DOI:10.1175/1520-0442(2004)017<1517:RBWVPA>2.0.CO;2

"""
function precip(w,p::Dict{Symbol, Any})
    @unpack wsat, a, b = p
    return exp(a*(w / wsat - b)) - exp(-a*b)
end

function precip_land(w::Any, p::Dict{Symbol, Any})
    @unpack wsat, a, bland = p
    return exp(a*(w / wsat - bland)) - exp(-a*bland)
end

function precip_oce(w,p::Dict{Symbol, Any})
    @unpack wsat, a, boce = p
    return exp(a*(w / wsat - boce)) - exp(-a*boce)
end




"""
    saturation_pressure(T:Float64)

Compute the saturation water vapor pressure in Pa with the August-Roche-Magnus formula 
which approximates the Clausius-Clapeyron relationship assuming a temperature-independent
specific latent heat of evaporation.
"""
function saturation_pressure(T::Float64)
    return hPa2Pa(6.1094 * exp(17.625*K2dC(T)/(K2dC(T)+243.04)))
end


"""
    saturation_specific_humidity(T::Float64, p::Dict{Symbol, Float64})

Compute the saturation specific humidity in kg/kg.
"""
function saturation_specific_humidity(T::Float64, p::Dict{Symbol, Float64})
    @unpack P0 = p
    return 0.622 * saturation_pressure(T) / P0
end




function soil_resistance(s::Float64, p::Dict{Symbol, Float64})
    @unpack rfc, sfc, spwp = p
    return rfc * (sfc - spwp)/(s - spwp)
end



"""
    specific_humidity(z::Float64, q_l0::Float64, p::Dict{Symbol, Float64})

Compute the specific humidity at height z, starting from a surface value q_0 and assuming a moisture scale height λ.
The shape of the assumed specific humidity profile is taken from Stevens (2007), DOI:10.1175/JAS3983.1.
"""
function specific_humidity(z::Float64, q_0::Float64, p::Dict{Symbol, Float64})
    @unpack λ = p
    return q_0 * exp(- z/λ)
end

"""
    surface_temperature(t::Float64, p::Dict{Symbol, Float64})   

Compute the surface temperature in K at a given time t (in hours since midnight). The parameters T_mean and T_amp,
contained in parameter dictionary p, determine the diurnal mean temperature and temperature amplitude.
"""
function surface_temperature(t::Float64, p::Dict{Symbol, Float64})
    @unpack T_mean, T_amp = p
    return T_mean + T_amp * cos(2π * (t - 0.5))
end


"""
    virtual_temperature(T::Float64, q::Float64)

Compute the virtual temperature for given normal temperature T in K and
specific humidity q in kg/kg.
"""
function virtual_temperature(T::Float64, q::Float64)
    return T * (1 + 0.61 * q)
end

function wind_DC(t::Float64, p::Dict{Symbol, Float64})
    @unpack u_max, t_shift = p
    return u_max * cos(2π * (t - 0.5 - t_shift))
end

function wind_DC_xt(x::Float64, t::Float64, p::Dict{Symbol, Any})
    @unpack α, L, u_max = p
    if (x < (1 - α) * L) & (x > 0)
        u_ocean = u_max * cos(2*pi*t) * cos(pi * x / ((1-α)*L))
        return u_ocean
    else
        u_land = u_max * cos(2*pi*t) * (-cos(pi * (x - (1-α)*L) / (α * L)))
        return u_land
    end
end

function dwind_dx_DC_xt(x::Float64, t::Float64, p::Dict{Symbol, Any})
    @unpack α, L, u_max = p
    if (x < (1 - α) * L) & (x > 0) # ocean box
        u_ocean = u_max * cos(2*pi*t) * pi / ((1-α) * L) * (-sin(pi * x / ((1-α)*L)))
        return u_ocean
    else # land box
        u_land = u_max * cos(2*pi*t) * pi / (α * L) * (sin(pi * (x - (1-α)*L) / (α * L)))
        return u_land
    end
end








