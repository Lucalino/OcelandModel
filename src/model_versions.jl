"""
    cm_τ(x, p, t)

This model verion is used in the paper. Closed model equations where wind speed u and 
domain length L are combined to form the atmospheric transport parameter τ. 
"""
function cm_τ(x, p, t)
        
    @unpack α, nZr, τ, eo = p
   
    ds = (precip(x[2], p) * infiltration(x[1], p) - El_tanh(x[1], p)) / nZr
    dwl = El_tanh(x[1], p) - precip(x[2], p) + (x[3] - x[2]) * τ /α
    dwo = eo - precip(x[3], p) - (x[3] - x[2]) * τ / (1-α)
    
    return SVector(ds, dwl, dwo)
end



"""
    cm_τ_Pl_neq_Po(x, p, t)

Closed model version with different precipitation parametrizations P_l(w) and P_o(w) over land and ocean, respectivly. 
Everything else is identical with the closed model version used in the paper (cm_τ(x, p, t)).
"""
function cm_τ_Pl_neq_Po(x, p, t)
        
    @unpack α, nZr, τ, eo = p
   
    ds = (precip_land(x[2], p) * infiltration(x[1], p) - El_tanh(x[1], p)) / nZr
    dwl = El_tanh(x[1], p) - precip_land(x[2], p) + (x[3] - x[2]) * τ /α
    dwo = eo - precip_oce(x[3], p) - (x[3] - x[2]) * τ / (1-α)
    
    return SVector(ds, dwl, dwo)
end




"""
    cm_uL(x, p, t)

Closed model equations where wind speed u and domain length L are kept as separate parameter
which can be varied independently of each other.
"""
function cm_uL(x, p, t)
    
    @unpack α, nZr, u, L, eo = p
   
    ds = (precip(x[2], p) * infiltration(x[1], p) - El_tanh(x[1], p)) / nZr
    dwl = El_tanh(x[1], p) - precip(x[2], p) + (x[3] - x[2]) * u / (α*L)
    dwo = eo - precip(x[3], p) - (x[3] - x[2]) * u / ((1-α) * L)
    
    return SVector(ds, dwl, dwo)
end


"""
    cm_DC_wind!(du, u, p, t)

Closed model equations where the state variables are box-averaged integrated water vapor path, u[2] and u[3], 
over land and ocean, respectively, and relative soil moisture saturation u[1]. The boundary wind value oscillates
between -u_max and u_max to mimick a diurnal cycle in the form of sea and land breezes.
"""
function cm_DC_wind!(du, u, p, t)
    @unpack nZr, eo, α, L = p
    du[1] = (precip(u[2], p) * infiltration(u[1], p) - evap_scaling(t, p) * El_tanh(u[1], p)) / nZr
    du[2] = evap_scaling(t, p) * El_tanh(u[1], p) - precip(u[2], p) +  wind_DC(t, p) * (u[3] - u[2]) / (α * L)
    du[3] = eo - precip(u[3], p) - wind_DC(t, p) * (u[3] - u[2]) / ((1-α) * L)
end

function cm_DC_w_xt!(dw, w, parr, t)
    #@unpack dx, x = p
    dx = parr[1]
    x = parr[2]
    # Once I can let parr be a dictionary again, I can use the unpack command again and just have
    # one parameter container p which is a dictionary defined before calling ODEProblem.
    pdict = cm_fixed_params(false)  
    E = 3.0

    # Set negative moisture values to zero before computing the next moisture derivatives
    for i = 1:length(w)
        if w[i] < 0
            w[i] = 0
        end
    end  

    dwdx = diff(vcat(w, [w[1]])) / dx
    wind = wind_DC_xt.(x,t,Ref(pdict))
    dwinddx = dwind_dx_DC_xt.(x,t,Ref(pdict))

    dw .= E .- precip.(w, Ref(pdict)) .- wind .* dwdx .- w .* dwinddx 

end



"""
    cm_lin(param)

Closed model version akin to the version used in the paper but with purely linear flux parametrizations. 
The equilibrium solution was found analytically.
"""
function cm_lin(param)
    @unpack eo, e, p, r, τ, α = param

    s  = (-1) * ((-eo * r * τ + e * (-1 + α) * τ + eo * r * α * τ + sqrt((-1 + α) * τ * (e^2 * (-1 + α) * τ + eo^2 * r^2 * (-1 + α) * τ +  2 * eo * e * r * (2 * p * (-1 + α) * α - (1 + α) \ τ)))) / (2 * e * r * α * (p * (-1 + α) - τ)))
    wl = (-1) * ((eo * r * τ + e * (-1 + α) * τ - eo * r * α * τ + sqrt((-1 + α) * τ * (e^2 * (-1 + α) * τ + eo^2 * r^2 * (-1 + α) * τ + 2 * eo * e * r * (2 * p * (-1 + α) * α - (1 + α) \ τ)))) / (2 * p * r * (p * (-1 + α) * α - τ)))
    wo = (eo * r * (-1 + α) * (2 * p^2 * (-1 + α) * α - 2 * p * τ - τ^2) + τ * (e * (-1 + α) * τ + sqrt((-1 + α) * τ * (e^2 * (-1 + α) * τ + eo^2 * r^2 * (-1 + α) * τ + 2 * eo * e * r * (2 * p * (-1 + α) * α - (1 + α) \ τ))))) / (2 * p * r * (p * (-1 + α) - τ) * (p * (-1 + α) * α - τ))

    #Formally equilibrium solutions but unphysical:
    #s2  = (eo + r * τ - eo * r * α * τ + e * (τ - α * τ) + sqrt((-1 + α) * τ * (e^2 * (-1 + α) * τ + eo^2 * r^2 * (-1 + α) * τ + 2 * eo * e * r * (2 * p * (-1 + α) * α - (1 + α) \ τ)))) / (2 * e * r * α * (p * (-1 + α) - τ))
    #wl2 = (-eo * r * τ + eo * r * α * τ + e * (τ - α * τ) + sqrt((-1 + α) * τ * (e^2 * (-1 + α) * τ + eo^2 * r^2 * (-1 + α) * τ + 2 * eo * e * r * (2 * p * (-1 + α) * α - (1 + α) \ τ)))) / (2 * p * r * (p * (-1 + α)* α - τ))
    #wo2 = (eo * r * (-1 + α) * (2 * p^2 * (-1 + α) * α - 2 * p * τ - τ^2) + τ * (e * (-1 + α) * τ - sqrt((-1 + α) * τ * (e^2 * (-1 + α) * τ + eo^2 * r^2 * (-1 + α) * τ + 2 * eo * e * r * (2 * p * (-1 + α) * α - (1 + α) \ τ))))) / (2 * p * r * (p * (-1 + α) - τ) * (p * (-1 + α) * α - τ))
    return s, wl, wo   
end



"""
    om_v_paper(x, p, t)

This model verion was used in the paper. It is the open analogy to the closed model
where advection fluxes are computed from the difference between the windward and leeward boxes' 
mean water vapor passes. 
The moisture variables are defined as follows: s = x[1], w1 = x[2], w2 = x[3], w3 = x[4].
"""
function om_v_paper(x, p, t)
    @unpack L1, L2, L3, nZr, u, eo, w0 = p

    ds  = (precip(x[3], p) * infiltration(x[1], p) - El_tanh(x[1], p)) / nZr
    dw1 = eo - precip(x[2], p) + (w0 - x[2]) * u/L1
    dw2 = El_tanh(x[1], p) - precip(x[3], p) + (x[2] - x[3]) * u/L2
    dw3 = eo - precip(x[4], p) + (x[3] - x[4]) * u/L3

    return SVector(ds, dw1, dw2, dw3)
end



"""
    om_v2_closed(x, p, t)

This model shares the box configuration with the open model used in the paper, but is modified 
such that moisture which leaves at the leeward boundary is fed back into the model through 
the windward boundary (w0 = w3). Thus, the model is closed and can be used to analyze whether 
we recover the closed model behavior if we close the open model.
The moisture variables are defined as follows: s = x[1], w1 = x[2], w2 = x[3], w3 = x[4].
"""
function om_closed(x, p, t)
    @unpack L1, L2, L3, nZr, u, eo = p

    ds  = (precip(x[3], p) * infiltration(x[1], p) - El_tanh(x[1], p)) / nZr
    dw1 = eo - precip(x[2], p) + (x[4] - x[2]) * u/L1
    dw2 = El_tanh(x[1], p) - precip(x[3], p) + (x[2] - x[3]) * u/L2
    dw3 = eo - precip(x[4], p) + (x[3] - x[4]) * u/L3

    return SVector(ds, dw1, dw2, dw3)
end




"""
    om_w_tracked(x, p, t)

Open model configuration with the assumption that the water vapor pass changes linearly across 
each model box in horizontal direction, depending on whether mean evapo(transpi)ration or 
precipitation is stronger. The model equations describe the time evolution of the mean water vapor 
pass and soil moisture in the different model boxes but atmospheric moisture that crosses the
boundary between two boxes (advection) is computed from the w-value at this boundary which is
generally not identical with the mean value of the box. 
The moisture variables are defined as follows: s = x[1], w1 = x[2], w2 = x[3], w3 = x[4].
"""
function om_w_tracked(x, p, t)
    @unpack L1, L2, L3, nZr, u, eo, w0 = p

    ds  = (precip(x[3], p) * infiltration(x[1], p) - El_tanh(x[1], p)) / nZr
    dw1 = 2 * (w0 - x[2]) * u/L1 + eo - precip(x[2], p)
    dw2 = 2 * (2*x[2] - w0 - x[3]) * u/L2 + El_tanh(x[1], p) - precip(x[3], p)
    dw3 = 2 * (2*x[3] - 2*x[2] - x[4] + w0) * u/L3 + eo - precip(x[4], p)

    return SVector(ds, dw1, dw2, dw3)
end




