"""
    closed_model_τ(x, p, t)

This model version is used in the paper. Closed model equations where wind speed u and 
domain length L are combined to form the atmospheric transport parameter τ. 
"""
function closed_model_τ(x, p, t)
        
    @unpack α, nZr, τ, eo = p
   
    ds = (precip(x[2], p) * infiltration(x[1], p) - El_tanh(x[1], p)) / nZr
    dwl = El_tanh(x[1], p) - precip(x[2], p) + (x[3] - x[2]) * τ /α
    dwo = eo - precip(x[3], p) - (x[3] - x[2]) * τ / (1-α)
    
    return SVector(ds, dwl, dwo)
end




"""
    closed_model_uL(x, p, t)

Closed model equations where wind speed u and domain length L are kept as separate parameters
which can be varied independently of each other.
"""
function closed_model_uL(x, p, t)
    
    @unpack α, nZr, u, L, eo = p
   
    ds = (precip(x[2], p) * infiltration(x[1], p) - El_tanh(x[1], p)) / nZr
    dwl = El_tanh(x[1], p) - precip(x[2], p) + (x[3] - x[2]) * u / (α*L)
    dwo = eo - precip(x[3], p) - (x[3] - x[2]) * u / ((1-α) * L)
    
    return SVector(ds, dwl, dwo)
end


function closed_model_DC_u!(du, u, p, t)
    @unpack nZr, eo, α, L = p
    du[1] = (precip(u[2], p) * infiltration(u[1], p) - evap_scaling(t, p) * El_tanh(u[1], p)) / nZr
    du[2] = evap_scaling(t, p) * El_tanh(u[1], p) - precip(u[2], p) + 2 * advected_moisture(u[2], u[3], t, p) * wind_DC(t, p) / (α * L)
    du[3] = eo - precip(u[3], p) - 2 * advected_moisture(u[2], u[3], t, p) * wind_DC(t, p) / ((1-α) * L)
    #a = eo - precip(u[3], p)
    #b = 2 * advected_moisture(u[2], u[3], t, p) * wind_DC(t, p) / ((1-α) * L)
    #@show a, b
end


"""
    open_model_v_paper(x, p, t)

This model version was used in the paper. It is the open analogy to the closed model
where advection fluxes are computed from the difference between the windward and leeward boxes' 
mean water vapor passes. 
The moisture variables are defined as follows: s = x[1], w1 = x[2], w2 = x[3], w3 = x[4].
"""
function open_model_v_paper(x, p, t)
    @unpack L1, L2, L3, nZr, u, eo, w0 = p

    ds  = (precip(x[3], p) * infiltration(x[1], p) - El_tanh(x[1], p)) / nZr
    dw1 = eo - precip(x[2], p) + (w0 - x[2]) * u/L1
    dw2 = El_tanh(x[1], p) - precip(x[3], p) + (x[2] - x[3]) * u/L2
    dw3 = eo - precip(x[4], p) + (x[3] - x[4]) * u/L3

    return SVector(ds, dw1, dw2, dw3)
end



"""
    open_model_v2_closed(x, p, t)

This model shares the box configuration with the open model used in the paper, but is modified 
such that moisture which leaves at the leeward boundary is fed back into the model through 
the windward boundary (w0 = w3). Thus, the model is closed and can be used to analyze whether 
we recover the closed model behavior if we close the open model.
The moisture variables are defined as follows: s = x[1], w1 = x[2], w2 = x[3], w3 = x[4].
"""
function open_model_closed(x, p, t)
    @unpack L1, L2, L3, nZr, u, eo = p

    ds  = (precip(x[3], p) * infiltration(x[1], p) - El_tanh(x[1], p)) / nZr
    dw1 = eo - precip(x[2], p) + (x[4] - x[2]) * u/L1
    dw2 = El_tanh(x[1], p) - precip(x[3], p) + (x[2] - x[3]) * u/L2
    dw3 = eo - precip(x[4], p) + (x[3] - x[4]) * u/L3

    return SVector(ds, dw1, dw2, dw3)
end




"""
    open_model_w_tracked(x, p, t)

Open model configuration with the assumption that the water vapor pass changes linearly across 
each model box in horizontal direction, depending on whether mean evapo(transpi)ration or 
precipitation is stronger. The model equations describe the time evolution of the mean water vapor 
pass and soil moisture in the different model boxes but atmospheric moisture that crosses the
boundary between two boxes (advection) is computed from the w-value at this boundary which is
generally not identical with the mean value of the box. 
The moisture variables are defined as follows: s = x[1], w1 = x[2], w2 = x[3], w3 = x[4].
"""
function open_model_w_tracked(x, p, t)
    @unpack L1, L2, L3, nZr, u, eo, w0 = p

    ds  = (precip(x[3], p) * infiltration(x[1], p) - El_tanh(x[1], p)) / nZr
    dw1 = 2 * (w0 - x[2]) * u/L1 + eo - precip(x[2], p)
    dw2 = 2 * (2*x[2] - w0 - x[3]) * u/L2 + El_tanh(x[1], p) - precip(x[3], p)
    dw3 = 2 * (2*x[3] - 2*x[2] - x[4] + w0) * u/L3 + eo - precip(x[4], p)

    return SVector(ds, dw1, dw2, dw3)
end




