
function closed_model_piecewise(x, p, t)
    
    #This format is required by ContinuousDynamicalSystem.jl algorithm
    @unpack nZr, eo, u, L, α = p
    ds = 1/nZr * (precip(x[2], p) * infiltration(x[1], p) - land_evap(x[1], p))
    dwl = land_evap(x[1], p) - precip(x[2], p) + (x[3] - x[2]) * u / (α*L)
    dwo = eo - precip(x[3], p) - (x[3] - x[2]) * u / ((1-α) * L)
    
    return SVector(ds, dwl, dwo)
end


function closed_model_smooth_tau(x, p, t)
    #This format is required by ContinuousDynamicalSystem.jl algorithm
    
    @unpack α, nZr, τ, eo = p
   
    ds = (precip(x[2], p) * infiltration(x[1], p) - evap_tanh(x[1], p)) / nZr
    dwl = evap_tanh(x[1], p) - precip(x[2], p) + (x[3] - x[2]) * τ /α
    dwo = eo - precip(x[3], p) - (x[3] - x[2]) * τ / (1-α)
    # ds = 1/(nZr) * (precip(x[2], wsat, a, b) * infiltration(x[1], ϵ, r) - evap_tanh(x[1], spwp, sfc, ep, pt))
    # dwl = evap_tanh(x[1], spwp, sfc, ep, pt) - precip(x[2], wsat, a, b) + (x[3] - x[2]) * u / (α*L)
    # dwo = eo - precip(x[3], wsat, a, b) - (x[3] - x[2]) * u / ((1-α) * L)
    
    return SVector(ds, dwl, dwo)
end


function closed_model_smooth(x, p, t)
    #This format is required by ContinuousDynamicalSystem.jl algorithm

    @unpack α, nZr, u, L, eo = p
   
    ds = (precip(x[2], p) * infiltration(x[1], p) - evap_tanh(x[1], p)) / nZr
    dwl = evap_tanh(x[1], p) - precip(x[2], p) + (x[3] - x[2]) * u / (α*L)
    dwo = eo - precip(x[3], p) - (x[3] - x[2]) * u / ((1-α) * L)
    # ds = 1/(nZr) * (precip(x[2], wsat, a, b) * infiltration(x[1], ϵ, r) - evap_tanh(x[1], spwp, sfc, ep, pt))
    # dwl = evap_tanh(x[1], spwp, sfc, ep, pt) - precip(x[2], wsat, a, b) + (x[3] - x[2]) * u / (α*L)
    # dwo = eo - precip(x[3], wsat, a, b) - (x[3] - x[2]) * u / ((1-α) * L)
    
    return SVector(ds, dwl, dwo)
end


function closed_model_linearised(x, p, t)

    @unpack α, nZr, u, L, eo = p
    ds = (lin_precip(x[2], p) * infiltration(x[1], p) - evap_tanh(x[1], p)) / nZr
    dwl = evap_tanh(x[1], p) - lin_precip(x[2], p) + (x[3] - x[2]) * u / (α*L)
    dwo = eo - lin_precip(x[3], p) - (x[3] - x[2]) * u / ((1-α) * L)
    return SVector(ds, dwl, dwo)

end


function closed_model_pw(dx, x, p, t)
    #This format is required by ODEProblem
    @unpack nZr, eo, u, L, α = p
    dx[1] = 1/(nZr) * (precip(x[2], p) * infiltration(x[1], p) - land_evap(x[1], p))
    dx[2] = land_evap(x[1], p) - precip(x[2], p) + (x[3] - x[2]) * u / (α*L)
    dx[3] = eo - precip(x[3], p) - (x[3] - x[2]) * u / ((1-α) * L)
end


function open_model_v1(x, p, t)
    # This model version tracks atmospheric moisture accumulation as a function of location x and computes advection 
    # as boundary value times windspeed rather than using the mean value for advection. The mean values that are the 
    # variables are computed under the assumption that wvp(x) is linear.

    @unpack L1, L2, L3, nZr, u, eo, w0 = p

    # Variable definitions: s = x[1], w1 = x[2], w2 = x[3], w3 = x[4]

    ds  = (precip(x[3], p) * infiltration(x[1], p) - evap_tanh(x[1], p)) / nZr
    dw1 = 2 * (w0 - x[2]) * u/L1 + eo - precip(x[2], p)
    dw2 = 2 * (2*x[2] - w0 - x[3]) * u/L2 + evap_tanh(x[1], p) - precip(x[3], p)
    dw3 = 2 * (2*x[3] - 2*x[2] - x[4] + w0) * u/L3 + eo - precip(x[4], p)

    return SVector(ds, dw1, dw2, dw3)
end


function open_model_v2(x, p, t)
    # This model version is the open analogy to the closed model where we don't track moisture accumulation over space 
    # and instead only consider the mean value for each atmsopheric box. Advection is then this mean value times wind
    # speed at each boundary.

    @unpack L1, L2, L3, nZr, u, eo, w0 = p

    # Variable definitions: s = x[1], w1 = x[2], w2 = x[3], w3 = x[4]

    ds  = (precip(x[3], p) * infiltration(x[1], p) - evap_tanh(x[1], p)) / nZr
    dw1 = eo - precip(x[2], p) + (w0 - x[2]) * u/L1
    dw2 = evap_tanh(x[1], p) - precip(x[3], p) + (x[2] - x[3]) * u/L2
    dw3 = eo - precip(x[4], p) + (x[3] - x[4]) * u/L3

    return SVector(ds, dw1, dw2, dw3)
end


function open_model_v2_closed(x, p, t)
    # This model version is the open analogy to the closed model where we don't track moisture accumulation over space 
    # and instead only consider the mean value for each atmsopheric box. Advection is then this mean value times wind
    # speed at each boundary. In contrast to the open_model_v2, we "close" the model by setting w0 = w3.

    @unpack L1, L2, L3, nZr, u, eo = p

    # Variable definitions: s = x[1], w1 = x[2], w2 = x[3], w3 = x[4]

    ds  = (precip(x[3], p) * infiltration(x[1], p) - evap_tanh(x[1], p)) / nZr
    dw1 = eo - precip(x[2], p) + (x[4] - x[2]) * u/L1
    dw2 = evap_tanh(x[1], p) - precip(x[3], p) + (x[2] - x[3]) * u/L2
    dw3 = eo - precip(x[4], p) + (x[3] - x[4]) * u/L3

    return SVector(ds, dw1, dw2, dw3)
end


function open_model_v2_closed_linearised(x, p, t)
    # This model version is the open analogy to the closed model where we don't track moisture accumulation over space 
    # and instead only consider the mean value for each atmsopheric box. Advection is then this mean value times wind
    # speed at each boundary. In contrast to the open_model_v2, we "close" the model by setting w0 = w3.

    @unpack L1, L2, L3, nZr, u, eo = p

    # Variable definitions: s = x[1], w1 = x[2], w2 = x[3], w3 = x[4]

    ds  = (lin_precip(x[3], p) * infiltration(x[1], p) - evap_tanh(x[1], p)) / nZr
    dw1 = eo - lin_precip(x[2], p) + (x[4] - x[2]) * u/L1
    dw2 = evap_tanh(x[1], p) - lin_precip(x[3], p) + (x[2] - x[3]) * u/L2
    dw3 = eo - lin_precip(x[4], p) + (x[3] - x[4]) * u/L3

    return SVector(ds, dw1, dw2, dw3)
end


function soilmoisture(s,p,t)
    @unpack pa, nZr = p
    #ds =( ( (1-α) * pa + land_evap(s, spwp, sfc, ep) ) * infiltration(s, ϵ,r) - land_evap(s, spwp, sfc, ep) ) / nZr
    ds = (pa * infiltration(s, p) + land_evap(s, p) * (infiltration(s, p) - 1)) / nZr
end