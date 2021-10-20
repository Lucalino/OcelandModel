
function closed_model_piecewise(x, p, t)
    
    #This format is required by ContinuousDynamicalSystem.jl algorithm
    @unpack nZr, eo, u, L, α = p
    ds = 1/nZr * (precip(x[2], p) * infiltration(x[1], p) - land_evap(x[1], p))
    dwl = land_evap(x[1], p) - precip(x[2], p) + (x[3] - x[2]) * u / (α*L)
    dwo = eo - precip(x[3], p) - (x[3] - x[2]) * u / ((1-α) * L)
    
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


function closed_model_pw(dx, x, p, t)
    #This format is required by ODEProblem
    @unpack nZr, eo, u, L, α = p
    dx[1] = 1/(nZr) * (precip(x[2], p) * infiltration(x[1], p) - land_evap(x[1], p))
    dx[2] = land_evap(x[1], p) - precip(x[2], p) + (x[3] - x[2]) * u / (α*L)
    dx[3] = eo - precip(x[3], p) - (x[3] - x[2]) * u / ((1-α) * L)
end

function soilmoisture(s,p,t)
    @unpack pa, nZr = p
    #ds =( ( (1-α) * pa + land_evap(s, spwp, sfc, ep) ) * infiltration(s, ϵ,r) - land_evap(s, spwp, sfc, ep) ) / nZr
    ds = (pa * infiltration(s, p) + land_evap(s, p) * (infiltration(s, p) - 1)) / nZr
end