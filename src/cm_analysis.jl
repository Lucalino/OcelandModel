# Functions used in the closed model analysis
# File created on 8th October 2021 by Luca Schmidt

"""
    derived_quantities!(df::DataFrame)

Compute P1, P2, P3, El, infilt, runoff and PR from the equilibrium solutions for s, wl and wo 
in a DataFrame and return them as new columns of the input DataFrame.
The input dataframe also needs to contain values for the relevant parameter.
"""
function derived_quantities!(df::DataFrame)

    d = Dict{Symbol, Float64}(
        :spwp => df.spwp, 
        :sfc  => df.sfc,
        :ep   => df.ep,    
        :eo   => df.eo,     
        :ϵ    => df.ϵ,     
        :r    => df.r, 
        :α    => df.α,
        :nZr  => df.nZr,
        :a    => df.a,
        :b    => df.b,
        :wsat => df.wsat,
        :u    => df.u,
        :L    => df.L,
        :pt   => df.pt,
    )

    df.Pl = precip.(df.wl, d)
    df.Po = precip.(df.wo, d)
    #df.El = land_evap.(df.s, d)
    df.El = evap_tanh(df.s, d)
    df.infilt = infiltration.(df.s, d)
    df.runoff = (1 .- df.infilt) .* df.P2
    df.PR = df.P2 .* (df.Lo1 + df.Lo2) ./ (df.Lo1 .* df.P1 + df.Lo2 .* df.P3)
    df.EminP1 = df.eo - df.P1
    df.EminP2 = df.El - df.P2
    df.EminP3 = df.eo - df.P3

    #convert to better units
    df.L = df.L * mm2km(1)
    df.Li = df.Li * mm2km(1)
    df.Lo1 = df.Lo1 * mm2km(1)
    df.Lo2 = df.Lo2 * mm2km(1)
    df.u = df.u * (mm2m(1) / day2s(1))
    return df
end



function cm_eq_nlsolve(x0)

    p = cm_rand_params()    
    system = (F, x) -> f!(F, x, p)    
    solution = nlsolve(system, x0)    
    sol_vec = transpose(solution.zero)
    p_vals = transpose(collect(values(p)))
    
    if converged(solution) == true
        convergence = 1
    else
        convergence = 0
    end

    sol_row = [p_vals sol_vec convergence]

    return sol_row
end

function cm_eq_fixedpoints(system, x0)
    p = cm_rand_params() 

    if system == "smooth"
        dynsys = ContinuousDynamicalSystem(closed_model_smooth, x0, p)
    elseif system == "piecewise"
        dynsys = ContinuousDynamicalSystem(closed_model_piecewise, x0, p)
    else
        println("Specification of ODE system missing (smooth or piecewise).")
    end

    box = cm_state_space(p)
    fp, eigs, stable = fixedpoints(dynsys, box)
    p_vals = transpose(collect(values(p)))

    if stable == [0]
        println("Unstable fixed point for parameteres ", col_names, "=", p_vals, "!")
    end

    fpm = Matrix(fp)

    if size(fpm, 1) == 0 
        println("No fixed point for parameters ", col_names, "=", p_vals, "!")
    elseif size(fpm, 1) > 1
        println("More than one fixed point for ", col_names, "=", p_vals, "!")
    else
        solrow = [p_vals transpose(fpm[1,:])]
    end
    return solrow
end


function cm_fixed_params(n::Int)
    d = Dict{Symbol, Float64}(
        :spwp => 0.3,     #permanent wilting point
        :ep   => 4.38,    #[mm/day] potential evaporation over land in mm/day, taken from [1]
        :pt   => 10.0,    #tuning parameter in tanh function
        :ptot => 3.0,     #[mm/day] average total precipitation rate over the full tropics
        :eo   => 3.0,     #[mm/day] ocean evaporation rate
        :ϵ    => 1.0,     #numerical parameter from Rodriguez-Iturbe et al. (1991)
        :r    => 2.0,     #numerical parameter from Rodriguez-Iturbe et al. (1991)
        :α    => 0.1,     #land fraction
        :nZr  => 100.0,   #[mm] reservoir depth/"field storage capacity of the soil"
        :a    => 15.6,    #numerical parameter from Bretherton et al. (2004)
        :b    => 0.603,   #numerical parameter from Bretherton et al. (2004)
        :wsat => 72.0,    #[mm] saturation water vapour pass derived from plots in Bretherton et al. ( 2004)
        :u    => 5.0 * m2mm(1.0)/s2day(1.0), #[mm/day] wind speed
        :L    => 1000.0 * km2mm(1.0), #[mm] domain size
    )
    d[:sfc] = d[:spwp] + 0.3 #field capacity
    pa = [0.5, 0.7, 0.9]
    d[:pa] = pa[n]     #[mm/day] advected precipitation component

    return d
end


function cm_state_space(p::Dict)
    @unpack wsat = p
    s_range = interval(0.0,1.0)
    wl_range = interval(0.0,wsat) 
    wo_range = interval(0.0,wsat)
    box = s_range × wl_range × wo_range
    return box
end






function cm_rand_params()

    d = Dict{Symbol, Float64}(
        :spwp => rand(Uniform(0.20, 0.54)), #permanent wilting point
        :ep   => rand(Uniform(4.1, 4.5)),     #[mm/day] potential evaporation over land in mm/day, taken from [1]
        :eo   => rand(Uniform(2.8, 3.2)),     #[mm/day] ocean evaporation rate
        :ϵ    => rand(Uniform(0.9, 1.1)),      #numerical parameter from Rodriguez-Iturbe et al. (1991)
        :r    => 2.0, #rand(Uniform(1.9, 2.1)),      #numerical parameter from Rodriguez-Iturbe et al. (1991)
        :α    => rand(Uniform(0.1, 0.7)),       #land fraction
        :nZr  => rand(Uniform(90.0, 110.0)), #[mm] reservoir depth/"field storage capacity of the soil" [mm] - NEEDS MORE RESEARCH
        :a    => rand(Uniform(11.4, 15.6)),    #numerical parameter from Bretherton et al. (2004)
        :b    => rand(Uniform(0.522, 0.603)),   #numerical parameter from Bretherton et al. (2004)
        :wsat => rand(Uniform(65.0, 80.0)), #[mm] saturation water vapour pass derived from plots in Bretherton et al.(2004) - NEEDS MORE RESEARCH
        :u    => rand(Uniform(5.0, 10.0)) * m2mm(1.0)/s2day(1.0), #[mm/day] wind speed
        :L    => rand(Uniform(500.0, 5000.0)) * km2mm(1.0), #[mm] domain size
        :pt   => 10.0, #tuning parameter for tanh function
    )

    d[:sfc] = d[:spwp] + 0.3

    return d
end

