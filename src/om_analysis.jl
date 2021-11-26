function om_derived_quantities!(df_name)

    df = CSV.read(datadir("sims", df_name * ".csv"), DataFrame)

    #remove parameter sets for which no fixed point was found
    dfclean = df[(df.s .> 0.0) .| (df.w1 .> 0.0) .| (df.w2 .> 0.0) .| (df.w3 .> 0.0), : ]
    
    #compute derived quantities
    dfclean.P1 = exp.(dfclean.a .* (dfclean.w1 ./ dfclean.wsat .- dfclean.b))
    dfclean.P2 = exp.(dfclean.a .* (dfclean.w2 ./ dfclean.wsat .- dfclean.b))
    dfclean.P3 = exp.(dfclean.a .* (dfclean.w3 ./ dfclean.wsat .- dfclean.b))
    dfclean.Ptot = dfclean.α .* dfclean.P2 .+ dfclean.L1 ./ dfclean.L .* dfclean.P1 .+ dfclean.L3 ./ dfclean.L .* dfclean.P3
    dfclean.El = dfclean.ep ./2 .* tanh.(dfclean.pt .* (dfclean.s .- (dfclean.spwp .+ dfclean.sfc)./2 ) ) .+ dfclean.ep ./ 2
    dfclean.Φ  = 1 .- dfclean.ϵ .* dfclean.s.^dfclean.r
    dfclean.R  = (1 .- dfclean.Φ) .* dfclean.P2
    dfclean.PR = (dfclean.L1 .+ dfclean.L3) .* dfclean.P2 ./ (dfclean.L1 .* dfclean.P1 + dfclean.L3 .* dfclean.P3)
    
    if occursin(r"v2_closed", df_name) == true
        dfclean = select!(dfclean, Not(:w0))
        dfclean.wo = (dfclean.L1 .* dfclean.w1 + dfclean.L3 .* dfclean.w3) ./ (dfclean.L1 .+ dfclean.L3)
    end

    #convert to better units
    dfclean.L = dfclean.L .* mm2km(1.0)
    dfclean.L1 = dfclean.L1 .* mm2km(1.0)
    dfclean.L2 = dfclean.L2 .* mm2km(1.0)
    dfclean.L3 = dfclean.L3 .* mm2km(1.0)
    dfclean.u = dfclean.u .* (mm2m(1.0) ./ day2s(1.0))
    
    CSV.write(datadir("sims", df_name * "_all_quantities.csv"), dfclean)
    return dfclean
end





"""
    derived_quantities_old!(df::DataFrame)

Compute P1, P2, P3, El, infilt, runoff and PR from the equilibrium solutions for w1, w2, w3 and s 
in a DataFrame and return them as new columns of the input DataFrame.
The input dataframe also needs to contain values for the relevant parameter.
"""
function derived_quantities_old!(df::DataFrame)
        
    df.P1 = precip.(df.w1, df.w_sat, df.a, df.b)
    df.P2 = precip.(df.w2, df.w_sat, df.a, df.b)
    df.P3 = precip.(df.w3, df.w_sat, df.a, df.b)
    df.El = land_evap.(df.s, df.spwp, df.sfc, df.Ep)
    df.infilt = infiltration.(df.s, df.ϵ, df.r)
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


function om_fixedpoints(x0, system)

    p = om_rand_params()
    if system == "v1"
        dynsys = ContinuousDynamicalSystem(open_model_v1, x0, p)
    elseif system == "v2"
        dynsys = ContinuousDynamicalSystem(open_model_v2, x0, p)
    elseif system == "v2_closed"
        dynsys = ContinuousDynamicalSystem(open_model_v2_closed, x0, p)
    else
        println("Assign valid system name: v1, v2 or v2_closed")
    end
    box = om_state_space(p)
    fp, eigs, stable = fixedpoints(dynsys, box)
    p_vals = transpose(collect(values(p)))

    if stable == [0]
        println("Unstable fixed point for parameteres ", collect(keys(p)), "=", p_vals, "!")
    end

    fpm = Matrix(fp)

    if size(fpm, 1) == 0 
        println("No fixed point for parameters ", collect(keys(p)), "=", p_vals, "!")
        solrow = [p_vals 0.0 0.0 0.0 0.0 0.0]
    elseif size(fpm, 1) > 1
        println("More than one fixed point for ", collect(keys(p)), "=", p_vals, "!")
        solrow = [p_vals 0.0 0.0 0.0 0.0 2.0]
    else
        solrow = [p_vals transpose(fpm[1,:]) 1.0]
    end
    return solrow
end





"""
    om_eq_solution(f1, f2, f3, d, nb_cols, w1_lower=0.0, w1_upper=100.0, w2_lower=0.0, w2_upper=100.0, w3_lower=0.0, w3_upper=100.0)

Find roots of the functions f1, f2, f3 and f4 (with variables w1, w2, w3, s) in given intervals 
using parameter values from dictionary d. Output is given in the form of a matrix with 
row(s) of the format [parameter values, solution values for w1 w2 w3 s].
"""
function om_eq_solution(f1, f2, f3, d, nb_cols, w1_lower=0.0, w1_upper=100.0, w2_lower=0.0, w2_upper=100.0, w3_lower=0.0, w3_upper=100.0)

    @unpack spwp, sfc, Ep, eo, ϵ, r, nZr, w0, L, Li, Lo1, Lo2, u, a, b, w_sat = d
    
    #search intervals
    w1_intv = w1_lower..w1_upper
    w2_intv = w2_lower..w2_upper
    w3_intv = w3_lower..w3_upper
    s_dry_intv  = 0.0..spwp
    s_trans_intv = spwp..sfc
    s_wet_intv = sfc..1.0

    #find roots with IntervalRootFinding.jl in the three soil moisture regimes
    
    input_functions = (f1, f2, f3)

    #anonymous function that reduces f1, f2, f3 to functions with one argument 
    #as required by roots() by giving d and t as constants
    #The first element is equivalent to: closure_one(u) = f1(u, d, 0.0)
    closures = [variables -> f(variables,d,0.0) for f in input_functions] 
    sol_dry = roots(closures[1], w1_intv × w2_intv × w3_intv × s_dry_intv)
    sol_trans = roots(closures[2], w1_intv × w2_intv × w3_intv × s_trans_intv)
    sol_wet = roots(closures[3], w1_intv × w2_intv × w3_intv × s_wet_intv)

    #initialise matrix of solutions
    sol_matrix = Array{Float64}(undef, 0, nb_cols)
    parameter_vec = [spwp sfc Ep eo ϵ r nZr w0 L Li Lo1 Lo2 u a b w_sat]

    #store solution(s) in rows of solution matrix

    for i=1:length(sol_dry) #loop through all dry solutions
        if isunique(sol_dry[i]) == true #solution unambiguous
            sol_vec = mid.(interval.(sol_dry))[i]
            w1, w2, w3, s = rou(sol_vec[1],2), rou(sol_vec[2],2), rou(sol_vec[3],2), rou(sol_vec[4],2)
            sol_matrix = [sol_matrix; parameter_vec w1 w2 w3 s]
        else 
            sol_matrix = [sol_matrix; parameter_vec NaN NaN NaN NaN ]
        end
    end

    for i=1:length(sol_trans) 
        if isunique(sol_trans[i]) == true 
            sol_vec = mid.(interval.(sol_trans))[i]
            w1, w2, w3, s = rou(sol_vec[1],2), rou(sol_vec[2],2), rou(sol_vec[3],2), rou(sol_vec[4],2)
            sol_matrix = [sol_matrix; parameter_vec w1 w2 w3 s]
        else 
            sol_matrix = [sol_matrix; parameter_vec NaN NaN NaN NaN ]
        end
    end

    for i=1:length(sol_wet) 
        if isunique(sol_wet[i]) == true 
            sol_vec = mid.(interval.(sol_wet))[i]
            w1, w2, w3, s = rou(sol_vec[1],2), rou(sol_vec[2],2), rou(sol_vec[3],2), rou(sol_vec[4],2)
            sol_matrix = [sol_matrix; parameter_vec w1 w2 w3 s]
        else 
            sol_matrix = [sol_matrix; parameter_vec NaN NaN NaN NaN ]
        end
    end

    if length(sol_matrix) == 0 #no solution at all for these parameter values
        sol_matrix = [sol_matrix; parameter_vec missing missing missing missing ]
    end

    return sol_matrix

    # P_o_mean = (Lo1 * precip(sol_trans_w1, w_sat, a, b) + Lo2 * precip(sol_trans_w3, w_sat, a, b)) / (Lo1 + Lo2)
    # precip_ratio =  precip(sol_trans_w2, w_sat, a, b) / P_o_mean
    # println("Precipitation ratio: ", rou(precip_ratio,4))
    # println("P_1: ", rou(precip(sol_trans_w1, w_sat, a, b),2))
    # println("P_2: ", rou(precip(sol_trans_w2, w_sat, a, b),2))
    # println("P_3: ", rou(precip(sol_trans_w3, w_sat, a, b),2))
    # println("E_l: ", rou(land_evap(sol_trans_s, spwp, sfc, Ep),2))
end



function om_rand_params()

    d = Dict{Symbol, Float64}(
        :spwp => rand(Uniform(0.20, 0.54)), #permanent wilting point
        :ep   => rand(Uniform(4.1, 4.5)),   #[mm/day] potential evaporation over land in mm/day, taken from [1]
        :eo   => rand(Uniform(2.8, 3.2)),   #[mm/day] ocean evaporation rate
        :ϵ    => rand(Uniform(0.9, 1.1)),   #numerical parameter from Rodriguez-Iturbe et al. (1991)
        :r    => 2.0,                       #numerical parameter from Rodriguez-Iturbe et al. (1991)
        :α    => rand(Uniform(0.0, 1.0)),   #land fraction
        :nZr  => rand(Uniform(90.0, 110.0)),#[mm] reservoir depth/"field storage capacity of the soil" [mm] - NEEDS MORE RESEARCH
        :a    => rand(Uniform(11.4, 15.6)), #numerical parameter from Bretherton et al. (2004)
        :b    => rand(Uniform(0.522, 0.603)),   #numerical parameter from Bretherton et al. (2004)
        :wsat => rand(Uniform(65.0, 80.0)), #[mm] saturation water vapour pass derived from plots in Bretherton et al.(2004) - NEEDS MORE RESEARCH
        :u    => rand(Uniform(5.0, 10.0)) * m2mm(1.0)/s2day(1.0), #[mm/day] wind speed
        :L    => 10000.0 * km2mm(1.0),       #rand(Uniform(500.0, 5000.0)) * km2mm(1.0), #[mm] domain size
        :pt   => 10.0,                      #tuning parameter for tanh function
    )

    d[:sfc] = d[:spwp] + 0.3
    d[:w0]  = rand(Uniform(0.0, d[:wsat]))  # water vapour pass at windward model boudnary
    d[:L2]  = d[:α] * d[:L]                 # island length
    d[:L1]  = d[:L]/2 - d[:L2]/2            # length of first ocean [mm], symmetric configuration
    d[:L3]  = d[:L] - d[:L2] - d[:L1]       # length of second ocean [mm]

    return d
end




"""
    om_rand_params_old()

Create random combinations of parameter values for {spwp, sfc, Ep, eo, ϵ, r, nZr, w0, L, Li, Lo1, Lo2, u, a, b, w_sat}.
Return dictionary of parameter names and values. Note that some parameter values are fixed and some are computed from
other randomised parameters.
"""
function om_rand_params_old(w0)
    spwp = rand(Uniform(0.35, 0.55)) #permanent wilting point
    sfc = rand(Uniform(spwp, 0.75))  #field capacity > pwp #NEEDS MORE RESEARCH
    Ep = 4.38  #potential evaporation over land [mm/day], taken from [1]
    eo = 3.0   #ocean evaporation rate [mm/day]
    ϵ = rand(Uniform(0.9, 1.1)) #1.0    #numerical parameter for infiltration
    r = rand(Uniform(1.9, 2.1)) #2.0    #numerical parameter for infiltration
    nZr = 100.0  #reservoir depth/"field storage capacity of the soil" [mm]
    w_sat = rand(Uniform(65.0, 75.0)) #72.0  #saturated water vapour pass for precipitation [mm]
    if w0 == Nothing
        w0 = rand(Uniform(0.0, w_sat)) #advection at outer model boundary [mm]
    end
    L   = km2mm(10000) #length of full model box [mm]
    Li  = km2mm(rand(Uniform(50.0, 300.0))) #length island [mm]
    Lo1 = L/2 - Li/2 #length ocean 1 [mm] - symmetric configuration
    Lo2 = L - Li - Lo1
    u = rand(Uniform(5.0, 10.0)) * m2mm(1)/s2day(1) #wind speed [mm/day]
    a = 15.6    #numerical parameter for precipitation
    b = 0.603   #numerical parameter for precipitation
    return @dict spwp sfc Ep eo ϵ r nZr w0 L Li Lo1 Lo2 u a b w_sat
end

function om_state_space(p::Dict)
    @unpack wsat = p
    s_range  = interval(0.0, 1.0)
    w1_range = interval(0.0, wsat) 
    w2_range = interval(0.0, wsat)
    w3_range = interval(0.0, wsat)
    box = s_range × w1_range × w2_range × w3_range
    return box
end





