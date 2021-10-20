"""
    derived_quantities!(df::DataFrame)

Compute P1, P2, P3, El, infilt, runoff and PR from the equilibrium solutions for w1, w2, w3 and s 
in a DataFrame and return them as new columns of the input DataFrame.
The input dataframe also needs to contain values for the relevant parameter.
"""
function derived_quantities!(df::DataFrame)
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

function mean_of_bins!(df::DataFrame, param::String, var::String, nb_bins::Int)
    df = sort(df, param)
    col_name = string("mean_",var)
    df[!,col_name] .= 1.0
    nb_rows = nrow(df)
    mod(nb_rows, nb_bins) == 0 || error("Number of bins must divide the data set's number of rows")
    bin_size = Int(nb_rows / nb_bins)
    bin_start = 1  

    for bin_end in Iterators.countfrom(bin_size, bin_size)
        bin_end > nb_rows && break
        df[bin_start:bin_end, col_name] .= mean(df[bin_start:bin_end, var])
        bin_start = bin_end+1
    end
    return df[!,col_name]
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



"""
    om_rand_params()

Create random combinations of parameter values for {spwp, sfc, Ep, eo, ϵ, r, nZr, w0, L, Li, Lo1, Lo2, u, a, b, w_sat}.
Return dictionary of parameter names and values. Note that some parameter values are fixed and some are computed from
other randomised parameters.
"""
function om_rand_params(w0)
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





