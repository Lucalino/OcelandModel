# Functions used in the closed model analysis
# File created on 8th October 2021 by Luca Schmidt


function cm_alphavar_params()

    d = Dict{Symbol, Float64}(
        :spwp => rand((0.3, 0.4)),     #permanent wilting point
        :ep   => 4.38,    #[mm/day] potential evaporation over land in mm/day, taken from [1]
        :pt   => 10.0,    #tuning parameter in tanh function
        :eo   => 3.0,     #[mm/day] ocean evaporation rate
        :ϵ    => 1.0,     #numerical parameter from Rodriguez-Iturbe et al. (1991)
        :r    => 2.0,     #numerical parameter from Rodriguez-Iturbe et al. (1991)
        :α    => rand(Uniform(0.0, 1.0)),     #land fraction
        :nZr  => 100.0,   #[mm] reservoir depth/"field storage capacity of the soil"
        :a    => 15.6,    #numerical parameter from Bretherton et al. (2004)
        :b    => 0.603,   #numerical parameter from Bretherton et al. (2004)
        :wsat => 72.0,    #[mm] saturation water vapour pass derived from plots in Bretherton et al. ( 2004)
        :u    => rand((5.0, 10.0)) * m2mm(1.0)/s2day(1.0), #[mm/day] wind speed
        :L    => rand((1000.0, 10000.0)) * km2mm(1.0), #[mm] domain size
    )

    d[:sfc] = d[:spwp] + 0.3 #field capacity
    return d
end



"""
    cm_derived_quantities!(df::DataFrame)

Compute P1, P2, P3, El, infilt, runoff and PR from the equilibrium solutions for s, wl and wo 
in a DataFrame and return them as new columns of the input DataFrame.
The input dataframe also needs to contain values for the relevant parameter.
"""
function cm_derived_quantities!(df_name)

    df = CSV.read(datadir("sims", df_name * ".csv"), DataFrame)
    
    df.Pl = exp.(df.a .* (df.wl ./ df.wsat .- df.b))
    df.Po = exp.(df.a .* (df.wo ./ df.wsat .- df.b))
    df.Ptot = df.α .* df.Pl .+ (1 .- df.α) .* df.Po
    df.El = df.ep./2 .* tanh.(df.pt .* (df.s .- (df.spwp .+ df.sfc)./2 ) ) .+ df.ep ./ 2
    df.Φ = 1 .- df.ϵ .* df.s.^df.r
    df.R = (1 .- df.Φ) .* df.Pl
    df.PR = df.Pl ./ df.Po
    df.dw = df.wo .- df.wl

    if occursin(r"tau", df_name) == true
        df.A  = (df.wo .- df.wl) .* df.τ ./ df.α
        df.B  = (df.wo .- df.wl) .* df.τ ./ (1.0 .- df.α)
        df.τ_time = 1 ./ df.τ

    elseif occursin(r"tau", df_name) == false
        df.A  = (df.wo .- df.wl) .* df.u ./ (df.α .* df.L)
        df.B  = (df.wo .- df.wl) .* df.u ./ ((1.0 .- df.α) .* df.L)
        
        #convert to better units
        df.Lkm = df.L .* mm2km(1.0)
        df.Likm = df.α .* df.Lkm
        df.ums = df.u .* (mm2m(1.0) ./ day2s(1.0))
    end
    
    df.sblc = (df.Pl .* df.Φ .- df.El) ./ df.nZr
    df.lblc = df.El .- df.Pl .+ df.A
    df.oblc = df.eo .- df.Po .- df.B

    CSV.write(datadir("sims", df_name * "_all_quantities.csv"), df)
    return df
end

function cm_derived_quantities_lin!(df_name)

    df = CSV.read(datadir("sims", df_name * ".csv"), DataFrame)
    
    df.Pl = 40.0 ./ df.wsat .* df.wl
    df.Po = 40.0 ./ df.wsat .* df.wo
    df.Ptot = df.α .* df.Pl .+ (1 .- df.α) .* df.Po
    df.El = df.ep./2 .* tanh.(df.pt .* (df.s .- (df.spwp .+ df.sfc)./2 ) ) .+ df.ep ./ 2
    df.Φ = 1 .- df.ϵ .* df.s.^df.r
    df.R = (1 .- df.Φ) .* df.Pl
    df.PR = df.Pl ./ df.Po
    df.dw = df.wo .- df.wl
    df.A  = (df.wo .- df.wl) .* df.u ./ (df.α .* df.L)
    df.B  = (df.wo .- df.wl) .* df.u ./ ((1.0 .- df.α) .* df.L)
    df.sblc = (df.Pl .* df.Φ .- df.El) ./ df.nZr
    df.lblc = df.El .- df.Pl .+ df.A
    df.oblc = df.eo .- df.Po .- df.B

    #convert to better units
    df.Lkm = df.L .* mm2km(1.0)
    df.Likm = df.α .* df.Lkm
    df.ums = df.u .* (mm2m(1.0) ./ day2s(1.0))
    CSV.write(datadir("sims", df_name * "_all_quantities.csv"), df)
    return df
end



function cm_eq_fixedpoints(system, x0, tau::Bool=false)
    #p = cm_alphavar_params()
    p = cm_rand_params(tau) 

    if system == "smooth"
        if tau == false
            dynsys = ContinuousDynamicalSystem(closed_model_smooth, x0, p)
            #dynsys = ContinuousDynamicalSystem(closed_model_linearised, x0, p)
        elseif tau == true
            dynsys = ContinuousDynamicalSystem(closed_model_smooth_tau, x0, p)
        end
    elseif system == "piecewise"
        dynsys = ContinuousDynamicalSystem(closed_model_piecewise, x0, p)
    else
        println("Specification of ODE system missing (smooth or piecewise).")
    end

    box = cm_state_space(p)
    fp, eigs, stable = fixedpoints(dynsys, box)
    p_vals = transpose(collect(values(p)))

    if stable == [0]
        println("Unstable fixed point for parameteres ", transpose(collect(keys(p))), "=", p_vals, "!")
    end

    fpm = Matrix(fp)

    if size(fpm, 1) == 0 
        println("No fixed point for parameters ", transpose(collect(keys(p))), "=", p_vals, "!")
    elseif size(fpm, 1) > 1
        println("More than one fixed point for ", transpose(collect(keys(p))), "=", p_vals, "!")
    else
        solrow = [p_vals transpose(fpm[1,:])]
    end
    return solrow
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


function cm_fixed_params(n::Int,tau=false)

    if tau == false

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

    elseif tau == true

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
            :τ    => 0.5
        )
    else
        println("You have not specified which parameter set to use. For tau-set, give second argument 'true'.")
    end

    d[:sfc] = d[:spwp] + 0.3 #field capacity
    pa = [0.5, 0.7, 0.9]
    d[:pa] = pa[n]     #[mm/day] advected precipitation component

    return d
end


function cm_mean_of_bins!(df::DataFrame, param::String, var::String, nb_bins::Int)
    df = sort!(df, param)
    col_name = string("mean_",var)
    df[!,col_name] .= 1.0
    nb_rows = nrow(df)
    mod(nb_rows, nb_bins) == 0 || error("Number of bins must divide the data set's number of rows")
    bin_size = Int(nb_rows / nb_bins)
    bin_start = 1  

    for bin_end in Iterators.countfrom(bin_size, bin_size)
        bin_end > nb_rows && break
        df[bin_start:bin_end, col_name] .= DataFrames.mean(df[bin_start:bin_end, var])
        bin_start = bin_end+1
    end
    return df[!,col_name]
end


function cm_rolling_average!(df::DataFrame, param::String, var::String, n::Int)
    df = sort!(df, param)
    col_name = string("rav_", var)
    df[!, col_name] = [i < n ? DataFrames.mean(df[begin:i, var]) : DataFrames.mean(df[i-n+1:i, var]) for i in 1:length(df[!,var])]
    return df[!,col_name]
end



function cm_state_space(p::Dict)
    @unpack wsat = p
    s_range = interval(0.0,1.0)
    wl_range = interval(0.0,wsat) 
    wo_range = interval(0.0,wsat)
    box = s_range × wl_range × wo_range
    return box
end



function cm_rand_params(tau::Bool=false)

    if tau == false

        d = Dict{Symbol, Float64}(
            :spwp => rand(Uniform(0.15, 0.54)), #permanent wilting point
            :ep   => rand(Uniform(4.0, 6.0)),     #[mm/day] potential evaporation over land in mm/day, taken from Entekhabi et al. (1992)
            :eo   => rand(Uniform(2.5, 3.5)),     #[mm/day] ocean evaporation rate, lower limit from Kumar et al. (2017), lower limit motivated by Zang et al. (1995)
            :ϵ    => rand(Uniform(0.9, 1.1)),      #numerical parameter from Rodriguez-Iturbe et al. (1991)
            :r    => 2, #rand(Uniform(1.9, 2.1)),      #numerical parameter from Rodriguez-Iturbe et al. (1991)
            :α    => rand(Uniform(0.0, 1.0)),       #land fraction
            :nZr  => rand(Uniform(50.0, 150.0)), #[mm] reservoir depth/"field storage capacity of the soil" [mm] - NEEDS MORE RESEARCH
            :a    => rand(Uniform(11.4, 15.6)),    #numerical parameter from Bretherton et al. (2004)
            :b    => rand(Uniform(0.522, 0.603)),   #numerical parameter from Bretherton et al. (2004)
            :wsat => rand(Uniform(65.0, 80.0)), #[mm] saturation water vapour pass derived from plots in Bretherton et al.(2004) - NEEDS MORE RESEARCH
            :u    => rand(Uniform(5.0, 10.0)) * m2mm(1.0)/s2day(1.0), #[mm/day] wind speed
            :L    => 10000.0 * km2mm(1.0), #[mm] domain size
            :pt   => 10.0, #tuning parameter for tanh function
        )
    
    elseif tau == true

        d = Dict{Symbol, Float64}(
            :spwp => rand(Uniform(0.15, 0.55)),  #permanent wilting point, taken from Hagemann & Stacke (2015)
            :ep   => rand(Uniform(4.0, 6.0)),    #[mm/day] potential evaporation over land in mm/day, taken from Entekhabi et al. (1992)
            :eo   => rand(Uniform(2.5, 3.5)),    #[mm/day] ocean evaporation rate, lower limit from Kumar et al. (2017), lower limit motivated by Zang et al. (1995)
            :ϵ    => rand(Uniform(0.9, 1.1)),    #numerical parameter from Rodriguez-Iturbe et al. (1991)
            :r    => rand(Uniform(2,6)),         #numerical parameter from Rodriguez-Iturbe et al. (1991) and Entekhabi et al. (1992)
            :α    => rand(Uniform(0.0, 1.0)),    #land fraction
            :nZr  => rand(Uniform(50.0, 120.0)), #[mm] reservoir depth/"field storage capacity of the soil" [mm] - taken from Entekhabi et al. (1992)
            :a    => rand(Uniform(11.4, 15.6)),  #numerical parameter from Bretherton et al. (2004)
            :b    => rand(Uniform(0.5, 0.6)),    #numerical parameter from Bretherton et al. (2004)
            :wsat => rand(Uniform(65.0, 80.0)),  #[mm] saturation water vapour pass derived from plots in Bretherton et al.(2004) - NEEDS MORE RESEARCH
            :τ    => rand(Uniform( (1.0 * m2mm(1)/s2day(1))/(40000.0 * km2mm(1)), (10.0 * m2mm(1)/s2day(1))/(1000.0 * km2mm(1)) )), #using u_min=1m/s, u_max=10m/s, L_min=1000km, L_max=40000km
            :pt   => 10.0, #tuning parameter for tanh function
        )
    
    else
        println("Indicate whether to use L-u parameter set (no input) or tau parameter set (true).")
    end

    d[:sfc] = d[:spwp] + 0.3 #inferred from Hagemann & Stacke (2015)

    return d
end


