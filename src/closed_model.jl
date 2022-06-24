"""
    cm_derived_quantities!(df::DataFrame)

Take a dataframe with closed model equilibrium solutions and corresponding aprameter values and
compute P1, P2, P3, El, infilt, runoff, PR (χ) and some other quantities from the equilibrium 
solutions for s, wl and wo and return them as new columns of the input DataFrame.
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

    #CSV.write(datadir("sims", df_name * "_all_quantities.csv"), df)
    return df
end

function cm_equilibrium_solution(x0, tau::Bool)
    
    p = cm_rand_params(tau) 

    if tau == true
        dynsys = ContinuousDynamicalSystem(closed_model_τ, x0, p)
    elseif tau == false
        dynsys = ContinuousDynamicalSystem(closed_model_uL, x0, p)
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

function cm_fixed_params(tau=true)

    d = Dict{Symbol, Float64}(
        :spwp => 0.3,     
        :ep   => 4.38,    
        :pt   => 10.0,         
        :eo   => 3.0,     
        :ϵ    => 1.0,     
        :r    => 2.0,     
        :α    => 0.1,     
        :nZr  => 100.0,   
        :a    => 15.6,    
        :b    => 0.603,  
        :wsat => 72.0,    
    )

    d[:sfc] = d[:spwp] + 0.3
    
    if tau == true
        d[:τ] = 0.5
    elseif tau == false
        d[:u] = 5.0 * m2mm(1.0)/s2day(1.0)
        d[:L] = 1000.0 * km2mm(1.0)
    else
        error("Choose whether to specify τ (true) or L and u (false).")
    end

    return d
end

function cm_rand_params(tau::Bool=true)

    if tau == true

        d = Dict{Symbol, Float64}(
            :spwp => rand(Uniform(0.15, 0.55)),  #permanent wilting point, taken from Hagemann & Stacke (2015)
            :ep   => rand(Uniform(4.0, 6.0)),    #[mm/day] potential evaporation over land in mm/day, taken from Entekhabi et al. (1992)
            :eo   => rand(Uniform(2.5, 3.5)),    #[mm/day] ocean evaporation rate, lower limit from Kumar et al. (2017), lower limit motivated by Zang et al. (1995)
            :ϵ    => rand(Uniform(0.9, 1.1)),    #numerical parameter from Rodriguez-Iturbe et al. (1991)
            :r    => rand(Uniform(2, 6)),         #numerical parameter from Rodriguez-Iturbe et al. (1991) and Entekhabi et al. (1992)
            :α    => rand(Uniform(0.0, 1.0)),    #land fraction
            :nZr  => rand(Uniform(50.0, 120.0)), #[mm] reservoir depth/"field storage capacity of the soil" [mm] - taken from Entekhabi et al. (1992)
            :a    => rand(Uniform(11.4, 15.6)),  #numerical parameter from Bretherton et al. (2004)
            :b    => rand(Uniform(0.5, 0.6)),    #numerical parameter from Bretherton et al. (2004)
            :wsat => rand(Uniform(65.0, 80.0)),  #[mm] saturation water vapour pass derived from plots in Bretherton et al.(2004) - NEEDS MORE RESEARCH
            :τ    => rand(Uniform( (1.0 * m2mm(1)/s2day(1))/(40000.0 * km2mm(1)), (10.0 * m2mm(1)/s2day(1))/(1000.0 * km2mm(1)) )), #using u_min=1m/s, u_max=10m/s, L_min=1000km, L_max=40000km
            :pt   => 10.0,                       #tuning parameter for tanh function
        )
        
    
    elseif tau == false

        d = Dict{Symbol, Float64}(
            :spwp => rand(Uniform(0.15, 0.55)),   #permanent wilting point
            :ep   => rand(Uniform(4.0, 6.0)),     #[mm/day] potential evaporation over land in mm/day, taken from Entekhabi et al. (1992)
            :eo   => rand(Uniform(2.5, 3.5)),     #[mm/day] ocean evaporation rate, lower limit from Kumar et al. (2017), lower limit motivated by Zang et al. (1995)
            :ϵ    => rand(Uniform(0.9, 1.1)),     #numerical parameter from Rodriguez-Iturbe et al. (1991)
            :r    => rand(Uniform(2, 6)),         #numerical parameter from Rodriguez-Iturbe et al. (1991)
            :α    => rand(Uniform(0.0, 1.0)),     #land fraction
            :nZr  => rand(Uniform(50.0, 120.0)),  #[mm] reservoir depth/"field storage capacity of the soil" [mm] - NEEDS MORE RESEARCH
            :a    => rand(Uniform(11.4, 15.6)),   #numerical parameter from Bretherton et al. (2004)
            :b    => rand(Uniform(0.5, 0.6)),     #numerical parameter from Bretherton et al. (2004)
            :wsat => rand(Uniform(65.0, 80.0)),   #[mm] saturation water vapour pass derived from plots in Bretherton et al.(2004) - NEEDS MORE RESEARCH
            :u    => rand(Uniform(1.0, 10.0)) * m2mm(1.0)/s2day(1.0), #[mm/day] wind speed
            :L    => 10000.0 * km2mm(1.0),        #[mm] domain size
            :pt   => 10.0,                        #tuning parameter for tanh function
        )
    
    else
        println("Indicate whether to use L-u parameter set (tau = false) or tau parameter set (tau = true).")
    end

    d[:sfc] = d[:spwp] + 0.3                      #inferred from Hagemann & Stacke (2015)

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






