"""
    cm_derived_quantities!(df::DataFrame)

Take a dataframe with closed model equilibrium solutions and corresponding aprameter values and
compute P1, P2, P3, El, infilt, runoff, PR (χ) and some other quantities from the equilibrium 
solutions for s, wl and wo and return them as new columns of the input DataFrame.
"""
function cm_derived_quantities!(df_name)

    df = CSV.read(datadir("closed model", df_name * ".csv"), DataFrame)
    
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

function cm_lin_derived_quantities!(df_name)
    df = CSV.read(datadir("closed model", df_name * ".csv"), DataFrame)

    df.Pl = df.p .* df.wl
    df.Po = df.p .* df.wo
    df.χ  = df.Pl ./ df.Po
    df.El = df.e .* df.s
    df.Rf = df.r .* df.s
    df.Φ  = 1 .- df.Rf
    df.R  = df.Rf .* df.Pl
    df.Al = (df.wo .- df.wl) .* df.τ ./ df.α
    df.Ao = (df.wo .- df.wl) .* df.τ ./ (1 .- df.α)
    df.B  = df.Ao
    df.sblc = (df.Pl .* df.Φ .- df.El) ./ df.nZr
    df.lblc = df.El .- df.Pl .+ df.Al
    df.oblc = df.eo .- df.Po .- df.Ao

    CSV.write(datadir("closed model", df_name * "_all_quantities.csv"), df)
    return df
end
    
function cm_DC_EQmeanvalues(ds::GeneralizedDynamicalSystem, x0, p)
    sg = range(0, 1; length = 100)
    wlg = wog = range(0, p[:wsat], length = 100)
    grid = (sg, wlg, wog)
    mapper = AttractorsViaRecurrences(ds, grid)
    mapper(x0)
    attrs = mapper.bsn_nfo.attractors
    if length(attrs) > 1
        @warn "More than one attractor found. Now you have to go implement a solution for this case."
    elseif length(attrs) == 0
        @warn "No attractor found."
        return nothing
    else 
        seq = round(mean(attrs[1][:,1]), digits=2)
        wleq = round(mean(attrs[1][:,2]), digits=1)
        woeq = round(mean(attrs[1][:,3]), digits=1)
        eq_mean_vals = [seq, wleq, woeq]
        return eq_mean_vals
    end
end


function cm_DC_derived_quantities(data::Dataset{4, Float64}, p)
    ## Naming for when solution to solve() was passed instead of a Dataset
    #t  = sol.t
    #s  = sol[1,:]
    #wl = sol[2,:]
    #wo = sol[3,:]

    t  = data[:,1]
    s  = data[:,2]
    wl = data[:,3]
    wo = data[:,4]
    wind = [wind_DC(elm, p) for elm in t]
    ad_mois = [advected_moisture(elm[1], elm[2], elm[3], p) for elm in zip(wl, wo, t)]
    land_size = p[:α] * p[:L]
    ocean_size = (1 .- p[:α]) * p[:L]
    El = [evap_scaling(elm, p) for elm in t] .* [El_tanh(elm, p) for elm in s]
    Eo = [p[:eo] for elm in t]
    Pl = [precip(elm, p) for elm in wl]
    Po = [precip(elm, p) for elm in wo]
    PR = Pl./Po
    PRmean = movingaverage(PR, 100)
    Φ = [infiltration(elm,p) for elm in s]
    R  = (1 .- Φ) .* Pl
    Al = 2 .* ad_mois .* wind ./ land_size
    Ao = - 2 .* ad_mois .* wind ./ ocean_size
    balance = Eo .+ El .- Po .- Pl
    df = DataFrame(t = t, s = s, wl = wl, wo = wo, wind = wind, El = El, Eo = Eo, Pl = Pl, Po = Po, PR = PR, PRmean = PRmean, Φ = Φ, R = R, Al = Al, Ao = Ao, watbal = balance)
    return df
end


function cm_equilibrium_solution(x0, tau::Bool)
    
    p = cm_rand_params(tau) 

    if tau == true
        #dynsys = ContinuousDynamicalSystem(cm_τ, x0, p)
        dynsys = ContinuousDynamicalSystem(cm_τ_Pl_neq_Po, x0, p)
    elseif tau == false
        dynsys = ContinuousDynamicalSystem(cm_uL, x0, p)
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


function cm_fixed_params(tau::Bool=true)

    d = Dict{Symbol, Any}(
        :spwp => 0.5,     
        :ep   => 4.38,    
        :pt   => 10.0,         
        :eo   => 3.0,     
        :ϵ    => 1.0,     
        :r    => 2.0,     
        :α    => 0.3,     
        :nZr  => 100.0,   
        :a    => 15.6,    
        :b    => 0.603,  
        :wsat => 72.0,    
        :λ    => 1500.00, #moisture scale height
        #:P0   => hPa2Pa(1000.00), #surface pressure
        :Rd   => 287, #specific gas constant for dry asymmetric
        :cα   => 4*10^(-3), #coefficient in air resistance term for land
        :rsfc => 0.0, #soil resistance at field capacity - NEEDS MORE RESEARCH
        :T_amp  => 5.0, #diurnal surface temperature amplitude
        :T_mean => 300, #diurnal mean surface temperature
        :t_shift=> 1/24, #time-shift of wind wrt surface temperature measured in days
        :u_max  => m2mm(1.0)/s2day(1.0), #maximum wind speed during the course of a day
        :f_a  => 0.0,
    )

    d[:sfc] = d[:spwp] + 0.3
    
    if tau == true
        d[:τ] = 0.5
    elseif tau == false
        d[:u] = 5.0 * m2mm(1.0)/s2day(1.0)
        d[:L] = 400.0 * km2mm(1.0)
    else
        error("Choose whether to specify τ (true) or L and u (false).")
    end

    #Create land-sea mask
    d[:dx] = 10.0 * km2mm(1.0)
    d[:nb_boxes] = Int(round(d[:L] / d[:dx]))
    if d[:L] / d[:dx] % 1 != 0 
        d[:L] = d[:nb_boxes] * d[:dx]
        @warn "Number of boxes was not an integer. Domain got adapted accordingly."
    end
    nb_ocean_boxes = Int(round((1 - d[:α]) * d[:L] / d[:dx]))
    nb_land_boxes = d[:nb_boxes] - nb_ocean_boxes
    ((1 - d[:α]) * d[:L] / d[:dx]) % 1 != 0 && @warn "Number of boxes did not perfectly match the land fraction. Ocean was enlarged and land fraction adapted accordingly."
    d[:α] = nb_land_boxes * d[:dx] / d[:L]
    d[:lsmask] = vcat(zeros(nb_ocean_boxes), ones(nb_land_boxes))

    # Create x-Vector
    # Grid points at interval centers
    x_centers = [
        d[:dx]/2 + d[:dx] * i
        for i = 0:d[:nb_boxes]-1
    ]

    # Grid points at left interval boundaries
    x_left = [
        d[:dx] * i for i = 0:d[:nb_boxes]-1
    ]

    d[:x] = x_centers

    return d
end



function cm_rand_params(tau::Bool=true)

    if tau == true

        d = Dict{Symbol, Any}(
            :spwp => rand(Uniform(0.15, 0.55)),  #permanent wilting point, taken from Hagemann & Stacke (2015)
            :ep   => rand(Uniform(4.0, 6.0)),    #[mm/day] potential evaporation over land in mm/day, taken from Entekhabi et al. (1992)
            :eo   => rand(Uniform(2.5, 3.5)),    #[mm/day] ocean evaporation rate, lower limit from Kumar et al. (2017), lower limit motivated by Zang et al. (1995)
            :ϵ    => rand(Uniform(0.9, 1.1)),    #numerical parameter from Rodriguez-Iturbe et al. (1991)
            :r    => rand(Uniform(2, 6)),         #numerical parameter from Rodriguez-Iturbe et al. (1991) and Entekhabi et al. (1992)
            :α    => rand(Uniform(0.0, 1.0)),    #land fraction
            :nZr  => rand(Uniform(50.0, 120.0)), #[mm] reservoir depth/"field storage capacity of the soil" [mm] - taken from Entekhabi et al. (1992)
            :a    => rand(Uniform(11.4, 15.6)),  #numerical parameter from Bretherton et al. (2004)
            :b    => rand(Uniform(0.5, 0.6)),    #numerical parameter from Bretherton et al. (2004)
            :bland=> 0.999*0.6,
            :boce => 0.6,  
            :wsat => rand(Uniform(65.0, 80.0)),  #[mm] saturation water vapour pass derived from plots in Bretherton et al.(2004) - NEEDS MORE RESEARCH
            :τ    => rand(Uniform( (1.0 * m2mm(1)/s2day(1))/(40000.0 * km2mm(1)), (10.0 * m2mm(1)/s2day(1))/(1000.0 * km2mm(1)) )), #using u_min=1m/s, u_max=10m/s, L_min=1000km, L_max=40000km
            :pt   => 10.0,                       #tuning parameter for tanh function
        )
        
    
    elseif tau == false

        d = Dict{Symbol, Any}(
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

function cm_lin_rand_params()
    d = Dict{Symbol, Float64}(
            :p    => rand(Uniform(0.06, 0.08)),#(0.07, 0.0875)),#0.025, 0.06)), #0.3, 0.5)),  #slope of precipitation as function of w
            :e    => rand(Uniform(4.0, 6.5)),#(4.0, 6.75)),    #slope of land evaporation as function of s
            :eo   => rand(Uniform(2.5, 3.5)),    #[mm/day] ocean evaporation rate, lower limit from Kumar et al. (2017), lower limit motivated by Zang et al. (1995)
            :r    => rand(Uniform(0.9, 1.0)),        #slope of runoff fraction as function of s
            :α    => rand(Uniform(0.0, 1.0)),    #land fraction
            :nZr  => rand(Uniform(50.0, 120.0)), #[mm] reservoir depth/"field storage capacity of the soil" [mm] - taken from Entekhabi et al. (1992)
            :τ    => rand(Uniform( (1.0 * m2mm(1)/s2day(1))/(40000.0 * km2mm(1)), (10.0 * m2mm(1)/s2day(1))/(1000.0 * km2mm(1)) )), #using u_min=1m/s, u_max=10m/s, L_min=1000km, L_max=40000km
        )
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
    dfclean.ds = (dfclean.P2 .* dfclean.Φ .- dfclean.El) ./ dfclean.nZr
    dfclean.dw2 = dfclean.El .- dfclean.P2 .+ (dfclean.w1 .- dfclean.w2) .* dfclean.u ./ dfclean.L2
    dfclean.dw3 = dfclean.eo .- dfclean.P3 .+ (dfclean.w2 .- dfclean.w3) .* dfclean.u ./ dfclean.L3
    dfclean.τ  = dfclean.u ./ dfclean.L
    dfclean.PR12 = dfclean.P2 ./ dfclean.P1
    dfclean.Δw1 = dfclean.w0 .- dfclean.w1
    dfclean.Δw2 = dfclean.w1 .- dfclean.w2
    dfclean.Δw3 = dfclean.w2 .- dfclean.w3
    dfclean.Δwtot = dfclean.w0 .- dfclean.w3


    if occursin(r"v2_closed", df_name) == true
        dfclean = select!(dfclean, Not(:w0))
        dfclean.wo = (dfclean.L1 .* dfclean.w1 + dfclean.L3 .* dfclean.w3) ./ (dfclean.L1 .+ dfclean.L3)
        dfclean.PR2 = dfclean.P2 ./ (exp.(dfclean.a .* (dfclean.wo ./ dfclean.wsat .- dfclean.b)))
        dfclean.dw1 = dfclean.eo .- dfclean.P1 .+ (dfclean.w3 .- dfclean.w1) .* dfclean.u ./ dfclean.L1
    elseif occursin(r"v2_closed", df_name) == false
        dfclean.dw1 = dfclean.eo .- dfclean.P1 .+ (dfclean.w0 .- dfclean.w1) .* dfclean.u ./ dfclean.L1
    end

    #convert to better units
    dfclean.L_km = dfclean.L .* mm2km(1.0)
    dfclean.L1_km = dfclean.L1 .* mm2km(1.0)
    dfclean.L2_km = dfclean.L2 .* mm2km(1.0)
    dfclean.L3_km = dfclean.L3 .* mm2km(1.0)
    dfclean.u_ms = dfclean.u .* (mm2m(1.0) ./ day2s(1.0))
    
    #CSV.write(datadir("sims", df_name * "_all_quantities.csv"), dfclean)
    return dfclean
end


function om_equilibrium_solution(x0, model_version::String)

    p = om_rand_params()

    if model_version == "paper"
        dynsys = ContinuousDynamicalSystem(om_v_paper, x0, p)
    elseif model_version == "closed" 
        dynsys = ContinuousDynamicalSystem(om_closed, x0, p)
    elseif model_version == "w_tracked"
        dynsys = ContinuousDynamicalSystem(om_w_tracked, x0, p)
    else
        error("Assign valid model_version name: paper, closed or w_tracked.")
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


function om_rand_params()

    d = Dict{Symbol, Float64}(
        :spwp => rand(Uniform(0.20, 0.55)), #permanent wilting point
        :ep   => rand(Uniform(2.0, 6.0)),   #[mm/day] potential evaporation over land in mm/day, taken from [1]
        :eo   => rand(Uniform(2.8, 3.2)),   #[mm/day] ocean evaporation rate
        :ϵ    => rand(Uniform(0.9, 1.1)),   #numerical parameter from Rodriguez-Iturbe et al. (1991)
        :r    => rand(Uniform(2,6)),         #numerical parameter from Rodriguez-Iturbe et al. (1991) and Entekhabi et al. (1992)
        :α    => rand(Uniform(0.0, 1.0)),   #land fraction
        :nZr  => rand(Uniform(50.0, 120.0)),#[mm] reservoir depth/"field storage capacity of the soil" [mm] - taken from Entekhabi et al. (1992)
        :a    => rand(Uniform(11.4, 15.6)), #numerical parameter from Bretherton et al. (2004)
        :b    => rand(Uniform(0.5, 0.6)),   #numerical parameter from Bretherton et al. (2004)
        :wsat => rand(Uniform(65.0, 80.0)), #[mm] saturation water vapour pass derived from plots in Bretherton et al.(2004) - NEEDS MORE RESEARCH
        :u    => rand(Uniform(0.0, 10.0)) * m2mm(1.0)/s2day(1.0), #[mm/day] wind speed
        :L    => rand(Uniform(200.0, 2000.0)) * km2mm(1.0), #[mm] domain size
        :pt   => 10.0,                      #tuning parameter for tanh function
    )

    d[:sfc] = d[:spwp] + 0.3
    d[:w0]  = rand(Uniform(0.0, d[:wsat]))  # water vapour pass at windward model boudnary
    d[:L2]  = d[:α] * d[:L]                 # island length
    d[:L1]  = d[:L]/2 - d[:L2]/2            # length of first ocean [mm], symmetric configuration
    #d[:L1] = (d[:L] - d[:L2])/4             # length of first ocean [mm], asymmetric configuration with 3*L1 = L3
    #d[:L1] = (d[:L] - d[:L2])*3/4           # length of first ocean [mm], asymmetric configuration with L1 = 3*L3
    d[:L3]  = d[:L] - d[:L2] - d[:L1]       # length of second ocean [mm]
    d[:τ]   = d[:u]/d[:L]                   # rate of atmospheric transport

    return d
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






