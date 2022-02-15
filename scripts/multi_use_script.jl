using DrWatson
@quickactivate "Oceland Model"

using CairoMakie
using GLMakie
#using Colors
#using Random
using DataFrames
#using Distributions
using DifferentialEquations
using DynamicalSystems
using InteractiveDynamics
using CSV
include(srcdir("parametrisations.jl"))
include(srcdir("utils.jl"))
include(srcdir("model_versions.jl"))
include(srcdir("cm_analysis.jl"))
include(srcdir("om_analysis.jl"))
include(srcdir("cm_plotting.jl"))


fs= 20.0 #fontsize
lw= 2.0  #linewidth
ms= 5.0  #markersize

"""
    evap_tanh_plot(mode::String)

Takes two different modes of computation: 
"pt" to produce a plot of smooth evapotranspiration curves with various tuning parameter values and
"pwp" to produce a plot of smooth evapotranspiration curves with various spwp and sfc and constant pt = 10.

"""
function evap_tanh_plot(mode::String)

    fig = Figure()
    ax  = Axis(fig[1,1])
    ax.xlabel = "soil moisture saturation"
    ax.ylabel = "evaporation rate [mm/day]"
    ax.xlabelsize = fs
    ax.ylabelsize = fs
    hidespines!(ax, :t, :r)
    
    c = [:dodgerblue, :chocolate1, :green, :orangered2]
    s    = collect(0.0:0.01:1.0)
    Ep   = 4.38
    
    if mode == "pt"
        for n = 1:4
            spwp = 0.4
            sfc  = spwp + 0.3
            A    = Ep/2
            B    = 16 - 2 * n
            C    = (spwp + sfc)/2
            parameters = @dict A B C      
            
            if n == 1
                piecewise = lines!(ax, s, [land_evap(el, spwp, sfc, Ep) for el in s], linestyle = :dot, color = :gray22, linewidth = 2.5, label = "piecewise parametrisation")
            end

            new = lines!(ax, s, [evap_tanh(el, parameters) for el in s], linewidth = lw, color = c[n], label = "tuning parameter = $(B)")
            ax.title = "s_pwp = 0.4, s_fc = 0.7"
        
            save(plotsdir("Sketches", "Evap_tanh_p_tuning_varied.png"), fig)

        end
    
    elseif mode == "pwp"
        for n = 1:4
            spwp = 0.1 + n * 0.1
            sfc  = spwp + 0.3
            A    = Ep/2
            B    = 10
            C    = (spwp + sfc)/2
            parameters = @dict A B C      
      
            new = lines!(ax, s, [evap_tanh(el, parameters) for el in s], linewidth = lw, color = c[n], label = "s_pwp = $(round(spwp, digits=1)), s_fc = $(round(sfc, digits=1))")
            piecewise = lines!(ax, s, [land_evap(el, spwp, sfc, Ep) for el in s], linestyle = :dash, color = c[n], linewidth = lw)
            ax.title = "Tuning parameter = 10"
            axislegend(position = :lt, labelsize = 16, framevisible = false)
            #save(plotsdir("Sketches", "Evap_tanh_pwp_varied.png"), fig)
        end

    else
        println("Which parameter do you want to vary?")
    
    end

    return fig
end



function diffeq_issue()
    #x0 = [0.3, 50.0, 40.0]
    x0 = [0.3, 0.5, 0.4]
    #p = [2.0, 70.0, 4.0e8]
    p = Dict(
        :a => 2.0,
        :b => 70.0,
        :c => 4.0e8,
    )

    function my_system(x, p, t)
        @unpack a, b, c = p
        dx = exp(x[2] - 0.5) * (1.0 - x[1]) - a * tanh(10.0 * (x[1] - 0.5)) - a
        dy = a * tanh(10.0 * (x[1] - 0.5)) + a - exp(x[2] - 0.5) + (x[3] - x[2]) * 2.0
        dz = 3.0 - exp(x[3] - 0.5) - (x[3] - x[2]) * 0.5
        # dx = (15.0 * (exp(x[2]/b - 0.6)) * (1.0 - x[1]^2) - a * tanh(10.0 * (x[1] - 0.5)) - a)/100.0
        # dy = a * tanh(10.0 * (x[1] - 0.5)) + a - exp(15.0 * (x[2]/b - 0.6)) + (x[3] - x[2]) * 2.0
        # dz = 3.0 - exp(15.0 * (x[3]/b - 0.6)) - (x[3] - x[2]) * 0.5
        return SVector(dx, dy, dz)
    end

    #p_fixed = cm_fixed_params()
    ds = ContinuousDynamicalSystem(my_system, x0, p)
    u0s =  [ds.u0]
    diffeq = (alg = Rodas5(), adaptive = false, dt = 0.001, reltol = 1e-8, abstol = 1e-8)
    tr = trajectory(ds, 100; diffeq...)
    ps = Dict(
        :a => 1.9:0.01:2.1
        #:spwp  => 0.2:0.01:0.54,   
        #:eo => 2.5:0.1:3.5,
    )

    fig, obs = interactive_evolution_timeseries(
        ds, u0s, ps; tail = 100000, diffeq, idxs = (1, 2, 3), 
        lims =((0.0,1.0), (0.0, 2.0), (0.0, 2.0))
    )
    #prob = ODEProblem(my_system, x0, tspan, p)
    #sol = solve(prob, alg = Rodas5(), reltol=1e-9, abstol=1e-9, saveat=1.0)
    return fig
end

#td = cm_rand_params()

function hm_plot()
    sr = range(0.0, 1.0, length = 100)
    wr = range(40.0, 60.0, length = 100)
    blc = [ exp(15.0 * (w / 72.0 - 0.58)) * (1 - s^2) - (4.3/2 * tanh(10 * (s - (0.3+0.6)/2 ) ) + 4.3/2) for s in sr, w in wr]
    fig = Figure()
    ax1  = Axis(fig[1,1])
    hm = heatmap!(ax1, sr, wr, blc)
    Colorbar(fig[1,2], hm)
    return fig
end

function line_plot()
    α = collect(0.0:0.001:1.0)
    s = collect(0.0:0.01:1.0)
    ep = 4.3
    pt = 10.0
    spwp = 0.3
    sfc = 0.6
    u = 4.32e8
    L1 = 1e9
    L2 = 1e10

    fig = Figure()
    ax  = Axis(fig[1,1])
    hidespines!(ax, :t, :r)
    #lines!(ax, s, (ep/2 .* tanh.( pt .* (s .- (spwp.+sfc)./2 ) ) .+ ep./2) ./ (1 .- s.^2), color = (:royalblue, 0.6), linewidth = 3.0, label = "Pl = El / Φ")
    lines!(ax, α, u ./ (α .* L1), color = :darkgreen, linewidth = 3.0, label = "L = 1000km")
    lines!(ax, α, u ./ (α .* L2), color = :chartreuse4, linewidth = 3.0, label = "L = 10000km")
    lines!(ax, α, u ./ ((1 .- α) .* L1), color = :dodgerblue4, linewidth = 3.0, label = "L = 1000km")
    lines!(ax, α, u ./ ((1 .- α) .* L2), color = :dodgerblue, linewidth = 3.0, label = "L = 10000km")
    hlines!(ax, 1.0, linestyle = :dot, color = :grey)
    ax.xlabel = "land fraction α"
    ax.ylabel = "u/(αL) in green and u/((1-α)L) in blue"
    ylims!(low = 0, high=15.0)
    axislegend(ax, position = :lt, frame_visible = false)
    save(plotsdir("Closed model", "u-α-L.png"),fig)
end


function cm_load_data()
    s1000 = CSV.read(datadir("sims", "closed model pmscan/cm_smooth_eq_MC_fixedpoints_10000_runs_domain1000_all_quantities" * ".csv"), DataFrame)
    s10000 = CSV.read(datadir("sims", "closed model pmscan/cm_smooth_eq_MC_fixedpoints_10000_runs_domain10000_all_quantities" * ".csv"), DataFrame)
    return s1000, s10000
end

function cm_tau_load_data()
    d = CSV.read(datadir("sims", "closed model pmscan/cm_smooth_tau_eq_MC_fixedpoints_runs10000_all_quantities" * ".csv"), DataFrame)
    return d
end

function om_load_data()
    s1000 = CSV.read(datadir("sims", "open model pmscan/om_v2_MC_fixedpoints_runs10000_domain1000_all_quantities" * ".csv"), DataFrame)
    #s10000 = CSV.read(datadir("sims", "open model pmscan/om_v2_closed_MC_fixedpoints_runs10000_domain10000_all_quantities" * ".csv"), DataFrame)
    s10000 = CSV.read(datadir("sims", "open model pmscan/om_v2_MC_fixedpoints_runs10000_domain10000_all_quantities" * ".csv"), DataFrame)

    return s1000, s10000
end

function load_lin_data()
    cd = CSV.read(datadir("sims", "closed model pmscan/cm_smooth_linearised_eq_MC_fixedpoints_runs5000_domain10000_all_quantities" * ".csv"), DataFrame)
    od = CSV.read(datadir("sims", "open model pmscan/om_v2_closed_lin_MC_fixedpoints_runs5000_domain10000_all_quantities" * ".csv"), DataFrame)
    return cd, od
end

function load_steprange()
    d = CSV.read(datadir("sims", "closed model pmscan/cm_smooth_Lvaried_eq_MC_fixedpoints_runs10000_all_quantities" * ".csv"), DataFrame)
    return d 
end


function p_break()
        d = Dict{Symbol, Float64}(
        :spwp => 0.2205204130513438,  #permanent wilting point
        :ep   => 4.179071665580645,   #[mm/day] potential evaporation over land in mm/day, taken from [1]
        :eo   => 3.063503944375833,   #[mm/day] ocean evaporation rate
        :ϵ    => 0.9976897025816163,  #numerical parameter from Rodriguez-Iturbe et al. (1991)
        :r    => 2.0,                 #numerical parameter from Rodriguez-Iturbe et al. (1991)
        :α    => 0.9925995771678071,  #land fraction
        :nZr  => 103.09181308321128,  #[mm] reservoir depth/"field storage capacity of the soil" [mm] - NEEDS MORE RESEARCH
        :a    => 11.714187532725438,  #numerical parameter from Bretherton et al. (2004)
        :b    => 0.5642281161109617,  #numerical parameter from Bretherton et al. (2004)
        :wsat => 66.73351939736625,   #[mm] saturation water vapour pass derived from plots in Bretherton et al.(2004) - NEEDS MORE RESEARCH
        :u    => 7.427281529064258e8, #[mm/day] wind speed
        :L    => 1000.0 * km2mm(1.0), #rand(Uniform(500.0, 5000.0)) * km2mm(1.0), #[mm] domain size
        :pt   => 10.0,                #tuning parameter for tanh function
        :L2   => 9.925995771678071e8,
        :L1   => 3.700211416096449e6,
        :L3   => 3.700211416096449e6,
        :sfc  => 0.5205204130513438,
        :w0   => 1.1089681429099199,        
    )
    return d
end


function parametrisation_plots()
    s = collect(0.0:0.01:1.0)
    w = collect(0.0:0.1:60.0)
    ep = 4.3
    pt = 10.0
    spwp = 0.3
    sfc = 0.6
    ϵ = 1.0
    r = 2.0
    a = 15.6
    b = 0.603
    wsat = 72.0
    E = [(ep/2 * tanh( pt * (el - (spwp+sfc)/2 ) ) + ep/2) for el in s]
    R  = [(ϵ * el^r) for el in s]
    P  = [exp(a*(el / wsat - b)) for el in w]

    f1 = Figure()
    ax1 = Axis(f1[1,1], xlabel = "Soil moisture saturation s", ylabel = "Evapotranspiration [mm/day]", xlabelsize = 28, ylabelsize = 28)
    lines!(ax1, s, E, linewidth = 3.0, color = :dodgerblue4)
    hidespines!(ax1, :t, :r)

    f2 = Figure()
    ax2 = Axis(f2[1,1], xlabel = "Soil moisture saturation s", ylabel = "Runoff fraction R_f", xlabelsize = 28, ylabelsize = 28)
    lines!(ax2, s, R, linewidth = 3.0, color = :dodgerblue4)
    hidespines!(ax2, :t, :r)

    f3 = Figure()
    ax3 = Axis(f3[1,1], xlabel = "Column water vapour pass w [mm]", ylabel = "Precipitation [mm/day]", xlabelsize = 28, ylabelsize = 28)
    lines!(ax3, w, P, linewidth = 3.0, color = :dodgerblue4)
    hidespines!(ax3, :t, :r)

    return f1, f2, f3
end


function smile(df1, df2)

    f = Figure(resolution = (1000, 800))
    #ax= Axis(f[1,1], xlabel = "Land fraction α", ylabel = "Precipitation ratio PR = Pl/Po", xlabelsize = 28, ylabelsize = 28, yticklabelsize = 20, xticklabelsize = 20)
    #scatter!(ax, df1.α, df1.PR, markersize = 10.0, color = :dodgerblue3, label = "closed model")
    #scatter!(ax, df2.α, df2.PR1, markersize = 10.0, color = :chartreuse4, label = "open model, closed version")
    #lines!(ax, sort(df1, "α")[!,"α"], cm_mean_of_bins!(df1,"α","PR",100), linewidth = 10.0, color = :midnightblue)
    #lines!(ax, sort(df2, "α")[!,"α"], cm_mean_of_bins!(df2,"α","PR1",100), linewidth = 10.0, color = :darkgreen)

    hlines!(ax, 1.0, color = :grey, linestyle = :dot)
    axislegend(ax, position = :lb, labelsize = 20)
    hidespines!(ax, :t, :r)
    hidedecorations!(ax, ticks = false, ticklabels = false, label = false)
    return f

end


function arbplot2(df1, df2)

    f = Figure(resolution = (1000, 800))
    ax= Axis(f[1,1], xlabel = "Land fraction α", ylabel = "Precipitation ratio PR = P2/P1", xlabelsize = 28, ylabelsize = 28, yticklabelsize = 20, xticklabelsize = 20, title = "w0 < 40mm", titlesize = 30)
    scatter!(ax, df1.α, df1.P2divP1, markersize = 10.0, color = (:dodgerblue3, 0.5), label = "L = 1000 km")
    scatter!(ax, df2.α, df2.P2divP1, markersize = 10.0, color = (:chartreuse4, 0.5), label = "L = 10000 km")
    #lines!(ax, sort(df1, "α")[!,"α"], cm_mean_of_bins!(df1,"α","P2divP1",200), linewidth = 8.0, color = :midnightblue)
    #lines!(ax, sort(df2, "α")[!,"α"], cm_mean_of_bins!(df2,"α","P2divP1",199), linewidth = 8.0, color = :darkgreen)

    hlines!(ax, 1.0, color = :grey, linestyle = :dot)
    axislegend(ax, position = :lb, labelsize = 20, framevisible = false)
    hidespines!(ax, :t, :r)
    hidedecorations!(ax, ticks = false, ticklabels = false, label = false)
    return f

end


function arbplot(df, x::String, y::String)

    f = Figure(resolution = (1000, 900))
    ax= Axis(f[1,1], xlabel = x, ylabel = y, xlabelsize = 38, ylabelsize = 38, yticklabelsize = 30, xticklabelsize = 30)
    scatter!(ax, df[!,x], df[!,y], markersize = 8.0, color = :chartreuse4) 
    lines!(ax, sort(df, x)[!,x], cm_mean_of_bins!(df,x,y,100), color = :darkgreen, linewidth = 4.0)
    #scatter!(ax, df.α, df[!,"PR"], markersize = 8.0, color = (:chocolate, 0.8), label = "Pl/Po_mean") #:chartreuse4
    #scatter!(ax, df.α, df[!,"PR2"], markersize = 8.0, color = (:dodgerblue3, 0.8), label = "Pl/P(wo_mean)")
    #scatter!(ax, df.α, df[!,"P3"], markersize = 8.0, color = (:chartreuse4, 0.8), label = "ocean 2")
    #lines!(ax, sort(df, "α")[!,"α"], cm_mean_of_bins!(df,"α","PR",100), color = :chocolate4, linewidth = 4.0)
    #lines!(ax, sort(df, "α")[!,"α"], cm_mean_of_bins!(df,"α","PR2",100), color = :navyblue, linewidth = 4.0)
    #lines!(ax, sort(df, "α")[!,"α"], cm_mean_of_bins!(df,"α","P3",100), color = :darkgreen, linewidth = 4.0)
    #axislegend(ax, position = :lb, labelsize = 38, framevisible = false)
    hidespines!(ax, :t, :r)
    hidedecorations!(ax, ticks = false, ticklabels = false, label = false)
    return f

end

    
function no_influx_ocean_plot()
    wsat = 72.0
    a = 15.6
    b = 0.603
    w = collect(0.0:1.0:72.0)
    eo = 3.0
    influx = [eo for elm in w]
    τ = [0.0864, 0.2160, 0.8640]
    lfs = 16
    outflux1 = [(exp.(a .* (elm ./ wsat .- b)) .* τ[1]) for elm in w]
    fig = Figure(resolution = (900, 500))
    ax = Axis(fig[1,1], xlabel = "Water vapor pass [mm]", ylabel = "Water flux rates [mm/day]", ylabelsize = lfs, xlabelsize = lfs)
    hidespines!(ax, :t, :r)
    lines!(ax, w, influx, label = "evaporative influx")
    
    for i=1:3
        outflux = [(exp.(a.*(elm ./ wsat .- b)).+ elm .* τ[i]) for elm in w]
        lines!(ax, w, outflux, label = "outflux for τ = $(τ[i]) per day")
    end
    xlims!(ax, (0, 40))
    ylims!(ax, (0, 20))
    axislegend(ax, position = :lt, framevisible = false, labelsize = 16)
    return fig
end

function mean_to_total(d::DataFrame)
    df = deepcopy(d)

    #Total advected water (needs to be computed with mean value!)
    df.A01 = df.w0 .* df.u
    df.A12 = df.w1 .* df.u
    df.A23 = df.w2 .* df.u
    df.A34 = df.w3 .* df.u
    df.A1 = df.A01 .- df.A12
    df.A2 = df.A12 .- df.A23
    df.A3 = df.A23 .- df.A34
    df.Atot = df.A01 .- df.A34


    #Total water content
    df.w1 = df.w1 .* df.L1
    df.w2 = df.w2 .* df.L2
    df.w3 = df.w3 .* df.L3
    df.s  = df.s .* df.L2 .* df.nZr

    #Total precipitation
    df.P1 = df.P1 .* df.L1
    df.P2 = df.P2 .* df.L2
    df.P3 = df.P3 .* df.L3
    
    #Total evapo(transpi)ration
    df.E1 = df.eo .* df.L1
    df.E2 = df.El .* df.L2
    df.E3 = df.eo .* df.L3

    #Total runoff
    df.R = df.R .* df.L2

    #Total precip ratio
    df.PRtot = df.P2 ./ (df.P1 .+ df.P3)

    return df
end