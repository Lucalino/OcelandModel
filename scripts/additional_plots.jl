using DrWatson
@quickactivate "Oceland Model"

using CairoMakie
using Colors 
using DataFrames
include(srcdir("parametrisations.jl"))


function El_tuning_param_plot()

    fig = Figure()
    ax  = Axis(fig[1,1])
    ax.xlabel = "soil moisture saturation"
    ax.ylabel = "evaporation rate [mm/day]"
    ax.xlabelsize = fs
    ax.ylabelsize = fs
    hidespines!(ax, :t, :r)
    
    c = [:dodgerblue, :chocolate1, :green, :orangered2]
    s    = collect(0.0:0.01:1.0)
    ep   = 4.38
    spwp = 0.4
    sfc  = spwp + 0.3
    p = @dict ep spwp sfc
    
    for n = 1:4
        p[:pt] = 16 - 2 * n     
        if n == 1
            lines!(ax, s, [El_piecewise(el, p) for el in s], linestyle = :dot, color = :gray22, linewidth = 2.5, label = "piecewise parametrisation")
        end
        lines!(ax, s, [El_tanh(el, p) for el in s], linewidth = lw, color = c[n], label = "tuning parameter = $(p[:pt])")
        ax.title = "s_pwp = " * string(p[:spwp]) * ", s_fc = " * string(p[:sfc]) * ", ep = " * string(p[:ep])
    end
    axislegend(ax, position = :lt, framevisible = false)
    #save(plotsdir("Sketches", "Evap_tanh_p_tuning_varied.png"), fig)
    return fig

end

function advection_efficiencies_plot()
    α = collect(0.0:0.001:1.0)
    u = 4.32e8
    L1 = 1e9
    L2 = 1e10

    fig = Figure()
    ax  = Axis(fig[1,1])
    hidespines!(ax, :t, :r)
    lines!(ax, α, u ./ (α .* L1), color = :darkgreen, linewidth = 3.0, label = "L = 1000km")
    lines!(ax, α, u ./ (α .* L2), color = :chartreuse4, linewidth = 3.0, label = "L = 10000km")
    lines!(ax, α, u ./ ((1 .- α) .* L1), color = :dodgerblue4, linewidth = 3.0, label = "L = 1000km")
    lines!(ax, α, u ./ ((1 .- α) .* L2), color = :dodgerblue, linewidth = 3.0, label = "L = 10000km")
    hlines!(ax, 1.0, linestyle = :dot, color = :grey)
    ax.xlabel = "land fraction α"
    ax.ylabel = "u/(αL) in green and u/((1-α)L) in blue"
    ylims!(low = 0, high=15.0)
    axislegend(ax, position = :ct, framevisible = false)
    #save(plotsdir("Closed model", "u-α-L.png"),fig)
    return fig
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
