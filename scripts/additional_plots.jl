using DrWatson
@quickactivate "Oceland Model"

using CairoMakie
using Colors 
using DataFrames
using LaTeXStrings
include(srcdir("parametrizations.jl"))
include(srcdir("figure_labels.jl"))
include(srcdir("create_model_output.jl"))
include(srcdir("utils.jl"))


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

function max_wind_plot(days::Float64)
    l = full_labels_dict()
    lfs = 20
    lw = 3.0
    p = cm_fixed_params(false)
    t = collect(0.0:0.01:days)
    u = [diurnal_wind(el, p) for el in t]
    fig = Figure(resolution = (500,400))
    ax = Axis(fig[1,1], xlabel = l["t"], ylabel = l["umax"], ylabelsize = lfs, xlabelsize = lfs, xgridcolor = :white, ygridcolor = :white)
    lines!(ax, t, u, linewidth = lw)
    hidespines!(ax, :t, :r)
    return fig
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

function cm_t_evolution_plot(s, wl, wo, t)
    c1 = "dodgerblue"
    c2 = "darkblue"
    c3 = "darkgreen"

    fig = Figure(resolution = (1000, 600))
    ax1 = Axis(fig[1,1], xlabelsize = 28, ylabel = "w [mm]", ylabelsize = 28, yticklabelsize = 20, xticklabelsize = 20)
    ax2 = Axis(fig[2,1], xlabel = "t [days]", xlabelsize = 28, ylabel = "s", ylabelsize = 28, yticklabelsize = 20, xticklabelsize = 20)
    hidespines!(ax1, :t, :r)
    hidespines!(ax2, :t, :r)
    lines!(ax1, t, wl, linewidth = 3.0, label="w_l", color=c1)
    lines!(ax1, t, wo, linewidth = 3.0, label="w_o", color=c2)
    lines!(ax2, t, s, linewidth = 3.0, color=c3)
    axislegend(ax1, framevisible = false, labelsize = 28)

    return fig
end

function PR_α_comparison_plot(df1, df2, lb1, lb2, x, y, xlab, ylab, ttl, comp)
    fig = Figure()
    ax  = Axis(fig[1,1], xlabel = xlab, ylabel = ylab, title = ttl)
    hidespines!(ax, :t, :r)
    scatter!(ax, df1[!,x], df1[!,y], color = (:royalblue, 0.6), label = lb1, markersize = 5.0)
    scatter!(ax, df2[!,x], df2[!,y], color = (:chocolate, 0.6), label = lb2, markersize = 5.0)
    vlines!(ax, df1[df1.PR .== minimum(df1.PR), "α"][1], linewidth = 3.0, linestyle = :dash, color = (:royalblue, 0.6))
    vlines!(ax, df2[df2.PR .== minimum(df2.PR), "α"][1], linewidth = 3.0, linestyle = :dash, color = (:chocolate, 0.6))
    axislegend(ax, position = :lb)
    save(plotsdir("Closed model", "Parameter scatter plots/", "smooth", "PR-alpha_parameter_influences", "PR-alpha_$(comp)_varied.png"),fig)
    return fig
end

function rel_mi_plot(rel_mi_data::DataFrame)
    sort!(rel_mi_data, "MI_rel", rev = true)
    l = short_labels_dict()
    n = nrow(rel_mi_data)
    xrange = collect(1:1:n)
    f = Figure(resolution = (600, 450))
    ax = Axis(f[1,1], yscale = log10, xlabel = L"Model parameters $p_i$", ylabel = L"Mutual information index $I_{MI}\, (p_i\, ,\chi)$", xlabelsize = 20, ylabelsize = 20, xgridcolor = :white, ygridcolor = :white)
    #ylims!(ax, 0.1, 500.0)
    hlines!(ax, 1.0, linestyle = :dash, color = :grey35, label = "3σ significance threshold")
    scatter!(ax, xrange, rel_mi_data[:, "MI_rel"], marker = '*', markersize = 25, color = :grey20)
    ax.xticks = (1:1:n, rel_mi_data[:,"pnames"])
    axislegend(ax, framevisible = false, labelsize = 16)
    hidespines!(ax, :t, :r)
    return f
end

function runoff_plot()
    l = full_labels_dict()
    lfs = 20
    legendfs = 16
    lw = 3.0
    c = [:grey20, :grey40, :grey60, :grey70]
    s = collect(0.0:0.01:1.0)
    r = [2, 4, 6]
    labls = [L"r = 2", L"r = 4", L"r = 6"]
    fig = Figure(resolution = (500, 400))
    ax = Axis(fig[1,1], xlabel = l["s"], ylabel = l["R"], ylabelsize = lfs, xlabelsize = lfs, xgridcolor = :white, ygridcolor = :white)

    for n=1:length(r)
        R = [el^r[n] for el in s]
        lines!(ax, s, R, linewidth = lw, color = c[n], label = labls[n])
    end
    hidespines!(ax, :t, :r)
    axislegend(ax, position = :lt, framevisible = false, labelsize = legendfs)
    return fig
end

function surface_temp_plot()
    l = full_labels_dict()
    lfs = 20
    lw = 3.0
    p = cm_fixed_params(false)
    t = collect(0.0:0.01:1.0)
    Tsrf = [surface_temperature(el, p) for el in t]
    fig = Figure(resolution = (500,400))
    ax = Axis(fig[1,1], xlabel = l["t"], ylabel = l["Tsrf"], ylabelsize = lfs, xlabelsize = lfs, xgridcolor = :white, ygridcolor = :white)
    lines!(ax, t, Tsrf, linewidth = lw)
    hidespines!(ax, :t, :r)
    return fig
end
    

function time_series(x::String = "t", y::String)
    l = full_labels_dict()
    lfs = 20
    lw = 3.0
    fig = Figure(resolution = (500,400))
    ax = Axis(fig[1,1], xlabel = l[x], ylabel = l[y], ylabelsize = lfs, xlabelsize = lfs, xgridcolor = :white, ygridcolor = :white)
    lines!(ax, t, Tsrf, linewidth = lw)
    hidespines!(ax, :t, :r)
    return fig
end
