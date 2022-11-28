using DrWatson
@quickactivate "Oceland Model"

using CairoMakie
using ColorSchemes
using DataFrames
using LaTeXStrings
using DynamicalSystems
using VegaLite
using VegaDatasets
include(srcdir("parametrizations.jl"))
include(srcdir("figure_labels.jl"))
include(srcdir("create_model_output.jl"))
include(srcdir("utils.jl"))

dcm      = CSV.read(datadir("closed model/cm_tau_fixedpoints_100000_runs_all_quantities" * ".csv"), DataFrame)

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

function fig_five_threeregimes(data::DataFrame = dcm, yquant::String = "PR", xquant1::String = "τ", xquant2::String = "α", xquant3::String = "spwp", xquant4::String = "r")
    l = labels_dict()
    lfs = 20
    ms = 5
    lw = 2.5
    nb_bins = 100
    c1 = (:grey35, 0.5)
    c_dry = (:orange3, 0.5)
    c_drymean = (:orange4)
    c_med = (:olivedrab, 0.3)
    c_medmean = (:darkolivegreen)
    c_wet = (:darkgreen, 0.5)
    c_wetmean = (:black) 
    c2 = :white#:lightgrey
    gc = :white

    ddry = data[data.s .< 0.36,:]
    dmed = data[(data.s .> 0.36) .& (data.s .< 0.61), :]
    dwet = data[data.s .> 0.61,:]

    fig = Figure(resolution = (800, 650))
    ax1 = Axis(fig[1, 1], xlabel = l[xquant1], ylabel = l[yquant], xlabelsize = lfs, ylabelsize = lfs, xgridcolor = gc, ygridcolor = gc)
    ax2 = Axis(fig[1, 2], xlabel = l[xquant2], xlabelsize = lfs, xgridcolor = gc, ygridcolor = gc)
    ax3 = Axis(fig[2, 1], xlabel = l[xquant3], ylabel = l[yquant], xlabelsize = lfs, ylabelsize = lfs, xgridcolor = gc, ygridcolor = gc) 
    ax4 = Axis(fig[2, 2], xlabel = l[xquant4], xlabelsize = lfs, ylabelsize = lfs, xgridcolor = gc, ygridcolor = gc)
        
    scatter!(ax1, dwet[!,xquant1], dwet[!,yquant], markersize = ms, color = c_wet)
    scatter!(ax1, dmed[!, xquant1], dmed[!,yquant], markersize = ms, color = c_med)
    scatter!(ax1, ddry[!,xquant1], ddry[!,yquant], markersize = ms, color = c_dry)
    lines!(ax1, sort(data, xquant1)[!,xquant1], mean_of_bins!(data, xquant1, yquant, nb_bins), color = c2, linewidth = lw)
    lines!(ax1, sort(ddry, xquant1)[!,xquant1], mean_of_bins!(ddry, xquant1, yquant, 60), color = c_drymean, linestyle = :dash, linewidth = lw)
    lines!(ax1, sort(dmed[1:70785,:], xquant1)[!,xquant1], mean_of_bins!(dmed[1:70785,:], xquant1, yquant, 65), color = c_medmean, linestyle = :dash, linewidth = lw)
    lines!(ax1, sort(dwet[1:15640,:], xquant1)[!,xquant1], mean_of_bins!(dwet[1:15640,:], xquant1, yquant, 68), color = c_wetmean, linestyle = :dash, linewidth = lw)

    scatter!(ax2, dwet[!,xquant2], dwet[!,yquant], markersize = ms, color = c_wet)
    scatter!(ax2, dmed[!, xquant2], dmed[!,yquant], markersize = ms, color = c_med)
    scatter!(ax2, ddry[!,xquant2], ddry[!,yquant], markersize = ms, color = c_dry)
    lines!(ax2, sort(data, xquant2)[!,xquant2], mean_of_bins!(data, xquant2, yquant, nb_bins), color = c2, linewidth = lw)
    lines!(ax2, sort(ddry, xquant2)[!,xquant2], mean_of_bins!(ddry, xquant2, yquant, 60), color = c_drymean, linestyle = :dash, linewidth = lw)
    lines!(ax2, sort(dmed[1:70785,:], xquant2)[!,xquant2], mean_of_bins!(dmed[1:70785,:], xquant2, yquant, 65), color = c_medmean, linestyle = :dash, linewidth = lw)
    lines!(ax2, sort(dwet[1:15640,:], xquant2)[!,xquant2], mean_of_bins!(dwet[1:15640,:], xquant2, yquant, 68), color = c_wetmean, linestyle = :dash, linewidth = lw)


    scatter!(ax3, dmed[!, xquant3], dmed[!,yquant], markersize = ms, color = c_med)
    scatter!(ax3, ddry[!,xquant3], ddry[!,yquant], markersize = ms, color = c_dry)
    scatter!(ax3, dwet[!,xquant3], dwet[!,yquant], markersize = ms, color = c_wet)
    lines!(ax3, sort(data, xquant3)[!,xquant3], mean_of_bins!(data, xquant3, yquant, nb_bins), color = c2, linewidth = lw)
    lines!(ax3, sort(ddry, xquant3)[!,xquant3], mean_of_bins!(ddry, xquant3, yquant, 60), color = c_drymean, linestyle = :dash, linewidth = lw)
    lines!(ax3, sort(dmed[1:70785,:], xquant3)[!,xquant3], mean_of_bins!(dmed[1:70785,:], xquant3, yquant, 65), color = c_medmean, linestyle = :dash, linewidth = lw)
    lines!(ax3, sort(dwet[1:15640,:], xquant3)[!,xquant3], mean_of_bins!(dwet[1:15640,:], xquant3, yquant, 68), color = c_wetmean, linestyle = :dash, linewidth = lw)

    scatter!(ax4, dwet[!,xquant4], dwet[!,yquant], markersize = ms, color = c_wet)
    scatter!(ax4, dmed[!, xquant4], dmed[!,yquant], markersize = ms, color = c_med)
    scatter!(ax4, ddry[!,xquant4], ddry[!,yquant], markersize = ms, color = c_dry)
    lines!(ax4, sort(data, xquant4)[!,xquant4], mean_of_bins!(data, xquant4, yquant, nb_bins), color = c2, linewidth = lw)
    lines!(ax4, sort(ddry, xquant4)[!,xquant4], mean_of_bins!(ddry, xquant4, yquant, 60), color = c_drymean, linestyle = :dash, linewidth = lw)
    lines!(ax4, sort(dmed[1:70785,:], xquant4)[!,xquant4], mean_of_bins!(dmed[1:70785,:], xquant4, yquant, 65), color = c_medmean, linestyle = :dash, linewidth = lw)
    lines!(ax4, sort(dwet[1:15640,:], xquant4)[!,xquant4], mean_of_bins!(dwet[1:15640,:], xquant4, yquant, 68), color = c_wetmean, linestyle = :dash, linewidth = lw)


    hidespines!(ax1, :t, :r)
    hidespines!(ax2, :t, :r)
    hidespines!(ax3, :t, :r)
    hidespines!(ax4, :t, :r)
    return fig
end
    

function fig_five_colormap(data::DataFrame = dcm, yquant::String = "PR", xquant1::String = "τ", xquant2::String = "α", xquant3::String = "spwp", xquant4::String = "r")
    fig = dcm |>  @vlplot(:point, x=:spwp, y=:PR, color=:s, width=400, height=400)
    return fig
end
    

function fig_five_lin(data::DataFrame = dcm, yquant::String = "χ", xquant1::String = "τ", xquant2::String = "α", xquant3::String = "e", xquant4::String = "r", xquant5::String = "p", xquant6::String = "eo")
    l = labels_lin_dict()
    lfs = 20
    ms = 5
    nb_bins = 100
    c1 = (:grey35, 0.5)
    c2 = :lightgrey
    gc = :white

    fig = Figure(resolution = (800, 900))
    ax1 = Axis(fig[1, 1], xlabel = l[xquant1], ylabel = l[yquant], xlabelsize = lfs, ylabelsize = lfs, xgridcolor = gc, ygridcolor = gc)
    ax2 = Axis(fig[1, 2], xlabel = l[xquant2], xlabelsize = lfs, xgridcolor = gc, ygridcolor = gc)
    ax3 = Axis(fig[2, 1], xlabel = l[xquant3], ylabel = l[yquant], xlabelsize = lfs, ylabelsize = lfs, xgridcolor = gc, ygridcolor = gc) 
    ax4 = Axis(fig[2, 2], xlabel = l[xquant4], xlabelsize = lfs, ylabelsize = lfs, xgridcolor = gc, ygridcolor = gc)
    ax5 = Axis(fig[3, 1], xlabel = l[xquant5], ylabel = l[yquant], xlabelsize = lfs, ylabelsize = lfs, xgridcolor = gc, ygridcolor = gc) 
    ax6 = Axis(fig[3, 2], xlabel = l[xquant6], xlabelsize = lfs, ylabelsize = lfs, xgridcolor = gc, ygridcolor = gc)
        
    scatter!(ax1, data[!,xquant1], data[!,yquant], markersize = ms, color = c1)
    lines!(ax1, sort(data, xquant1)[!,xquant1], mean_of_bins!(data, xquant1, yquant, nb_bins), color = c2)
    scatter!(ax2, data[!,xquant2], data[!,yquant], markersize = ms, color = c1)
    lines!(ax2, sort(data, xquant2)[!,xquant2], mean_of_bins!(data, xquant2, yquant, nb_bins), color = c2)
    scatter!(ax3, data[!,xquant3], data[!,yquant], markersize = ms, color = c1)
    lines!(ax3, sort(data, xquant3)[!,xquant3], mean_of_bins!(data, xquant3, yquant, nb_bins), color = c2)
    scatter!(ax4, data[!,xquant4], data[!,yquant], markersize = ms, color = c1)
    lines!(ax4, sort(data, xquant4)[!,xquant4], mean_of_bins!(data, xquant4, yquant, nb_bins), color = c2)
    scatter!(ax5, data[!,xquant5], data[!,yquant], markersize = ms, color = c1)
    lines!(ax5, sort(data, xquant5)[!,xquant5], mean_of_bins!(data, xquant5, yquant, nb_bins), color = c2)
    scatter!(ax6, data[!,xquant6], data[!,yquant], markersize = ms, color = c1)
    lines!(ax6, sort(data, xquant6)[!,xquant6], mean_of_bins!(data, xquant6, yquant, nb_bins), color = c2)
    hidespines!(ax1, :t, :r)
    hidespines!(ax2, :t, :r)
    hidespines!(ax3, :t, :r)
    hidespines!(ax4, :t, :r)
    hidespines!(ax5, :t, :r)
    hidespines!(ax6, :t, :r)
    return fig
end



function fig_two_lin(data::DataFrame = dcm)
    l = full_labels_dict()
    lfs = 20
    lw = 2.5
    c1 = :grey30
    c2 = :grey40

    f1 = Figure(resolution = (900, 400))
    ax1 = Axis(f1[1, 1], xlabel = l["s"], ylabel = L"\mathrm{Probability}\, \, \mathrm{density}\, \, \mathrm{function}", ylabelsize = lfs, xlabelsize = lfs, xgridvisible = false, ygridvisible = false)
    ax3 = Axis(f1[1, 2], xlabel = L"Mean water vapor pass $w$", xlabelsize = lfs, xgridvisible = false, ygridvisible = false)
    ax4 = Axis(f1[1,3], xlabel = L"Precipitation ratio $\chi$", xlabelsize = lfs, xgridvisible = false, ygridvisible = false)

    lines!(ax1, kde(data[!,"s"]), color = c1, linewidth = lw)    
    lines!(ax3, kde(data[!,"wo"]), color = c1, linewidth = lw, label = L"\mathrm{ocean}")
    lines!(ax3, kde(data[!,"wl"]), color = c1, linestyle = :dot, linewidth = lw, label = L"\mathrm{land}")
    lines!(ax4, kde(data[!,"χ"]), color = c1, linewidth = lw)

    axislegend(ax3, position = :rt, framevisible = false, labelsize = 20)
    hidespines!(ax1, :t, :r)
    hidespines!(ax3, :t, :r)
    hidespines!(ax4, :t, :r)
    #xlims!(ax3, (0,60))
    #xlims!(ax1, (-0.02, 0.82))
    colsize!(f1.layout, 1, Relative(1/3))
    colsize!(f1.layout, 2, Relative(1/3))
    colsize!(f1.layout, 3, Relative(1/3))
    return f1
end
    