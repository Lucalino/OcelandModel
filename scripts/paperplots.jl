using DrWatson
@quickactivate "Oceland Model"
using CairoMakie
using DataFrames
using CSV
using LaTeXStrings
using Colors
using Statistics
using KernelDensity
include(srcdir("utils.jl"))
include(srcdir("figure_labels.jl"))


dcm     = CSV.read(datadir("sims", "closed model/final/cm_smooth_tau_fixedpoints_runs100000_updated_ranges_final_all_quantities" * ".csv"), DataFrame)
dcm_mi   = CSV.read(datadir("sims", "mutual information/final/cm_rel_mi_100000_runs_final" * ".csv"), DataFrame)
dom     = CSV.read(datadir("sims", "open model/final/om_v2_sym_fixedpoints_runs100000_updated_ranges_final_all_quantities" * ".csv"), DataFrame)


function fig_two(data::DataFrame = dcm)
    l = full_labels_dict()
    lfs = 20
    lw = 2.5
    c1 = :grey30
    c2 = :grey40
    data.sresc = (data.s .- data.spwp) ./ (data.sfc .- data.spwp)

    f1 = Figure(resolution = (1100, 400))
    ax1 = Axis(f1[1, 1], xlabel = l["s"], ylabel = L"\mathrm{Probability}\, \, \mathrm{density}\, \, \mathrm{function}", ylabelsize = lfs, xlabelsize = lfs, xgridvisible = false, ygridvisible = false)
    ax2 = Axis(f1[1,2], xlabel = L"\tilde{s} = \frac{s - s_\mathrm{pwp}}{s_\mathrm{sfc}-s_\mathrm{pwp}}", xlabelsize = lfs, xgridvisible = false, ygridvisible = false)
    ax3 = Axis(f1[1, 3], xlabel = L"Mean water vapor pass $w$", xlabelsize = lfs, xgridvisible = false, ygridvisible = false)
    ax4 = Axis(f1[1,4], xlabel = L"Precipitation ratio $\chi$", xlabelsize = lfs, xgridvisible = false, ygridvisible = false)

    lines!(ax1, kde(data[!,"s"]), color = c1, linewidth = lw)    
    lines!(ax2, kde(data[!,"sresc"]), color = c1, linewidth = lw)
    vlines!(ax2, 0.0, linewidth = 2.0, linestyle = :dash, color = c2)
    vlines!(ax2, 1.0, linewidth = 2.0, linestyle = :dash, color = c2)    
    lines!(ax3, kde(data[!,"wo"]), color = c1, linewidth = lw, label = L"\mathrm{ocean}")
    lines!(ax3, kde(data[!,"wl"]), color = c1, linestyle = :dot, linewidth = lw, label = L"\mathrm{land}")
    lines!(ax4, kde(data[!,"PR"]), color = c1, linewidth = lw)

    axislegend(ax3, position = :lt, framevisible = false, labelsize = 20)
    hidespines!(ax1, :t, :r)
    hidespines!(ax2, :t, :r)
    hidespines!(ax3, :t, :r)
    hidespines!(ax4, :t, :r)
    xlims!(ax3, (0,60))
    xlims!(ax1, (-0.02, 0.82))
    colsize!(f1.layout, 1, Relative(5/24))
    colsize!(f1.layout, 2, Relative(5/24))
    colsize!(f1.layout, 3, Relative(1/3))
    colsize!(f1.layout, 4, Relative(1/4))
    return f1
end

function fig_three(data::DataFrame = dcm, statevar::String = "s", window_length::Int = 20000)
    fl = full_labels_dict()
    lfs = 20
    legendfs = 16

    colors = [colorant"#ff7f00", colorant"#33a02c", colorant"#6a3d9a",
              colorant"#73acd0", colorant" #215cc3", colorant"#e31a1c"]
    fluxes = ["eo", "El", "R", "Pl", "Po", "B"]
    fluxnames = [L"E_\mathrm{o}", L"E_\ell", L"R=A_\ell", L"P_\ell", L"P_\mathrm{o}", L"-A_\mathrm{o}"]
    
    fig = Figure(resolution = (600, 400))
    ax = Axis(fig[1, 1], xlabel = fl[statevar], ylabel = L"Mean water fluxes $F_i$ [mm/day]", 
                ylabelsize = lfs, xlabelsize = lfs, xgridcolor = :white, ygridcolor = :white)
    for i in 1:length(fluxes)
        lines!(ax, sort!(data, statevar)[!,statevar], movingaverage(data[:,fluxes[i]], window_length), label = fluxnames[i], color = colors[i])
    end
    #vlines!(ax, 0.61)
    #vlines!(ax, 0.36)
    axislegend(ax, position = :lc, framevisible = false, labelsize = legendfs)
    hidespines!(ax, :t, :r)
    ylims!(ax, (0,3.5))
    return fig
end

function fig_four(rel_mi_data::DataFrame = dcm_mi)
    sort!(rel_mi_data, "MI_rel", rev = true)
    l = short_labels_dict()
    n = nrow(rel_mi_data)
    xrange = collect(1:1:n)  

    f = Figure(resolution = (600, 450))
    ax = Axis(f[1,1], yscale = log10, xlabel = L"Model parameters $p_i$", ylabel = L"Mutual information index $I_{MI}\, (p_i\, ,\chi)$", xlabelsize = 20, ylabelsize = 20, xgridcolor = :white, ygridcolor = :white)
    ylims!(ax, 0.1, 500.0) # separate
    hlines!(ax, 1.0, linestyle = :dash, color = :grey35, label = "3σ significance threshold")
    scatter!(ax, xrange, rel_mi_data[:, "MI_rel"], marker = '*', markersize = 25, color = :grey20)
    ax.xticks = (1:1:n, rel_mi_data[:,"pnames"])
    axislegend(ax, framevisible = false, labelsize = 16)
    hidespines!(ax, :t, :r)
    return f

end

function fig_five(data::DataFrame = dcm, yquant::String = "PR", xquant1::String = "τ", xquant2::String = "α", xquant3::String = "spwp", xquant4::String = "r")
    l = labels_dict()
    lfs = 20
    ms = 5
    nb_bins = 100
    c1 = (:grey35, 0.5)
    c2 = :lightgrey
    gc = :white

    fig = Figure(resolution = (800, 650))
    ax1 = Axis(fig[1, 1], xlabel = l[xquant1], ylabel = l[yquant], xlabelsize = lfs, ylabelsize = lfs, xgridcolor = gc, ygridcolor = gc)
    ax2 = Axis(fig[1, 2], xlabel = l[xquant2], xlabelsize = lfs, xgridcolor = gc, ygridcolor = gc)
    ax3 = Axis(fig[2, 1], xlabel = l[xquant3], ylabel = l[yquant], xlabelsize = lfs, ylabelsize = lfs, xgridcolor = gc, ygridcolor = gc) 
    ax4 = Axis(fig[2, 2], xlabel = l[xquant4], xlabelsize = lfs, ylabelsize = lfs, xgridcolor = gc, ygridcolor = gc)
        
    scatter!(ax1, data[!,xquant1], data[!,yquant], markersize = ms, color = c1)
    lines!(ax1, sort(data, xquant1)[!,xquant1], mean_of_bins!(data, xquant1, yquant, nb_bins), color = c2)
    scatter!(ax2, data[!,xquant2], data[!,yquant], markersize = ms, color = c1)
    lines!(ax2, sort(data, xquant2)[!,xquant2], mean_of_bins!(data, xquant2, yquant, nb_bins), color = c2)
    scatter!(ax3, data[!,xquant3], data[!,yquant], markersize = ms, color = c1)
    lines!(ax3, sort(data, xquant3)[!,xquant3], mean_of_bins!(data, xquant3, yquant, nb_bins), color = c2)
    scatter!(ax4, data[!,xquant4], data[!,yquant], markersize = ms, color = c1)
    lines!(ax4, sort(data, xquant4)[!,xquant4], mean_of_bins!(data, xquant4, yquant, nb_bins), color = c2)
    hidespines!(ax1, :t, :r)
    hidespines!(ax2, :t, :r)
    hidespines!(ax3, :t, :r)
    hidespines!(ax4, :t, :r)
    return fig
end

function fig_six(data::DataFrame = dcm, param::String = "α", quant1::String = "Pl", quant2::String = "Po")
    l = labels_dict()
    lfs = 20
    ms = 5
    nb_bins = 100
    c1 = (:grey35, 0.5)
    c2 = :lightgrey

    fig = Figure(resolution = (800, 400))
    ax1 = Axis(fig[1, 1], xlabel = l[param], ylabel = l[quant1], xlabelsize = lfs, ylabelsize = lfs, xgridcolor = :white, ygridcolor = :white)
    ax2 = Axis(fig[1, 2], xlabel = l[param], ylabel = l[quant2], xlabelsize = lfs, ylabelsize = lfs, xgridcolor = :white, ygridcolor = :white)
    scatter!(ax1, data[!,param], data[!,quant1], markersize = ms, color = c1)
    lines!(ax1, sort(data, param)[!,param], mean_of_bins!(data, param, quant1, nb_bins), color = c2)
    scatter!(ax2, data[!,param], data[!,quant2], markersize = ms, color = c1)
    lines!(ax2, sort(data, param)[!,param], mean_of_bins!(data, param, quant2, nb_bins), color = c2)
    hidespines!(ax1, :t, :r)
    hidespines!(ax2, :t, :r)
    return fig

end

function fig_seven(data::DataFrame = dcm)
    l = full_labels_dict()
    lfs = 20
    legendfs = 16
    ms = 5
    ds = 12
    lw = 3.0
    n = 20000
    c1 = (:grey35, 0.5)
    c2 = :lightgrey
    gc = :white
    c = [:grey20, :grey60]
    s = collect(0.0:0.01:1.0)
    ep = 5.0
    spwp = [0.3, 0.4] 
    labels = [L"$s_\mathrm{pwp} = 0.3$", L"$s_\mathrm{pwp} = 0.4$"]
    s_bdot = mean(data[(data.spwp .> 0.29) .& (data.spwp .< 0.31), "s"])
    E_bdot = ep/2 * tanh( 10.0 * (s_bdot - (2*spwp[1]+0.3)/2 ) ) + ep/2
    s_odot = mean(data[(data.spwp .> 0.39) .& (data.spwp .< 0.41), "s"])
    E_odot = ep/2 * tanh( 10.0 * (s_odot - (2*spwp[2]+0.3)/2 ) ) + ep/2
    s_gdot = s_bdot
    E_gdot = ep/2 * tanh( 10.0 * (s_gdot - (2*spwp[2]+0.3)/2 ) ) + ep/2

    fig = Figure(resolution = (800, 400))
    ax1 = Axis(fig[1,1], xlabel = l["s"], ylabel = l["El"], ylabelsize = lfs, xlabelsize = lfs, xgridcolor = gc, ygridcolor = gc)   
    for n=1:length(spwp)
        E = [(ep/2 * tanh( 10.0 * (el - (2*spwp[n]+0.3)/2 ) ) + ep/2) for el in s]
        lines!(ax1, s, E, linewidth = lw, color = c[n], label = labels[n])
    end
    scatter!(ax1, Point2f(s_bdot, E_bdot), marker = :circle, markersize = ds, color = :blue)
    scatter!(ax1, Point2f(s_gdot, E_gdot), marker = :circle, markersize = ds, color = :green)
    scatter!(ax1, Point2f(s_odot, E_odot), marker = :circle, markersize = ds, color = :darkorange2)
    axislegend(ax1, position = :lt, framevisible = false, labelsize = legendfs)
    hidespines!(ax1, :t, :r)
    
    ax2 = Axis(fig[1,2], xlabel = l["spwp"], ylabel = l["s"], xlabelsize = lfs, ylabelsize = lfs, xgridcolor = :white, ygridcolor = :white)
    scatter!(ax2, data[!,"spwp"], data[!,"s"], markersize = ms, color = c1)
    lines!(ax2, sort(data, "spwp")[!,"spwp"], mean_of_bins!(data, "spwp", "s", 100), color = c2)
    #lines!(ax2, sort!(data, "spwp")[!,"spwp"], movingaverage(data[:,"s"],n), color = c2)
    scatter!(ax2, Point2f(spwp[1], s_bdot), marker = :circle, markersize = ds, color = :blue)
    scatter!(ax2, Point2f(spwp[2], s_odot), marker = :circle, markersize = ds, color = :darkorange2)
    hidespines!(ax2, :t, :r)
    return fig

end

function fig_eight(data::DataFrame = dom)
    l = full_labels_dict()
    r = ranges_dict()
    lfs = 20
    χlarger1 = data[data.PR .> 1,:]
    χsmaller1 = data[data.PR .< 1,:]
    label1 = L"\chi\, >\, 1\, \, (6.7\, %)"
    label2 = L"\chi\, <\, 1\, \, (93.3\,%)"
    lw1 = 2.0
    lw2 = 2.0
    col1 = :darkorange3
    col2 = :dodgerblue4
    yl = L"\mathrm{Probability}\, \, \mathrm{density}\, \, \mathrm{function}"

    f = Figure(resolution = (1200, 400))
    ax1 = Axis(f[1, 1], xlabel = l["w0"], xlabelsize = lfs, ylabel = yl, ylabelsize = lfs,  
                yticklabelcolor = col2, xgridvisible = false, ygridvisible = false)
    ax2 = Axis(f[1,1], yticklabelcolor = col1, yaxisposition = :right, 
                xgridvisible = false, ygridvisible = false)
    ax3 = Axis(f[1, 2], xlabel = l["α"], xlabelsize = lfs, yticklabelcolor = col2,
                xgridvisible = false, ygridvisible = false)    
    ax4 = Axis(f[1,2], yticklabelcolor = col1, yaxisposition = :right, 
                xgridvisible = false, ygridvisible = false)  
    ax5 = Axis(f[1, 3], xlabel = l["τ"], xlabelsize = lfs, yticklabelcolor = col2,
                xgridvisible = false, ygridvisible = false)
    ax6 = Axis(f[1,3], yticklabelcolor = col1, yaxisposition = :right, 
                xgridvisible = false, ygridvisible = false)

    l2 = lines!(ax1, kde(χsmaller1[!,"w0"]), color = col2, linewidth = lw2)
    l1 = lines!(ax2, kde(χlarger1[!,"w0"]), color = col1, linewidth = lw1)
    lines!(ax3, kde(χsmaller1[!,"α"]), color = col2, linewidth = lw2)
    lines!(ax4, kde(χlarger1[!,"α"]), color = col1, linewidth = lw1)
    lines!(ax5, kde(χsmaller1[!,"τ"]), color = col2, linewidth = lw2)
    lines!(ax6, kde(χlarger1[!,"τ"]), color = col1, linewidth = lw1)
 
    axislegend(ax5, [l2, l1], [label2, label1], framevisible = false, labelsize = 20)
    hidespines!(ax1, :t, :r)
    hidespines!(ax2, :t, :l, :b)
    hidespines!(ax3, :t, :r)
    hidespines!(ax4, :t, :l, :b)
    hidespines!(ax5, :t, :r)
    hidespines!(ax6, :t, :l, :b)
    hidexdecorations!(ax2)
    hidexdecorations!(ax4) 
    hidexdecorations!(ax6)
    xlims!(ax1, r["w0"])
    xlims!(ax2, r["w0"])
    xlims!(ax3, r["α"])
    xlims!(ax4, r["α"])
    xlims!(ax5, r["τ"])
    xlims!(ax6, r["τ"])

    return f
end







