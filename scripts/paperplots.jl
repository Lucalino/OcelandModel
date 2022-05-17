using DrWatson
@quickactivate "Oceland Model"

using CairoMakie
using DataFrames
using CSV
using LaTeXStrings
using Colors
using Statistics
using DynamicalSystems
using Random
using KernelDensity

include(srcdir("cm_analysis.jl"))
include(srcdir("om_analysis.jl"))
include(srcdir("utils.jl"))

dcm     = CSV.read(datadir("sims", "closed model pmscan/Final runs/cm_smooth_tau_eq_MC_fixedpoints_runs50000_updated_ranges_final-all_all_quantities" * ".csv"), DataFrame)
#dcmold   = CSV.read(datadir("sims", "closed model pmscan/cm_smooth_tau_eq_MC_fixedpoints_runs50000_updated_ranges_all_quantities" * ".csv"), DataFrame)
cm_mi   = CSV.read(datadir("sims", "mutual information/final/cm_rel_mi_50000_runs_final" * ".csv"), DataFrame)
#dcmoldold  = CSV.read(datadir("sims", "closed model pmscan/cm_smooth_tau_eq_MC_fixedpoints_runs10000_all_quantities" * ".csv"), DataFrame)
#cmpsens = CSV.read(datadir("sims", "closed model pmscan/cm_tau_10000runs_parameter_sensitivities_MI" * ".csv"), DataFrame)
dom = CSV.read(datadir("sims", "open model pmscan/om_v2_fixedpoints_runs50000_sym_updatedparams_all_quantities" * ".csv"), DataFrame)
#domasl = CSV.read(datadir("sims", "open model pmscan/asym/om_v2_MC_fixedpoints_runs10000_asym3L1=L3_all_quantities" * ".csv"), DataFrame)
#domasr = CSV.read(datadir("sims", "open model pmscan/asym/om_v2_MC_fixedpoints_runs10000_asymL1=3L3_all_quantities" * ".csv"), DataFrame)

#Used for PR(param)-plots and s(spwp) plot
function one_tile_plot(data::DataFrame, xquant::String, yquant::String)
    l = full_labels_dict()
    lfs = 20
    ms = 5
    nb_bins = 1000
    c1 = (:grey35, 0.8)
    c2 = :grey20
    fig = Figure(resolution = (500, 400))
    ax = Axis(fig[1, 1], xlabel = l[xquant], ylabel = l[yquant], ylabelsize = lfs, xlabelsize = lfs, xgridcolor = :white, ygridcolor = :white)
    scatter!(ax, data[!,xquant], data[!,yquant], markersize = ms, color = c1)
    #hlines!(ax, 1.0, linewidth = 2.0, linestyle = :dot, color = :grey80)
    #hlines!(ax, 0.51, linewidth = 2.0, linestyle = :dot, color = :grey80)
    #vlines!(ax, 0.3, linewidth = 2.0, linestyle = :dot, color = :grey80)
    #vlines!(ax, 0.4, linewidth = 2.0, linestyle = :dot, color = :grey80)
    #lines!(ax, sort(data, xquant)[!,xquant], cm_mean_of_bins!(data, xquant, yquant, nb_bins), color = c2)
    #lines!(ax, sort(data, xquant)[!,xquant], cm_rolling_average!(data, xquant, yquant, 50), color = :forestgreen, label = "bin size 50")
    #lines!(ax, sort(data, xquant)[!,xquant], cm_rolling_average!(data, xquant, yquant, nb_bins), color = :blue, label = "bin size 1000")
    lines!(ax, sort!(data, xquant)[!,xquant], movingaverage(data[!,yquant], 10000), color = c2)
    vlines!(ax,0.3)
    vlines!(ax,0.4)
    #hlines!(ax,0.88)
    #hlines!(ax,0.98)
    #lines!(ax, data[!,xquant], data[!,xquant], linestyle = :dot, color = :orange)
    hidespines!(ax, :t, :r)
    #axislegend(ax, position = :lb, framevisible = false, labelsize = 16)
    return fig
end

function two_tiles_plot(data::DataFrame, param::String, quant1::String, quant2::String)
    l = labels_dict()
    t = titles_dict()
    lfs = 20
    tfs = 16
    ms = 5
    nb_bins = 100
    c1 = (:grey35, 0.5)
    c2 = :grey25
    fig = Figure(resolution = (800, 400))
    ax1 = Axis(fig[1, 1], xlabel = l[param], ylabel = l[quant1], xlabelsize = lfs, ylabelsize = lfs, xgridcolor = :white, ygridcolor = :white)
    ax2 = Axis(fig[1, 2], xlabel = l[param], ylabel = l[quant2], xlabelsize = lfs, ylabelsize = lfs, xgridcolor = :white, ygridcolor = :white)
    scatter!(ax1, data[!,param], data[!,quant1], markersize = ms, color = c1)
    lines!(ax1, sort(data, param)[!,param], cm_mean_of_bins!(data, param, quant1, nb_bins), color = c2)
    scatter!(ax2, data[!,param], data[!,quant2], markersize = ms, color = c1)
    lines!(ax2, sort(data, param)[!,param], cm_mean_of_bins!(data, param, quant2, nb_bins), color = c2)
    hidespines!(ax1, :t, :r)
    hidespines!(ax2, :t, :r)
    return fig
end


#Used for StateVariables(param)-plots
function three_tiles_plot(data::DataFrame, param::String, quant1::String, quant2::String, quant3::String)
    l = labels_dict()
    t = titles_dict()
    lfs = 20
    tfs = 16
    ms = 5
    nb_bins = 100
    c1 = (:grey35, 0.8)
    c2 = :grey20
    fig = Figure(resolution = (600, 800))
    ax1 = Axis(fig[1, 1], ylabel = l[quant1], ylabelsize = lfs, title = t[quant1], titlesize = tfs, xgridcolor = :white, ygridcolor = :white)
    ax2 = Axis(fig[2, 1], ylabel = l[quant2], ylabelsize = lfs, title = t[quant2], titlesize = tfs, xgridcolor = :white, ygridcolor = :white)
    ax3 = Axis(fig[3, 1], xlabel = l[param], ylabel = l[quant3], xlabelsize = lfs, ylabelsize = lfs, title = t[quant3], titlesize = tfs, xgridcolor = :white, ygridcolor = :white)
    scatter!(ax1, data[!,param], data[!,quant1], markersize = ms, color = c1)
    lines!(ax1, sort(data, param)[!,param], cm_mean_of_bins!(data, param, quant1, nb_bins), color = c2)
    scatter!(ax2, data[!,param], data[!,quant2], markersize = ms, color = c1)
    lines!(ax2, sort(data, param)[!,param], cm_mean_of_bins!(data, param, quant2, nb_bins), color = c2)
    scatter!(ax3, data[!,param], data[!,quant3], markersize = ms, color = c1)
    lines!(ax3, sort(data, param)[!,param], cm_mean_of_bins!(data, param, quant3, nb_bins), color = c2)
    hidespines!(ax1, :t, :r)
    hidespines!(ax2, :t, :r)
    hidespines!(ax3, :t, :r)
    return fig
end

function four_tiles_plot(data::DataFrame, yquant::String, xquant1::String, xquant2::String, xquant3::String, xquant4::String)
    l = labels_dict()
    lfs = 20
    ms = 5
    nb_bins = 100
    c1 = (:grey35, 0.5)
    c2 = :grey25
    fig = Figure(resolution = (800, 650))
    ax1 = Axis(fig[1, 1], xlabel = l[xquant1], ylabel = l[yquant], xlabelsize = lfs, ylabelsize = lfs, xgridcolor = :white, ygridcolor = :white)
    ax2 = Axis(fig[1, 2], xlabel = l[xquant2], xlabelsize = lfs, xgridcolor = :white, ygridcolor = :white)
    ax3 = Axis(fig[2, 1], xlabel = l[xquant3], ylabel = l[yquant], xlabelsize = lfs, ylabelsize = lfs, xgridcolor = :white, ygridcolor = :white)
    ax4 = Axis(fig[2, 2], xlabel = l[xquant4], xlabelsize = lfs, ylabelsize = lfs, xgridcolor = :white, ygridcolor = :white)
    scatter!(ax1, data[!,xquant1], data[!,yquant], markersize = ms, color = c1)
    lines!(ax1, sort(data, xquant1)[!,xquant1], cm_mean_of_bins!(data, xquant1, yquant, nb_bins), color = c2)
    #lines!(ax1, sort!(data, xquant1)[!,xquant1], movingaverage(data[:,yquant],10000), color = c2)
    scatter!(ax2, data[!,xquant2], data[!,yquant], markersize = ms, color = c1)
    lines!(ax2, sort(data, xquant2)[!,xquant2], cm_mean_of_bins!(data, xquant2, yquant, nb_bins), color = c2)
    #lines!(ax2, sort!(data, xquant2)[!,xquant2], movingaverage(data[:,yquant],10000), color = c2)
    scatter!(ax3, data[!,xquant3], data[!,yquant], markersize = ms, color = c1)
    lines!(ax3, sort(data, xquant3)[!,xquant3], cm_mean_of_bins!(data, xquant3, yquant, nb_bins), color = c2)
    #lines!(ax3, sort!(data, xquant3)[!,xquant3], movingaverage(data[:,yquant],10000), color = c2)
    scatter!(ax4, data[!,xquant4], data[!,yquant], markersize = ms, color = c1)
    lines!(ax4, sort(data, xquant4)[!,xquant4], cm_mean_of_bins!(data, xquant4, yquant, nb_bins), color = c2)
    #lines!(ax4, sort!(data, xquant4)[!,xquant4], movingaverage(data[:,yquant],10000), color = c2)
    hidespines!(ax1, :t, :r)
    hidespines!(ax2, :t, :r)
    hidespines!(ax3, :t, :r)
    hidespines!(ax4, :t, :r)
    return fig
end


#Used for Fluxes(param)-plots
function six_tiles_plot(data::DataFrame, param::String, quant1::String, quant2::String, quant3::String, quant4::String, quant5::String, quant6::String, nb_bins::Integer = 100) 
    l = labels_dict()
    t = titles_dict()
    lfs = 20
    tfs = 16
    ms = 5
    c1 = (:grey35, 0.8)
    c2 = :grey20
    fig = Figure(resolution = (800, 1000))
    ax1 = Axis(fig[1, 1], ylabel = l[quant1], xlabelsize = lfs, ylabelsize = lfs, title = t[quant1], titlesize = tfs, xgridcolor = :white, ygridcolor = :white)
    ax2 = Axis(fig[1, 2], ylabel = l[quant2], xlabelsize = lfs, ylabelsize = lfs, title = t[quant2], titlesize = tfs, xgridcolor = :white, ygridcolor = :white)
    ax3 = Axis(fig[2, 1], ylabel = l[quant3], xlabelsize = lfs, ylabelsize = lfs, title = t[quant3], titlesize = tfs, xgridcolor = :white, ygridcolor = :white)
    ax4 = Axis(fig[2, 2], ylabel = l[quant4], xlabelsize = lfs, ylabelsize = lfs, title = t[quant4], titlesize = tfs, xgridcolor = :white, ygridcolor = :white)
    ax5 = Axis(fig[3, 1], xlabel = l[param], ylabel = l[quant5], xlabelsize = lfs, ylabelsize = lfs, title = t[quant5], titlesize = tfs, xgridcolor = :white, ygridcolor = :white)
    ax6 = Axis(fig[3, 2], xlabel = l[param], ylabel = l[quant6], xlabelsize = lfs, ylabelsize = lfs, title = t[quant6], titlesize = tfs, xgridcolor = :white, ygridcolor = :white)
    scatter!(ax1, data[!,param], data[!,quant1], markersize = ms, color = c1)
    lines!(ax1, sort(data, param)[!,param], cm_mean_of_bins!(data, param, quant1, nb_bins), color = c2)
    scatter!(ax2, data[!,param], data[!,quant2], markersize = ms, color = c1)
    lines!(ax2, sort(data, param)[!,param], cm_mean_of_bins!(data, param, quant2, nb_bins), color = c2)
    scatter!(ax3, data[!,param], data[!,quant3], markersize = ms, color = c1)
    lines!(ax3, sort(data, param)[!,param], cm_mean_of_bins!(data, param, quant3, nb_bins), color = c2)
    scatter!(ax4, data[!,param], data[!,quant4], markersize = ms, color = c1)
    lines!(ax4, sort(data, param)[!,param], cm_mean_of_bins!(data, param, quant4, nb_bins), color = c2)
    scatter!(ax5, data[!,param], data[!,quant5], markersize = ms, color = c1)
    lines!(ax5, sort(data, param)[!,param], cm_mean_of_bins!(data, param, quant5, nb_bins), color = c2)
    scatter!(ax6, data[!,param], data[!,quant6], markersize = ms, color = c1)
    lines!(ax6, sort(data, param)[!,param], cm_mean_of_bins!(data, param, quant6, nb_bins), color = c2)
    hidespines!(ax1, :t, :r)
    hidespines!(ax2, :t, :r)
    hidespines!(ax3, :t, :r)
    hidespines!(ax4, :t, :r)
    hidespines!(ax5, :t, :r)
    hidespines!(ax6, :t, :r)
    return fig
end


function sensitivity_plot(data::DataFrame, sens_version::Int, model_version::String, nb_bins::Int = 100)
    l = labels_norm_dict()
    lfs = 20

    if model_version == "cm"
        p = ["ϵ", "α", "eo", "spwp", "ep",  "τ", "a", "b", "wsat", "nZr"]
        n = length(p)
        colors = [colorant"#a6cee3", colorant"#1f78b4", colorant"#b2df8a", 
                  colorant"#33a02c", colorant"#fb9a99", colorant"#e31a1c",
                  colorant"#fdbf6f", colorant"#ff7f00", colorant"#cab2d6",
                  colorant"#6a3d9a"]

    elseif model_version == "om"
        p = ["ϵ", "α", "eo", "spwp", "ep", "τ", "a", "b", "wsat", "nZr", "w0"]
        n = length(p)
        colors = [colorant"#a6cee3", colorant"#1f78b4", colorant"#b2df8a", 
                  colorant"#33a02c", colorant"#fb9a99", colorant"#e31a1c",
                  colorant"#fdbf6f", colorant"#ff7f00", colorant"#cab2d6",
                  colorant"#6a3d9a", colorant"#ffff99"]
    end

    fig = Figure(resolution = (800, 600))
    ax = Axis(fig[1, 1], xlabel = L"Normalized parameter values $\hat{p}$", ylabel = L"Precipitation ratio $PR$", ylabelsize = lfs, xlabelsize = lfs)
    
    if sens_version == 1
        measure = "σ_"
        ax.title = "σ method (lower values -> higher sensitivity)"
        for i = 1:n
            sensitivity_v1!(data, p[i])
            data  = sort!(data, p[i])
            pnorm = (data[:,p[i]] .- minimum(data[:,p[i]])) ./ (maximum(data[:,p[i]]) .- minimum(data[:,p[i]]))
            lines!(ax, pnorm, data[:,string(measure,p[i])], color = (colors[i], 0.6), linewidth = 2)
        end

        for i = 1:n
            hlines!(ax, DataFrames.mean(data[:,string(measure, p[i])]), color = colors[i], linewidth = 3, label = l[p[i]])
        end
    
    elseif sens_version == 2
        measure = "Δμ_"
        ax.title = "Δμ method (higher values -> higher sensitivity)"
        for i = 1:n
            sensitivity_v2!(data, p[i])
            data = sort!(data, p[i])
            pnorm = (data[:,p[i]] .- minimum(data[:,p[i]])) ./ (maximum(data[:,p[i]]) .- minimum(data[:,p[i]]))
            lines!(ax, pnorm, data[:, string(measure,p[i])], color = (colors[i], 0.6), linewidth = 2)
        end

        for i = 1:n
            hlines!(ax, DataFrames.mean(data[:,string(measure, p[i])]), color = colors[i], linewidth = 3, label = l[p[i]])
        end
    
    else
        println("Please specify which sensitivity version should be used (1 or 2).")
    end

    fig[1, 2] = Legend(fig, ax, "Parameters", framevisible = false)
    hidespines!(ax, :t, :r)
    return fig
end

function rel_mi_plot(rel_mi_data::DataFrame)
    sort!(rel_mi_data, "MI_rel", rev = true)
    l = short_labels_dict()
    n = nrow(rel_mi_data)

    l_arr = Vector{String}()
    for i = 1:n
        push!(l_arr, l[rel_mi_data[i,"pnames"]])
    end
    
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

function fluxes_plot(data::DataFrame, statevar::String, bin_length::Int = 10000, nb_bins::Int = 200)
    fl = full_labels_dict()
    lfs = 20
    legendfs = 16
    s = collect(minimum(data.s):0.01:maximum(data.s)+0.01)
    wo = collect(minimum(data.wo):1.0:maximum(data.wo))
    wl = collect(minimum(data.wl):1.0:maximum(data.wl))
    colors = [colorant"#ff7f00", colorant"#33a02c", colorant"#6a3d9a",
              colorant"#73acd0", colorant" #215cc3", colorant"#e31a1c"] #  #1f78b4 #a6cee3

    fig = Figure(resolution = (600, 400))
    ax = Axis(fig[1, 1], xlabel = fl[statevar], ylabel = L"Mean water fluxes $F_i$ [mm/day]", ylabelsize = lfs, xlabelsize = lfs,
                xgridcolor = :white, ygridcolor = :white) #title = "bin length: $(bin_length)",
    fluxes = ["eo", "El", "R", "Pl", "Po", "B"]
    fluxnames = [L"E_\mathrm{o}", L"E_\ell", L"R=A_\ell", L"P_\ell", L"P_\mathrm{o}", L"-A_\mathrm{o}"]
    ls = [:solid, :solid, :solid, :solid, :solid, :solid]

    # if statevar == "s"
    #     band!(ax, s, [(2.0/2 * tanh( 10 * (el - 1.38/2 ) ) + 2.0/2) for el in s], [(6.0/2 * tanh( 10 * (el - 0.7/2 ) ) + 6.0/2) for el in s], color = (:lightgreen, 0.4))
    #     scatter!(ax, data.s, data.El, color = (colorant"#33a02c", 0.4), markersize = 2.0)
    # elseif statevar == "wl"
    #     band!(ax, wl, [exp(15.6*(el / 65.0 - 0.522)) for el in wl], [exp(11.4*(el / 80.0 - 0.603)) for el in wl], color = (:lightblue1, 0.8))
    #     scatter!(ax, data.wl, data.Pl, color = (colorant"#1f78b4", 0.6), markersize = 2.0)
    # elseif statevar == "wo"
    #     band!(ax, wo, [exp(15.6*(el / 65.0 - 0.522)) for el in wo], [exp(11.4*(el / 80.0 - 0.603)) for el in wo], color = (:lightblue1, 0.8))
    #     scatter!(ax, data.wo, data.Po, color = (colorant"#1f78b4", 0.6), markersize = 2.0)
    # else
    #     println("Some problem has occured. Check the code.")
    # end

    
    for i in 1:length(fluxes)
        #lines!(ax, sort(data, statevar)[!,statevar], cm_mean_of_bins!(data, statevar, fluxes[i], nb_bins), label = fluxes[i], color = colors[i], linestyle = ls[i])
        #lines!(ax, sort(data, statevar)[!,statevar], cm_rolling_average!(data, statevar, fluxes[i], nb_bins), label = fluxes[i], color = colors[i], linestyle = ls[i])
        lines!(ax, sort!(data, statevar)[!,statevar], movingaverage(data[:,fluxes[i]], bin_length), label = fluxnames[i], color = colors[i], linestyle = ls[i])
        #scatter!(ax, dcm[:, statevar], dcm[:,fluxes[i]], color = (colors[i], 0.5), markersize = 2.0)
    end
    #lines!(ax, s, [(4.5/2 * tanh( 10 * (el - 0.7/2 ) ) + 4.5/2) for el in s], linestyle = :dot, color = :green)
    #lines!(ax, s, [(4.1/2 * tanh( 10 * (el - 1.38/2 ) ) + 4.1/2) for el in s], linestyle = :dot, color = :green)
    #lines!(ax, w, [exp(15.6*(x / 65.0 - 0.522)) for el in w], linestyle = :dot, color = :blue)
    #lines!(ax, w, [exp(11.4*(el / 80.0 - 0.603)) for el in w], linestyle = :dot, color = :blue)
    #vlines!(ax, 0.61)
    #vlines!(ax, 0.36)
    axislegend(ax, position = :lc, framevisible = false, labelsize = legendfs)
    hidespines!(ax, :t, :r)
    ylims!(ax, (0,3.5))
    return fig
end

function pdf_plots(data::DataFrame, type::String = "kde", norm::Symbol = :pdf)
    l = full_labels_dict()
    lfs = 20
    ms = 5
    lw = 2.5
    f1 = Figure(resolution = (1100, 400))
    ax1 = Axis(f1[1, 1], xlabel = l["s"], ylabel = L"\mathrm{Probability}\, \, \mathrm{density}\, \, \mathrm{function}", ylabelsize = lfs, xlabelsize = lfs, xgridvisible = false, ygridvisible = false)
    
    if type == "kde"
        # density!(data[!,"s"], color = (:darkorange, 0.8))
        # density!(data[!,"s"], color = (:grey30, 0.8), strokearound = true, strokewidth = 1) #(:grey20, 0.8))
        lines!(kde(data[!,"s"]), color = :grey30, linewidth = lw)
    elseif type == "hist"
        hist!(ax1, data[!,"s"], bins = 100, normalization = norm, color = (:darkorange, 0.8))
    end


    ax2 = Axis(f1[1,2], xlabel = L"\tilde{s} = \frac{s - s_\mathrm{pwp}}{s_\mathrm{sfc}-s_\mathrm{pwp}}", xlabelsize = lfs, xgridvisible = false, ygridvisible = false)
    data.sresc = (data.s .- data.spwp) ./ (data.sfc .- data.spwp)

    if type == "kde"
        # density!(data[!,"sresc"], color = (:white, 0.8), strokearound = true, strokewidth = 1) # color = (:grey20, 0.8))
        lines!(kde(data[!,"sresc"]), color = :grey30, linewidth = lw)
        vlines!(ax2, 0.0, linewidth = 2.0, linestyle = :dash, color = :grey40)
        vlines!(ax2, 1.0, linewidth = 2.0, linestyle = :dash, color = :grey40)
    elseif type == "hist"
        hist!(ax2, data[!,"sresc"], bins = 100, normalization = norm, color = (:darkorange, 0.8))
    end

    ax3 = Axis(f1[1, 3], xlabel = L"Mean water vapor pass $w$", xlabelsize = lfs, xgridvisible = false, ygridvisible = false)
    
    if type == "kde"
        # density!(data[!,"wo"], color = (colorant" #215cc3", 0.8), label = "ocean")
        # density!(data[!,"wl"], color = (colorant"#73acd0", 0.8), label = "land")
        #density!(data[!,"wo"], color = (:white, 0.8), strokearound = true, strokewidth = 1, label = L"\mathrm{ocean}") # color = (:grey20, 0.8))
        #density!(data[!,"wl"], color = (:white, 0.6), strokearound = true, strokewidth = 1, label = L"\mathrm{land}") # color = (:grey50, 0.8))
        lines!(kde(data[!,"wo"]), color = :grey30, linewidth = lw, label = L"\mathrm{ocean}")
        lines!(kde(data[!,"wl"]), color = :grey30, linestyle = :dot, linewidth = lw, label = L"\mathrm{land}")
    elseif type == "hist"
        hist!(ax3, data[!,"wo"], bins = 100, normalization = norm, label = L"\mathrm{ocean}", color = (colorant" #215cc3", 0.8))
        hist!(ax3, data[!,"wl"], bins = 100, normalization = norm, label = L"\mathrm{land}", color = (colorant"#73acd0", 0.8))   
    end

    ax4 = Axis(f1[1,4], xlabel = L"\chi", xlabelsize = lfs, xgridvisible = false, ygridvisible = false)

    if type == "kde"
        #density!(data[!,"PR"], color = (:grey30, 0.8), strokearound = true, strokewidth = 1)
        lines!(kde(data[!,"PR"]), color = :grey30, linewidth = lw)
    end

    axislegend(ax3, position = :lt, framevisible = false, labelsize = 20)
    hidespines!(ax1, :t, :r)
    hidespines!(ax2, :t, :r)
    hidespines!(ax3, :t, :r)
    hidespines!(ax4, :t, :r)
    xlims!(ax3, (8,60))

    colsize!(f1.layout, 1, Relative(5/24))
    colsize!(f1.layout, 2, Relative(5/24))
    colsize!(f1.layout, 3, Relative(4/12))
    colsize!(f1.layout, 4, Relative(3/12))

    return f1
end

function om_pdf_statevars(data::DataFrame, type::String = "kde", norm::Symbol = :pdf)
    l = full_labels_dict()
    lfs = 20
    ms = 5
    f1 = Figure(resolution = (1000, 400))
    ax1 = Axis(f1[1, 1], xlabel = l["s"], ylabel = L"\mathrm{Probability}\, \, \mathrm{density}\, \, \mathrm{function}", ylabelsize = lfs, xlabelsize = lfs, xgridvisible = false, ygridvisible = false)
    
    if type == "kde"
        # density!(data[!,"s"], color = (:darkorange, 0.8))
        density!(data[!,"s"], color = (:grey20, 0.8))
    elseif type == "hist"
        hist!(ax1, data[!,"s"], bins = 100, normalization = norm, color = (:darkorange, 0.8))
    end
    xlims!(ax1, (-0.1,1.1))

    ax2 = Axis(f1[1,2], xlabel = L"\tilde{s} = \frac{s - s_\mathrm{pwp}}{s_\mathrm{sfc}-s_\mathrm{pwp}}", xlabelsize = lfs, xgridvisible = false, ygridvisible = false)
    data.sresc = (data.s .- data.spwp) ./ (data.sfc .- data.spwp)

    if type == "kde"
        density!(data[!,"sresc"], color = (:grey20, 0.8))
        vlines!(ax2, 0.0, linewidth = 2.0, linestyle = :dot, color = :grey40)
        vlines!(ax2, 1.0, linewidth = 2.0, linestyle = :dot, color = :grey40)
    elseif type == "hist"
        hist!(ax2, data[!,"sresc"], bins = 100, normalization = norm, color = (:darkorange, 0.8))
    end
    xlims!(ax2, (-1.4,2.9))

    ax3 = Axis(f1[1, 3], xlabel = L"Mean water vapor pass $w$", xlabelsize = lfs, xgridvisible = false, ygridvisible = false)
    
    if type == "kde"
        # density!(data[!,"wo"], color = (colorant" #215cc3", 0.8), label = "ocean")
        # density!(data[!,"wl"], color = (colorant"#73acd0", 0.8), label = "land")
        oce2 = density!(data[!,"w3"], color = (:grey50, 0.5), strokecolor = :grey30, strokewidth = 3.5)
        land = density!(data[!,"w2"], color = (:grey50, 0.5), strokecolor = :grey30, strokewidth = 2)
        oce1 = density!(data[!,"w1"], color = (:grey50, 0.5), strokecolor = :grey25, strokewidth = 0.75)


    elseif type == "hist"
        hist!(ax3, data[!,"wo"], bins = 100, normalization = norm, label = L"\mathrm{ocean}", color = (colorant" #215cc3", 0.8))
        hist!(ax3, data[!,"wl"], bins = 100, normalization = norm, label = L"\mathrm{land}", color = (colorant"#73acd0", 0.8))   
    end

    axislegend(ax3, [oce2, land, oce1], [L"\mathrm{ocean}\, \,  2",  L"\mathrm{land}", L"\mathrm{ocean}\, \,  1"], position = :lt, framevisible = false, labelsize = 20)
    hidespines!(ax1, :t, :r)
    hidespines!(ax2, :t, :r)
    hidespines!(ax3, :t, :r)
    xlims!(ax3, (-10,90))

    return f1
end

function om_pdf_params(data1::DataFrame, label1::LaTeXString, data2::DataFrame, label2::LaTeXString, p::String)
    l = full_labels_dict()
    r = ranges_dict()
    lfs = 20
    ms = 5

    f = Figure(resolution = (600, 400))

    ax1 = Axis(f[1, 1], xlabel = l[p], xlabelsize = lfs,
                ylabel = L"\mathrm{Probability}\, \, \mathrm{density}\, \, \mathrm{function}", ylabelsize = lfs,  
                yticklabelcolor = :grey20,
                xgridvisible = false, ygridvisible = false)
    d1 = density!(data2[!,p], color = (:grey20, 0.8))

    ax2 = Axis(f[1,1], yticklabelcolor = :grey50, yaxisposition = :right, xgridvisible = false, ygridvisible = false)
    d2 = density!(data1[!,p], color = (:grey50, 0.8))
    hidexdecorations!(ax2)  

    Legend(f[1,2], [d1, d2], [label2, label1], framevisible = false, labelsize = 20)
    #axislegend(ax1, position = :lt, framevisible = false, labelsize = 20)
    hidespines!(ax1, :t, :r)
    hidespines!(ax2, :t, :l, :b)
    xlims!(ax1, r[p])
    xlims!(ax2, r[p])

    return f
end

function om_pdf_w0_alpha_tau(data1::DataFrame, label1::LaTeXString, data2::DataFrame, label2::LaTeXString)
    l = full_labels_dict()
    r = ranges_dict()
    lfs = 20
    ms = 5
    d1lw = 2.0
    d2lw = 2.0
    d1col = :darkorange3
    d2col = :dodgerblue4
    d1ls  = :solid

    f = Figure(resolution = (1200, 400))

    ax1 = Axis(f[1, 1], xlabel = l["w0"], xlabelsize = lfs,
                ylabel = L"\mathrm{Probability}\, \, \mathrm{density}\, \, \mathrm{function}", ylabelsize = lfs,  
                yticklabelcolor = d2col,
                xgridvisible = false, ygridvisible = false)
    d1 = lines!(kde(data2[!,"w0"]), color = d2col, linewidth = d2lw)

    ax2 = Axis(f[1,1], yticklabelcolor = d1col, yaxisposition = :right, xgridvisible = false, ygridvisible = false)
    d2 = lines!(kde(data1[!,"w0"]), color = d1col, linestyle = d1ls, linewidth = d1lw)
    hidexdecorations!(ax2)  

    ax3 = Axis(f[1, 2], xlabel = l["α"], xlabelsize = lfs,
                yticklabelcolor = d2col,
                xgridvisible = false, ygridvisible = false)
    d3 = lines!(kde(data2[!,"α"]), color = d2col, linewidth = d2lw)

    ax4 = Axis(f[1,2], yticklabelcolor = d1col, yaxisposition = :right, xgridvisible = false, ygridvisible = false)
    d4 = lines!(kde(data1[!,"α"]), color = d1col, linestyle = d1ls, linewidth = d1lw)
    hidexdecorations!(ax4) 

    ax5 = Axis(f[1, 3], xlabel = l["τ"], xlabelsize = lfs,
                yticklabelcolor = d2col,
                xgridvisible = false, ygridvisible = false)
    d5 = lines!(kde(data2[!,"τ"]), color = d2col, linewidth = d2lw)

    ax6 = Axis(f[1,3], yticklabelcolor = d1col, yaxisposition = :right, xgridvisible = false, ygridvisible = false)
    d6 = lines!(kde(data1[!,"τ"]), color = d1col, linestyle = d1ls, linewidth = d1lw)
    hidexdecorations!(ax6) 
    axislegend(ax5, [d1, d2], [label2, label1], framevisible = false, labelsize = 20)
    #Legend(f[1,3], [d1, d2], [label2, label1], framevisible = false, labelsize = 20)
    #axislegend(ax1, position = :lt, framevisible = false, labelsize = 20)
    hidespines!(ax1, :t, :r)
    hidespines!(ax2, :t, :l, :b)
    hidespines!(ax3, :t, :r)
    hidespines!(ax4, :t, :l, :b)
    hidespines!(ax5, :t, :r)
    hidespines!(ax6, :t, :l, :b)
    xlims!(ax1, r["w0"])
    xlims!(ax2, r["w0"])
    xlims!(ax3, r["α"])
    xlims!(ax4, r["α"])
    xlims!(ax5, r["τ"])
    xlims!(ax6, r["τ"])

    return f
end


function rescaled_s_pdf_plot(data::DataFrame)
    data.sresc = (data.s .- data.spwp) ./ (data.sfc .- data.spwp)
    lfs = 20
    ms = 5
    f = Figure(resolution = (500, 400))
    ax = Axis(f[1, 1], xlabel = L"Rescaled soil moisture $\frac{s - s_\mathrm{pwp}}{s_\mathrm{sfc}-s_\mathrm{pwp}}$", ylabel = L"\mathrm{Probability}\, \, \mathrm{density}\, \, \mathrm{function}", ylabelsize = lfs, xlabelsize = lfs, xgridvisible = false, ygridvisible = false)
    density!(data[!,"sresc"], color = (:grey20, 0.8))
    hidespines!(ax, :t, :r)
    vlines!(ax, 0.0, linewidth = 2.0, linestyle = :dot, color = :grey40)
    vlines!(ax, 1.0, linewidth = 2.0, linestyle = :dot, color = :grey40)
    return f
end


function El_plot()
    l = full_labels_dict()
    lfs = 20
    legendfs = 16
    lw = 3.0
    c = [:grey20, :grey60] #[:grey20, :grey40, :grey60, :grey70]
    s = collect(0.0:0.01:1.0)
    ep = 5.0
    spwp = [0.3, 0.4] #[0.2, 0.3, 0.4, 0.5]
    labls = [L"$s_\mathrm{pwp} = 0.3$", L"$s_\mathrm{pwp} = 0.4$"] #[L"$s_\mathrm{pwp} = 0.2$", L"$s_\mathrm{pwp} = 0.3$", L"$s_\mathrm{pwp} = 0.4$", L"$s_\mathrm{pwp} = 0.5$"]
    fig = Figure(resolution = (500, 400))
    ax = Axis(fig[1,1], xlabel = l["s"], ylabel = l["El"], ylabelsize = lfs, xlabelsize = lfs, xgridcolor = :white, ygridcolor = :white)

    for n=1:length(spwp)
        E = [(ep/2 * tanh( 10.0 * (el - (2*spwp[n]+0.3)/2 ) ) + ep/2) for el in s]
        lines!(ax, s, E, linewidth = lw, color = c[n], label = labls[n])
    end
    vlines!(ax, 0.448, linewidth = 2.0, linestyle = :dot, color = :grey80)
    vlines!(ax, 0.535, linewidth = 2.0, linestyle = :dot, color = :grey80)
    hidespines!(ax, :t, :r)
    axislegend(ax, position = :lt, framevisible = false, labelsize = legendfs)
    return fig
end


function R_plot()
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
    #vlines!(ax, 0.43, linewidth = 2.0, linestyle = :dot, color = :grey80)
    #vlines!(ax, 0.51, linewidth = 2.0, linestyle = :dot, color = :grey80)
    hidespines!(ax, :t, :r)
    axislegend(ax, position = :lt, framevisible = false, labelsize = legendfs)
    return fig
end

function ranges_dict()
    d = Dict(
        "α"  => (-0.2, 1.2),
        "w0" => (-10, 90),
        "τ"  => (-0.2, 4.5),
        "PR" => (-0.2, 2.0),
        "u_ms"=>(-1, 11),
        "L_km"=>(100, 2100)
    )
    return d
end

function short_labels_dict()
    d = Dict{String, LaTeXString}(
        "spwp" => "s_pwp",
        "sfc"  => "s_fc",
        "ep"   => "e_p", 
        "eo"   => "e_o",
        "ϵ"    => "ϵ",
        "r"    => "r",  
        "α"    => "α", 
        "nZr"  => "nZ_r",   
        "a"    => "a",  
        "b"    => "b",
        "wsat" => "w_sat", 
        "u"    => "u",
        "L"    => "L",
        "Pl"   => "P_l",
        "Po"   => L"$P_\mathrm{o}$",
        "El"   => L"$E_\mathrm{l}$",
        "R"    => L"$R$",
        "PR"   => L"$\chi$",
        "Ptot" => L"$P_\mathrm{mean}$",
        "A"    => L"$A_\mathrm{l}$",
        "B"    => L"$-A_\mathrm{o}$", 
        "τ"    => "τ",
        "τ_time"=> L"\tau^{-1}",
        "s"    => L"$s$",
        "wl"   => L"$w_\mathrm{l}$",
        "wo"   => L"$w_\mathrm{o}$]",
        "L1_km"   => L"$L_1$",
        "L2_km"   => L"$L_2$",
        "L3_km"   => L"$L_3$",
        "u_ms"    => L"$u$",
        "w0"   => L"$w_0$",
        "PR12" => L"$PR_{12}$",
        "P1"   => L"$P_1$",
        "P2"   => L"$P_2$",
        "P3"   => L"$P_3$",
        "Δwtot"=> L"$Δw_\mathrm{tot}$",
    )
    return d
end

function labels_dict()
    d = Dict{String, LaTeXString}(
        "spwp" => L"$s_\mathrm{pwp}$",
        "ep"   => L"$e_\mathrm{p}$ [mm/day]", 
        "eo"   => L"$e_\mathrm{o}$ [mm/day]",
        "ϵ"    => L"$\epsilon$",
        "r"    => L"$r$",  
        "α"    => L"$\alpha$", 
        "nZr"  => L"$nZ_\mathrm{r}$ [mm]",   
        "a"    => L"$a$",  
        "b"    => L"$b$",
        "wsat" => L"$w_\mathrm{sat}$ [mm]", 
        "u"    => L"$u$ [m/s]",
        "L"    => L"$L$ [km]",
        "Pl"   => L"$P_\mathrm{l}$ [mm/day]",
        "Po"   => L"$P_\mathrm{o}$ [mm/day]",
        "El"   => L"$E_\mathrm{l}$ [mm/day]",
        "R"    => L"$R$ [mm/day]",
        "PR"   => L"$\chi$",
        "Φ"    => L"Infiltration $\Phi$",
        "Ptot" => L"$P_\mathrm{mean}$ [mm/day]",
        "A"    => L"$A_\mathrm{l}$ [mm/day]",
        "B"    => L"$-A_\mathrm{o}$ [mm/day]", 
        "τ"    => L"$\tau$ [1/day]",
        "τ_time" => L"$\tau^{-1}$ [day]",
        "s"    => L"$s$",
        "wl"   => L"$w_\mathrm{l}$ [mm]",
        "wo"   => L"$w_\mathrm{o}$ [mm]",
        "L1_km"   => L"$L_1$ [km]",
        "L2_km"   => L"$L_2$ [km]",
        "L3_km"   => L"$L_3$ [km]",
        "u_ms"    => L"$u$ [m/s]",
        "w0"   => L"$w_0$ [mm]",
        "PR12" => L"$PR_{12}$",
        "P1"   => L"$P_1$ [mm/day]",
        "P2"   => L"$P_2$ [mm/day]",
        "P3"   => L"$P_3$ [mm/day]",
        "Δwtot"=> L"$Δw_\mathrm{tot} [mm]$",
    )
    return d
end

function labels_norm_dict()
    d = Dict{String, LaTeXString}(
        "spwp" => L"$s_\mathrm{pwp}$",
        "ep"   => L"$e_\mathrm{p}$", 
        "eo"   => L"$e_\mathrm{o}$",
        "ϵ"    => L"$\epsilon$",
        "r"    => L"$r$",  
        "α"    => L"$\alpha$", 
        "nZr"  => L"$nZ_\mathrm{r}$",   
        "a"    => L"$a$",  
        "b"    => L"$b$",
        "wsat" => L"$w_\mathrm{sat}$", 
        "u"    => L"$u$ [m/s]",
        "L"    => L"$L$ [km]",
        "Pl"   => L"$P_\mathrm{l}$",
        "Po"   => L"$P_\mathrm{o}$",
        "El"   => L"$E_\mathrm{l}$",
        "R"    => L"$R$",
        "PR"   => L"$\chi$",
        "Ptot" => L"$P_\mathrm{mean}$",
        "A"    => L"$A_\mathrm{l}$",
        "B"    => L"$-A_\mathrm{o}$", 
        "τ"    => L"$\tau$",
        "τ_time" => L"$\tau^{-1}$",
        "s"    => L"$s$",
        "wl"   => L"$w_\mathrm{l}$",
        "wo"   => L"$w_\mathrm{o}$",
        "L1_km"   => L"$L_1$",
        "L2_km"   => L"$L_2$",
        "L3_km"   => L"$L_3$",
        "u_ms"    => L"$u$",
        "w0"   => L"$w_0$",
        "PR12" => L"$PR_{12}$",
        "P1"   => L"$P_1$",
        "P2"   => L"$P_2$",
        "P3"   => L"$P_3$",
        "Δwtot"=> L"$Δw_\mathrm{tot}$",
    )
    return d
end


function full_labels_dict()
    d = Dict{String, LaTeXString}(
        "spwp" => L"Permanent wilting point $s_\mathrm{pwp}$",
        "α"    => L"Land fraction $\alpha$", 
        "r"    => L"Runoff exponent $r$", 
        "ϵ"    => L"Runoff parameter $\epsilon$",
        "a"    => L"Precipitation parameter $a$",
        "b"    => L"Precipitation parameter $b$",
        "wsat" => L"Saturation water vapor pass $w_\mathrm{sat}$ [mm]",
        "eo"   => L"Ocean evaporation rate $E_\mathrm{o}$ [mm/day]",
        "PR"   => L"Precipitation ratio $\chi$",
        "τ"    => L"Atmospheric transport parameter $\tau$ [1/day]",
        "u_ms" => L"Wind speed $u$ [m/s]",
        "τ_time" => L"Timescale of atmospheric transport $\tau$ [day]",
        "Atot" => L"Advection rate $A\, \,  \left[10^8\,\mathrm{mm}^2\mathrm{/day}\right]$",
        "s"    => L"Soil moisture saturation $s$",
        "El"   => L"Evapotranspiration $E_\ell$ [mm/day]",
        "L_km" => L"Full domain length $L$ [km]",
        "L1"   => L"First ocean length $L_1$ [km]",
        "L2"   => L"Land length $L_2$ [km]",
        "L3"   => L"Second ocean length $L_3$ [km]",
        "w0"   => L"Boundary water vapor pass $w_0$ [mm]",
        "PR12" => L"Precipitation ratio $P_1/P_2$",
        "Δwtot"=> L"Total advection $Δw_\mathrm{tot}$ [mm]",
        "dw"   => L"Moisture difference $w_\mathrm{o}-w_\ell$ [mm]",
        "w1"   => L"Mean water vapor pass $w_{\mathrm{o},1}$ [mm]",
        "w2"   => L"Mean water vapor pass $w_{\ell}$ [mm]",
        "w3"   => L"Mean water vapor pass $w_{\mathrm{o},2}$ [mm]",
        "mean_wo"=> L"Mean ocean water vapor pass $w_\mathrm{o}$",
        "P1"   => L"Precipitation rate $P_1$ [mm/day]",
        "P2"   => L"Precipitation rate $P_2$ [mm/day]",
        "P3"   => L"Precipitation rate $P_3$ [mm/day]",
        "Po"   => L"Ocean precipitation rate $P_\mathrm{o}$ [mm/day]",
        "Pl"   => L"Land precipitation rate $P_\ell$ [mm/day]",
        "R"    => L"Runoff rate $R$ [mm/day]",
        "Φ"    => L"Infiltration function $\Phi$",
        "wl"   => L"Mean land water vapor pass $w_\ell$ [mm]",
        "wo"   => L"Mean ocean water vapor pass $w_\mathrm{o}$ [mm]",
        "ep"   => L"Potential ET $e_\mathrm{p}$ [mm/day]",
        "shat" => L"$\hat{s}$ indepdent of $s_\mathrm{pwp}$ and $s_\mathrm{fc}$",
        "Ehat" => L"Evapotranspiration $E$ [mm/day]",
        "A"    => L"Land advection rate $A$ [mm/day]",
        "B"    => L"Ocean advection rate $B$ [mm/day]",
    )
    return d
end


function titles_dict()
    d = Dict{String, String}(
        "spwp" => "Permanent wilting point",
        "ep"   => "Potential evapotranspiration", 
        "eo"   => "Mean ocean evaporation rate",
        "ϵ"    => "Runoff parameter",
        "r"    => "Runoff parameter",  
        "α"    => "Land fraction", 
        "nZr"  => "Reservoir depth",   
        "a"    => "Precipitation parameter",  
        "b"    => "Precipitation parameter",
        "wsat" => "Saturation water vapor pass", 
        "u"    => "Wind speed",
        "L"    => "Full domain size",
        "Pl"   => "Land precipitation rate",
        "Po"   => "Ocean precipitation rate",
        "El"   => "Land evapotranspiration rate",
        "R"    => "Runoff rate",
        "PR"   => "Precipitation ratio",
        "Ptot" => "Mean precipitation rate",
        "A"    => "Land advection rate",
        "B"    => "Ocean advection rate", 
        "Φ"  => "Infiltration function",
        "s"    => "Rel. soil moisture saturation",
        "wl"   => "Mean land water vapour pass",
        "wo"   => "Mean ocean water vapour pass",
        "w0"   => "Boundary water vapor pass",
        "L1_km"=> "First ocean length",
        "L2_km"=> "Land length",
        "L3_km"=> "Second ocean length",
        "PR12" => "PR without second ocean",
        "P1"   => "1st ocean precipitation rate",
        "P2"   => "Land precipitation rate",
        "P3"   => "2nd ocean precipitation rate",
        "Δwtot"=> "Total advection",
    )
    # d = Dict{String, LaTeXString}(
    #     "spwp" => LaTeXString("Permanent wilting point"),
    #     "ep"   => LaTeXString("Potential evapotranspiration"), 
    #     "eo"   => LaTeXString("Mean ocean evaporation rate"),
    #     "ϵ"    => LaTeXString("Runoff parameter"),
    #     "r"    => LaTeXString("Runoff parameter"),  
    #     "α"    => LaTeXString("Land fraction"), 
    #     "nZr"  => LaTeXString("Reservoir depth"),   
    #     "a"    => LaTeXString("Precipitation parameter"),  
    #     "b"    => LaTeXString("Precipitation parameter"),
    #     "wsat" => LaTeXString("Saturation water vapor pass"), 
    #     "u"    => LaTeXString("Wind speed"),
    #     "L"    => LaTeXString("Full domain size"),
    #     "Pl"   => LaTeXString("Land precipitation rate"),
    #     "Po"   => LaTeXString("Ocean precipitation rate"),
    #     "El"   => LaTeXString("Land evapotranspiration rate"),
    #     "R"    => LaTeXString("Runoff rate"),
    #     "PR"   => LaTeXString("Precipitation ratio"),
    #     "Ptot" => LaTeXString("Mean precipitation rate"),
    #     "A"    => LaTeXString("Land advection rate"),
    #     "B"    => LaTeXString("Ocean advection rate"), 
    # )
    return d
end



