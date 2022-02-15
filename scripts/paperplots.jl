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

include(srcdir("cm_analysis.jl"))
include(srcdir("om_analysis.jl"))
include(srcdir("utils.jl"))

dcm  = CSV.read(datadir("sims", "closed model pmscan/cm_smooth_tau_eq_MC_fixedpoints_runs10000_all_quantities" * ".csv"), DataFrame)
cmpsens = CSV.read(datadir("sims", "closed model pmscan/cm_tau_10000runs_parameter_sensitivities_MI" * ".csv"), DataFrame)
#doms = CSV.read(datadir("sims", "open model pmscan/sym/om_v2_MC_fixedpoints_runs10000_sym_all_quantities" * ".csv"), DataFrame)
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
    ax = Axis(fig[1, 1], xlabel = l[xquant], ylabel = l[yquant], ylabelsize = lfs, xlabelsize = lfs)
    scatter!(ax, data[!,xquant], data[!,yquant], markersize = ms, color = c1)
    #hlines!(ax, 0.43, linewidth = 2.0, linestyle = :dot, color = :grey80)
    #hlines!(ax, 0.51, linewidth = 2.0, linestyle = :dot, color = :grey80)
    #vlines!(ax, 0.3, linewidth = 2.0, linestyle = :dot, color = :grey80)
    #vlines!(ax, 0.4, linewidth = 2.0, linestyle = :dot, color = :grey80)
    #lines!(ax, sort(data, xquant)[!,xquant], cm_mean_of_bins!(data, xquant, yquant, nb_bins), color = c2)
    #lines!(ax, sort(data, xquant)[!,xquant], cm_rolling_average!(data, xquant, yquant, 50), color = :forestgreen, label = "bin size 50")
    #lines!(ax, sort(data, xquant)[!,xquant], cm_rolling_average!(data, xquant, yquant, nb_bins), color = :blue, label = "bin size 1000")
    lines!(ax, sort!(data, xquant)[!,xquant], movingaverage(data[!,yquant], 100))
    #lines!(ax, data[!,xquant], data[!,xquant], linestyle = :dot, color = :orange)
    hidespines!(ax, :t, :r)
    #axislegend(ax, position = :lb, framevisible = false, labelsize = 16)
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
    ax1 = Axis(fig[1, 1], ylabel = l[quant1], ylabelsize = lfs, title = t[quant1], titlesize = tfs)
    ax2 = Axis(fig[2, 1], ylabel = l[quant2], ylabelsize = lfs, title = t[quant2], titlesize = tfs)
    ax3 = Axis(fig[3, 1], xlabel = l[param], ylabel = l[quant3], xlabelsize = lfs, ylabelsize = lfs, title = t[quant3], titlesize = tfs)
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
    ax1 = Axis(fig[1, 1], ylabel = l[quant1], xlabelsize = lfs, ylabelsize = lfs, title = t[quant1], titlesize = tfs)
    ax2 = Axis(fig[1, 2], ylabel = l[quant2], xlabelsize = lfs, ylabelsize = lfs, title = t[quant2], titlesize = tfs)
    ax3 = Axis(fig[2, 1], ylabel = l[quant3], xlabelsize = lfs, ylabelsize = lfs, title = t[quant3], titlesize = tfs)
    ax4 = Axis(fig[2, 2], ylabel = l[quant4], xlabelsize = lfs, ylabelsize = lfs, title = t[quant4], titlesize = tfs)
    ax5 = Axis(fig[3, 1], xlabel = l[param], ylabel = l[quant5], xlabelsize = lfs, ylabelsize = lfs, title = t[quant5], titlesize = tfs)
    ax6 = Axis(fig[3, 2], xlabel = l[param], ylabel = l[quant6], xlabelsize = lfs, ylabelsize = lfs, title = t[quant6], titlesize = tfs)
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

function rel_mi_plot(data::DataFrame)
    sort!(data, "MI_rel", rev = true)
    l = short_labels_dict()
    n = nrow(data)

    l_arr = Vector{String}()
    for i = 1:n
        push!(l_arr, l[data[i,"pnames"]])
    end
    
    xrange = collect(1:1:n)
    f = Figure(resolution = (700, 550))
    ax = Axis(f[1,1], yscale = log10, xlabel = L"Model parameters $p_i$", ylabel = L"Relative mututal information between $PR$ and $p_i$", xlabelsize = 20, ylabelsize = 20, xgridcolor = :white, ygridcolor = :white)
    ylims!(ax, 0.1, 100.0) # separate
    hlines!(ax, 1.0, linestyle = :dash, color = :grey35, label = "3σ significance threshold")
    scatter!(ax, xrange, data[:, "MI_rel"], marker = '*', markersize = 25, color = :grey20)
    ax.xticks = (1:1:n, data[:,"pnames"])
    axislegend(ax, framevisible = false, labelsize = 16)
    hidespines!(ax, :t, :r)
    return f
end

function fluxes_plot(data::DataFrame, statevar::String, bin_length::Int = 100, nb_bins::Int = 200)
    fl = full_labels_dict()
    lfs = 20
    legendfs = 16
    s = collect(minimum(dcm.s):0.01:maximum(dcm.s))
    w = collect(minimum(dcm.wo):1.0:maximum(dcm.wo))
    colors = [colorant"#ff7f00", colorant"#33a02c", colorant"#6a3d9a",
              colorant"#1f78b4", colorant"#a6cee3", colorant"#e31a1c"]

    fig = Figure(resolution = (600, 400))
    ax = Axis(fig[1, 1], xlabel = fl[statevar], title = "bin length: $(bin_length)", ylabel = L"Mean water fluxes $F_i$ [mm/day]", ylabelsize = lfs, xlabelsize = lfs,
                xgridcolor = :white, ygridcolor = :white)
    fluxes = ["eo", "El", "R", "Po", "Pl", "A"]
    ls = [:solid, :solid, :solid, :solid, :solid, :dot]
    band!(ax, w, [exp(15.6*(el / 65.0 - 0.522)) for el in w], [exp(11.4*(el / 80.0 - 0.603)) for el in w], color = (:lightblue1, 0.8))
    scatter!(ax, data.wo, data.Po, color = (colorant"#1f78b4", 0.6), markersize = 2.0)
    #band!(ax, s, [(4.1/2 * tanh( 10 * (el - 1.38/2 ) ) + 4.1/2) for el in s], [(4.5/2 * tanh( 10 * (el - 0.7/2 ) ) + 4.5/2) for el in s], color = (:lightgreen, 0.5))
    #scatter!(ax, data.s, data.El, color = (colorant"#33a02c", 0.4), markersize = 2.0)
    
    for i in 1:length(fluxes)
        #lines!(ax, sort(data, statevar)[!,statevar], cm_mean_of_bins!(data, statevar, fluxes[i], nb_bins), label = fluxes[i], color = colors[i], linestyle = ls[i])
        #lines!(ax, sort(data, statevar)[!,statevar], cm_rolling_average!(data, statevar, fluxes[i], nb_bins), label = fluxes[i], color = colors[i], linestyle = ls[i])
        lines!(ax, sort!(data, statevar)[!,statevar], movingaverage(data[:,fluxes[i]], bin_length), label = fluxes[i], color = colors[i], linestyle = ls[i])
        #scatter!(ax, dcm[:, statevar], dcm[:,fluxes[i]], color = (colors[i], 0.5), markersize = 2.0)
    end
    #lines!(ax, s, [(4.5/2 * tanh( 10 * (el - 0.7/2 ) ) + 4.5/2) for el in s], linestyle = :dot, color = :green)
    #lines!(ax, s, [(4.1/2 * tanh( 10 * (el - 1.38/2 ) ) + 4.1/2) for el in s], linestyle = :dot, color = :green)
    #lines!(ax, w, [exp(15.6*(x / 65.0 - 0.522)) for el in w], linestyle = :dot, color = :blue)
    #lines!(ax, w, [exp(11.4*(el / 80.0 - 0.603)) for el in w], linestyle = :dot, color = :blue)
    axislegend(ax, position = :lc, framevisible = false, labelsize = legendfs)
    hidespines!(ax, :t, :r)
    ylims!(ax, (0,3.5))
    return fig
end


function El_plot()
    l = full_labels_dict()
    lfs = 20
    legendfs = 16
    lw = 3.0
    c = [:grey20, :grey40, :grey60, :grey70]
    s = collect(0.0:0.01:1.0)
    ep = 4.3
    spwp = [0.2, 0.3, 0.4, 0.5]
    labls = [L"$s_\mathrm{pwp} = 0.2$", L"$s_\mathrm{pwp} = 0.3$", L"$s_\mathrm{pwp} = 0.4$", L"$s_\mathrm{pwp} = 0.5$"]
    fig = Figure(resolution = (500, 400))
    ax = Axis(fig[1,1], xlabel = l["s"], ylabel = l["El"], ylabelsize = lfs, xlabelsize = lfs)

    for n=1:length(spwp)
        E = [(ep/2 * tanh( 10.0 * (el - (2*spwp[n]+0.3)/2 ) ) + ep/2) for el in s]
        lines!(ax, s, E, linewidth = lw, color = c[n], label = labls[n])
    end
    vlines!(ax, 0.43, linewidth = 2.0, linestyle = :dot, color = :grey80)
    vlines!(ax, 0.51, linewidth = 2.0, linestyle = :dot, color = :grey80)
    hidespines!(ax, :t, :r)
    axislegend(ax, position = :lt, framevisible = false, labelsize = legendfs)
    return fig
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
        "PR"   => L"$PR$",
        "Ptot" => L"$P_\mathrm{mean}$",
        "A"    => L"$A_\mathrm{l}$",
        "B"    => L"$-A_\mathrm{o}$", 
        "τ"    => "τ",
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
        "PR"   => L"$PR$",
        "Ptot" => L"$P_\mathrm{mean}$ [mm/day]",
        "A"    => L"$A_\mathrm{l}$ [mm/day]",
        "B"    => L"$-A_\mathrm{o}$ [mm/day]", 
        "τ"    => L"$\tau$ [1/day]",
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
        "PR"   => L"$PR$",
        "Ptot" => L"$P_\mathrm{mean}$",
        "A"    => L"$A_\mathrm{l}$",
        "B"    => L"$-A_\mathrm{o}$", 
        "τ"    => L"$\tau$",
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
        "PR"   => L"Precipitation ratio $PR$",
        "τ"    => L"Rate of transport $\tau$ [1/day]",
        "Atot" => L"Advection rate $A\, \,  \left[10^8\,\mathrm{mm}^2\mathrm{/day}\right]$",
        "s"    => L"Rel. soil moisture saturation $s$",
        "El"   => L"Evapotranspiration $E_\mathrm{l}$ [mm/day]",
        "L1"   => L"First ocean length $L_1$ [km]",
        "L2"   => L"Land length $L_2$ [km]",
        "L3"   => L"Second ocean length $L_3$ [km]",
        "w0"   => L"Boundary water vapor pass $w_0$ [mm]",
        "PR12" => L"Precipitation ratio $P_1/P_2$",
        "Δwtot"=> L"Total advection $Δw_\mathrm{tot}$ [mm]",
        "w1"   => L"Mean water vapor pass $w_1$ [mm]",
        "w2"   => L"Mean water vapor pass $w_2$ [mm]",
        "w3"   => L"Mean water vapor pass $w_3$ [mm]",
        "P1"   => L"Precipitation rate $P_1$ [mm/day]",
        "P2"   => L"Precipitation rate $P_2$ [mm/day]",
        "P3"   => L"Precipitation rate $P_3$ [mm/day]",
        "Pl"   => L"Land precipitation rate $P_l$ [mm/day]",
        "R"    => L"Runoff rate $R$ [mm/day]",
        "Φ"    => L"Infiltration function $\Phi$",
        "wl"   => L"Mean land water vapor pass $w_l$ [mm]",
        "wo"   => L"Mean ocean water vapor pass $w_o$ [mm]",
        "ep"   => L"Potential ET $e_p$ [mm/day]",
        "shat" => L"$\hat{s}$ indepdent of $s_{pwp}$ and $s_{fc}$",
        "Ehat" => L"Evapotranspiration $E$ [mm/day]",
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



