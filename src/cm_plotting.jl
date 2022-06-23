
dcm     = CSV.read(datadir("sims", "closed model pmscan/final/cm_smooth_tau_fixedpoints_runs100000_updated_ranges_final_all_quantities" * ".csv"), DataFrame)
#dcmold   = CSV.read(datadir("sims", "closed model pmscan/cm_smooth_tau_eq_MC_fixedpoints_runs50000_updated_ranges_all_quantities" * ".csv"), DataFrame)
cm_mi   = CSV.read(datadir("sims", "mutual information/final/cm_rel_mi_100000_runs_final" * ".csv"), DataFrame)
#dcmoldold  = CSV.read(datadir("sims", "closed model pmscan/cm_smooth_tau_eq_MC_fixedpoints_runs10000_all_quantities" * ".csv"), DataFrame)
#cmpsens = CSV.read(datadir("sims", "closed model pmscan/cm_tau_10000runs_parameter_sensitivities_MI" * ".csv"), DataFrame)

function cm_t_evolution_plot(s, wl, wo, t, 
                            fs = 14.0,
                            c1 = "dodgerblue", 
                            c2 = "darkblue",
                            c3 = "darkgreen")
    fig = Figure(resolution = (1000, 600))
    ax1 = Axis(fig[1,1], xlabelsize = 28, ylabelsize = 28, yticklabelsize = 20, xticklabelsize = 20)
    ax2 = Axis(fig[2,1], xlabelsize = 28, ylabelsize = 28, yticklabelsize = 20, xticklabelsize = 20)
    ax2.xlabel = "t [days]"
    ax1.ylabel = "w [mm]"
    ax2.ylabel = "s"
    #ax2.ylabelsize = fs
    #ax2.xlabelsize = fs
    #ax2.ylabelsize = fs
    hidespines!(ax1, :t, :r)
    hidespines!(ax2, :t, :r)
    lines!(ax1, t, wl, linewidth = 3.0, label="w_l", color=c1)
    lines!(ax1, t, wo, linewidth = 3.0, label="w_o", color=c2)
    lines!(ax2, t, s, linewidth = 3.0, color=c3)
    axislegend(ax1, framevisible = false, labelsize = 28)
    return fig
end


function cm_basins_plot(b = basins, s = sg, wl = wlg)
    #cmap = colormap("Blues"; logscale=false)
    fig = Figure()
    ax = Axis(fig[1,1])
    hm = heatmap!(ax, s, wl, b)
    Colorbar(fig[1, 2], hm)
    return fig
end


function cm_parameter_scatter_plots(df::DataFrame, system, idtext, nb_bins = 100)
    ms = 5.0

    x = ["spwp"]#, "sfc", "ϵ", "τ", "wsat", "a", "b", "ep", "nZr", "α", "eo"]
    x_units = [""]#, "", "", " [1/day]", " [mm]", "", "", " [mm/day]", " [mm]", "", " [mm/day]"]
    #x = ["u"]
    #x_units = [" [m/s]"]
    y = ["PR", "Pl", "Po", "El", "Φ", "R"] 
    y_units = ["", " [mm/day]", " [mm/day]", " [mm/day]", "", " [mm/day]"]
    
    for n = 1:length(x)
   
        fig = Figure(resolution = (1200, 1500))
        ax1 = Axis(fig[1, 1])
        ax2 = Axis(fig[1, 2])
        ax3 = Axis(fig[2, 1])
        ax4 = Axis(fig[2, 2])
        ax5 = Axis(fig[3, 1])
        ax6 = Axis(fig[3, 2])
        
        scatter!(ax1, df[!,x[n]], df[!,y[1]], markersize = ms)
        lines!(ax1, sort(df, x[n])[!,x[n]], cm_mean_of_bins!(df,x[n],y[1],nb_bins), color = :red)
        #lines!(ax1, sort(df, x[n])[!,x[n]], cm_mean_of_bins!(df,x[n],y[2],nb_bins) ./ cm_mean_of_bins!(df, x[n], y[3], nb_bins), color = :red)
        #hlines!(ax1, mean(df[:,y[1]]), color = :orange)
        scatter!(ax2, df[!,x[n]], df[!,y[2]], markersize = ms)
        lines!(ax2, sort(df, x[n])[!,x[n]], cm_mean_of_bins!(df,x[n],y[2],nb_bins), color = :red)
        #hlines!(ax2, mean(df[:,y[2]]), color = :orange)
        scatter!(ax3, df[!,x[n]], df[!,y[3]], markersize = ms)
        lines!(ax3, sort(df, x[n])[!,x[n]], cm_mean_of_bins!(df,x[n],y[3],nb_bins), color = :red)
        #hlines!(ax3, mean(df[:,y[3]]), color = :orange)
        scatter!(ax4, df[!,x[n]], df[!,y[4]], markersize = ms)
        lines!(ax4, sort(df, x[n])[!,x[n]], cm_mean_of_bins!(df,x[n],y[4],nb_bins), color = :red)
        #hlines!(ax4, mean(df[:,y[4]]), color = :orange)
        scatter!(ax5, df[!,x[n]], df[!,y[5]], markersize = ms)
        lines!(ax5, sort(df, x[n])[!,x[n]], cm_mean_of_bins!(df,x[n],y[5],nb_bins), color = :red)
        #hlines!(ax5, mean(df[:,y[5]]), color = :orange)
        scatter!(ax6, df[!,x[n]], df[!,y[6]], markersize = ms)
        lines!(ax6, sort(df, x[n])[!,x[n]], cm_mean_of_bins!(df,x[n],y[6],nb_bins), color = :red)
        #hlines!(ax6, mean(df[:,y[6]]), color = :orange)
        
        ax5.xlabel = string(x[n], x_units[n])
        ax6.xlabel = string(x[n], x_units[n])
        
        ax1.ylabel = string(y[1],y_units[1])
        ax2.ylabel = string(y[2],y_units[2])
        ax3.ylabel = string(y[3],y_units[3])
        ax4.ylabel = string(y[4],y_units[4])
        ax5.ylabel = string(y[5],y_units[5])
        ax6.ylabel = string(y[6],y_units[6])
                
        hidespines!(ax1, :t, :r)
        hidespines!(ax2, :t, :r)
        hidespines!(ax3, :t, :r)
        hidespines!(ax4, :t, :r)
        hidespines!(ax5, :t, :r)
        hidespines!(ax6, :t, :r)

        #supertitle = fig[0, :] = Label(fig, "Domain size = $(round(df[1, "L"], digits = 0)) km", textsize = 22)

        save(plotsdir("Closed model","Parameter scatter plots/", system * "/cm_s_10000runs_" * idtext * "_$(x[n])_scatter.png"),fig)
        #return fig
    end
end

function cm_param_vars_scatter_plots(df::DataFrame, system, idtext, nb_bins = 100)
    ms = 5.0

    x = ["spwp", "ϵ", "u", "wsat", "a", "α", "eo"]
    x_units = ["", "", " [m/s]", " [mm]", "", "", " [mm/day]"]
    y = ["s", "wl", "wo"] 
    y_units = ["", " [mm]", " [mm]"]
    
    for n = 1:length(x)
   
        fig = Figure(resolution = (800, 1000))
        ax1 = Axis(fig[1, :], title = "Domain size = $(round(df[1, "L"], digits=0)) km", textsize = 22)
        ax2 = Axis(fig[2, :])
        ax3 = Axis(fig[3, :])
        
        scatter!(ax1, df[!,x[n]], df[!,y[1]], markersize = ms)
        lines!(ax1, sort(df, x[n])[!,x[n]], cm_mean_of_bins!(df,x[n],y[1],nb_bins), color = :red)
        hlines!(ax1, mean(df[:,y[1]]), color = :orange)
        scatter!(ax2, df[!,x[n]], df[!,y[2]], markersize = ms)
        lines!(ax2, sort(df, x[n])[!,x[n]], cm_mean_of_bins!(df,x[n],y[2],nb_bins), color = :red)
        hlines!(ax2, mean(df[:,y[2]]), color = :orange)
        scatter!(ax3, df[!,x[n]], df[!,y[3]], markersize = ms)
        lines!(ax3, sort(df, x[n])[!,x[n]], cm_mean_of_bins!(df,x[n],y[3],nb_bins), color = :red)
        hlines!(ax3, mean(df[:,y[3]]), color = :orange)
        
        ax3.xlabel = string(x[n], x_units[n])
        
        ax1.ylabel = string(y[1],y_units[1])
        ax2.ylabel = string(y[2],y_units[2])
        ax3.ylabel = string(y[3],y_units[3])

        hidespines!(ax1, :t, :r)
        hidespines!(ax2, :t, :r)
        hidespines!(ax3, :t, :r)

        #supertitle = fig[0, :] = Label(fig, , textsize = 22)
        save(plotsdir("Closed model","State variables scatter plots", system * "/cm_10000runs_" * idtext * "_$(x[n])_scatter.png"),fig)
        #return fig
    end
end

function cm_domainsize_influence_scatter_plots(df1::DataFrame, df2::DataFrame, system, nb_bins = 100)
    ms = 5.0
    lw = 3.0

    x = ["spwp", "ϵ", "u", "wsat", "a", "b", "ep", "nZr", "α", "eo"]
    x_units = ["", "", " [m/s]", " [mm]", "", "", " [mm/day]", " [mm]", "", " [mm/day]"]
    y = ["PR", "Pl", "Po", "El", "infilt", "runoff"] 
    y_units = ["", " [mm/day]", " [mm/day]", " [mm/day]", "", " [mm/day]"]
    
    for n = 1:length(x)
   
        fig = Figure(resolution = (1200, 1500))
        ax1 = Axis(fig[1, 1])
        ax2 = Axis(fig[1, 2])
        ax3 = Axis(fig[2, 1])
        ax4 = Axis(fig[2, 2])
        ax5 = Axis(fig[3, 1])
        ax6 = Axis(fig[3, 2])
        
        scatter!(ax1, df1[!,x[n]], df1[!,y[1]], color = (:royalblue, 0.6), markersize = ms)
        scatter!(ax1, df2[!,x[n]], df2[!,y[1]], color = (:chartreuse4, 0.6), markersize = ms)
        lines!(ax1, sort(df1, x[n])[!,x[n]], cm_mean_of_bins!(df1,x[n],y[1],nb_bins), linewidth = lw, color = :midnightblue)
        #lines!(ax1, sort(df1, x[n])[!,x[n]], cm_mean_of_bins!(df1,x[n],y[2],nb_bins) ./ cm_mean_of_bins!(df1, x[n], y[3], nb_bins), linewidth = lw, color = :midnightblue)
        hlines!(ax1, mean(df1[:,y[1]]), linewidth = lw, color = :midnightblue)
        lines!(ax1, sort(df2, x[n])[!,x[n]], cm_mean_of_bins!(df2,x[n],y[1],nb_bins), linewidth = lw, color = :darkgreen)
        #lines!(ax1, sort(df2, x[n])[!,x[n]], cm_mean_of_bins!(df2,x[n],y[2],nb_bins) ./cm_mean_of_bins!(df2,x[n],y[3],nb_bins), linewidth = lw, color = :darkgreen)
        hlines!(ax1, mean(df2[:,y[1]]), linewidth = lw, color = :darkgreen)

        scatter!(ax2, df1[!,x[n]], df1[!,y[2]], color = (:royalblue, 0.6), markersize = ms)
        scatter!(ax2, df2[!,x[n]], df2[!,y[2]], color = (:chartreuse4, 0.6), markersize = ms)
        lines!(ax2, sort(df1, x[n])[!,x[n]], cm_mean_of_bins!(df1,x[n],y[2],nb_bins), linewidth = lw, color = :midnightblue)
        hlines!(ax2, mean(df1[:,y[2]]), linewidth = lw, color = :midnightblue)
        lines!(ax2, sort(df2, x[n])[!,x[n]], cm_mean_of_bins!(df2,x[n],y[2],nb_bins), linewidth = lw, color = :darkgreen)
        hlines!(ax2, mean(df2[:,y[2]]), linewidth = lw, color = :darkgreen)

        scatter!(ax3, df1[!,x[n]], df1[!,y[3]], color = (:royalblue, 0.6), markersize = ms)
        scatter!(ax3, df2[!,x[n]], df2[!,y[3]], color = (:chartreuse4, 0.6), markersize = ms)
        lines!(ax3, sort(df1, x[n])[!,x[n]], cm_mean_of_bins!(df1,x[n],y[3],nb_bins), linewidth = lw, color = :midnightblue)
        hlines!(ax3, mean(df1[:,y[3]]), linewidth = lw, color = :midnightblue)
        lines!(ax3, sort(df2, x[n])[!,x[n]], cm_mean_of_bins!(df2,x[n],y[3],nb_bins), linewidth = lw, color = :darkgreen)
        hlines!(ax3, mean(df2[:,y[3]]), linewidth = lw, color = :darkgreen)

        scatter!(ax4, df1[!,x[n]], df1[!,y[4]], color = (:royalblue, 0.6), markersize = ms)
        scatter!(ax4, df2[!,x[n]], df2[!,y[4]], color = (:chartreuse4, 0.6), markersize = ms)
        lines!(ax4, sort(df1, x[n])[!,x[n]], cm_mean_of_bins!(df1,x[n],y[4],nb_bins), linewidth = lw, color = :midnightblue)
        hlines!(ax4, mean(df1[:,y[4]]), linewidth = lw, color = :midnightblue)
        lines!(ax4, sort(df2, x[n])[!,x[n]], cm_mean_of_bins!(df2,x[n],y[4],nb_bins), linewidth = lw, color = :darkgreen)
        hlines!(ax4, mean(df2[:,y[4]]), linewidth = lw, color = :darkgreen)

        scatter!(ax5, df1[!,x[n]], df1[!,y[5]], color = (:royalblue, 0.6), markersize = ms)
        scatter!(ax5, df2[!,x[n]], df2[!,y[5]], color = (:chartreuse4, 0.6), markersize = ms)
        lines!(ax5, sort(df1, x[n])[!,x[n]], cm_mean_of_bins!(df1,x[n],y[5],nb_bins), linewidth = lw, color = :midnightblue)
        hlines!(ax5, mean(df1[:,y[5]]), linewidth = lw, color = :midnightblue)
        lines!(ax5, sort(df2, x[n])[!,x[n]], cm_mean_of_bins!(df2,x[n],y[5],nb_bins), linewidth = lw, color = :darkgreen)
        hlines!(ax5, mean(df2[:,y[5]]), linewidth = lw, color = :darkgreen)

        scatter!(ax6, df1[!,x[n]], df1[!,y[6]], color = (:royalblue, 0.6), markersize = ms)
        scatter!(ax6, df2[!,x[n]], df2[!,y[6]], color = (:chartreuse4, 0.6), markersize = ms)
        lines!(ax6, sort(df1, x[n])[!,x[n]], cm_mean_of_bins!(df1,x[n],y[6],nb_bins), linewidth = lw, color = :midnightblue)
        hlines!(ax6, mean(df1[:,y[6]]), linewidth = lw, color = :midnightblue)
        lines!(ax6, sort(df2, x[n])[!,x[n]], cm_mean_of_bins!(df2,x[n],y[6],nb_bins), linewidth = lw, color = :darkgreen)
        hlines!(ax6, mean(df2[:,y[6]]), linewidth = lw, color = :darkgreen)
        
        ax5.xlabel = string(x[n], x_units[n])
        ax6.xlabel = string(x[n], x_units[n])
        
        ax1.ylabel = string(y[1],y_units[1])
        ax2.ylabel = string(y[2],y_units[2])
        ax3.ylabel = string(y[3],y_units[3])
        ax4.ylabel = string(y[4],y_units[4])
        ax5.ylabel = string(y[5],y_units[5])
        ax6.ylabel = string(y[6],y_units[6])
                
        hidespines!(ax1, :t, :r)
        hidespines!(ax2, :t, :r)
        hidespines!(ax3, :t, :r)
        hidespines!(ax4, :t, :r)
        hidespines!(ax5, :t, :r)
        hidespines!(ax6, :t, :r)

        #supertitle = fig[0, :] = Label(fig, "Domain size = $(round(df[1, "L"], digits = 0)) km", textsize = 22)

        save(plotsdir("Closed model","Parameter scatter plots/", system * "/cm_10000runs_L1000+L10000_$(x[n])_scatter.png"),fig)
        #return fig
    end
end

function simple_scatter_plot(system, domain, df1, df2, x, y, xlabel, ylabel)
    fig = Figure()
    ax  = Axis(fig[1,1])
    hidespines!(ax, :t, :r)
    scatter!(ax, df1[!,x], df1[!,y], color = (:royalblue, 0.6), markersize = 5.0)
    #lines!(ax, sort(df1, x)[!,x], cm_mean_of_bins!(df1,x,y,100), color = :black)
    #scatter!(ax, df2[!,x], df2[!,y], color = (:chartreuse4, 0.6), markersize = 5.0)
    ax.xlabel = xlabel
    ax.ylabel = ylabel
    #save(plotsdir("Closed model", "Parameter scatter plots/", system * "/cm_10000runs_domain" * domain * "_" * x * "_" * y * ".png"),fig)
    return fig
end

function PR_alpha_comp_plot(df1, df2, lb1, lb2, x, y, xlab, ylab, ttl, comp)
    fig = Figure()
    ax  = Axis(fig[1,1], xlabel = xlab, ylabel = ylab, title = ttl)
    hidespines!(ax, :t, :r)
    scatter!(ax, df1[!,x], df1[!,y], color = (:royalblue, 0.6), label = lb1, markersize = 5.0)
    #scatter!(ax, df2[!,x], df2[!,y], color = (:chartreuse4, 0.6), label = lb2, markersize = 5.0)
    scatter!(ax, df2[!,x], df2[!,y], color = (:chocolate, 0.6), label = lb2, markersize = 5.0)
    vlines!(ax, df1[df1.PR .== minimum(df1.PR), "α"][1], linewidth = 3.0, linestyle = :dash, color = (:royalblue, 0.6))
    vlines!(ax, df2[df2.PR .== minimum(df2.PR), "α"][1], linewidth = 3.0, linestyle = :dash, color = (:chocolate, 0.6))
    axislegend(ax, position = :lb)
    save(plotsdir("Closed model", "Parameter scatter plots/", "smooth","PR-alpha_parameter_influences", "PR-alpha_$(comp)_varied.png"),fig)
    return fig
end



function PR_u_plots()
    nb_bins = 10
    pw1000lw = sort(pw1000, "wl")[1:1000, :]
    pw10000lw = sort(pw10000, "wl")[1:1000, :]
    pw1000hw = sort(pw1000, "wl")[9001:10000, :]
    pw10000hw = sort(pw10000, "wl")[9001:10000, :]
    fig = Figure(resolution = (1000, 600))
    ax1 = Axis(fig[1,1])
    ax2 = Axis(fig[1,2])
    hidespines!(ax1, :t, :r)
    hidespines!(ax2, :t, :r)
    scatter!(ax1, pw1000lw[!,"u"], pw1000lw[!,"PR"], color = (:lightseagreen, 0.6), markersize = 5.0)
    scatter!(ax1, pw1000hw[!,"u"], pw1000hw[!,"PR"], color = (:royalblue, 0.6), markersize = 5.0)
    lines!(ax1, sort(pw1000lw, "u")[!,"u"], cm_mean_of_bins!(pw1000lw,"u","PR",nb_bins), label = "low wl-values, dw_mean = 0.246", linewidth = 3.0, color = :lightseagreen)
    lines!(ax1, sort(pw1000hw, "u")[!,"u"], cm_mean_of_bins!(pw1000hw,"u","PR",nb_bins), label = "high wl-values, dw_mean = 0.251", linewidth = 3.0, color = :royalblue)
    axislegend(ax1, position = :rb, frame_visible = false)
    scatter!(ax2, pw10000lw[!,"u"], pw10000lw[!,"PR"], color = (:greenyellow, 0.6), markersize = 5.0)
    scatter!(ax2, pw10000hw[!,"u"], pw10000hw[!,"PR"], color = (:chartreuse4, 0.6), markersize = 5.0)
    lines!(ax2, sort(pw10000lw, "u")[!,"u"], cm_mean_of_bins!(pw10000lw,"u","PR",nb_bins), label = "low wl-values, dw_mean = 2.136", linewidth = 3.0, color = :greenyellow)
    lines!(ax2, sort(pw10000hw, "u")[!,"u"], cm_mean_of_bins!(pw10000hw,"u","PR",nb_bins), label = "high wl-values, dw_mean = 1.486", linewidth = 3.0, color = :chartreuse4)
    axislegend(ax2, position = :rb, frame_visible = false)

    ax1.xlabel = "u [m/s]"
    ax2.xlabel = "u [m/s]"
    ax1.ylabel = "PR"
    save(plotsdir("Closed model", "Parameter scatter plots/", "piecewise/cm_10000runs_domain_both_PR-u-regimes.png"),fig)
end

function two_tiles_plot(data::DataFrame, param::String, quant1::String, quant2::String)
    l = labels_dict()
    lfs = 20
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

function rel_mi_plot(rel_mi_data::DataFrame)
    sort!(rel_mi_data, "MI_rel", rev = true)
    l = short_labels_dict()
    n = nrow(rel_mi_data)

    #rel_mi_data.labels = [l[i] for i in rel_mi_data.pnames]
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

function four_tiles_plot(data::DataFrame, yquant::String, xquant1::String, xquant2::String, xquant3::String, xquant4::String)
    l = labels_dict()
    lfs = 20
    ms = 5
    nb_bins = 100
    c1 = (:grey35, 0.5)
    c2 = :grey25
    fig = Figure(resolution = (800, 650))
    ax1 = Axis(fig[1, 1], xlabel = l[xquant1], ylabel = l[yquant], xlabelsize = lfs, ylabelsize = lfs, xgridcolor = :white, ygridcolor = :white)
    text!(L"\mathrm{a})", position = (0.6, 0.2))
    
    ax2 = Axis(fig[1, 2], xlabel = l[xquant2], xlabelsize = lfs, xgridcolor = :white, ygridcolor = :white)
    text!(L"\mathrm{b})", position = (0, 0.1))
   
    ax3 = Axis(fig[2, 1], xlabel = l[xquant3], ylabel = l[yquant], xlabelsize = lfs, ylabelsize = lfs, xgridcolor = :white, ygridcolor = :white)
    text!(L"\mathrm{c})", position = (0.2, 0.1))
    
    ax4 = Axis(fig[2, 2], xlabel = l[xquant4], xlabelsize = lfs, ylabelsize = lfs, xgridcolor = :white, ygridcolor = :white)
    text!(L"\mathrm{d})", position = (5.5, 0.1))
    
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

function fluxes_plot(data::DataFrame, statevar::String, bin_length::Int = 20000, nb_bins::Int = 200)
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
    text!(L"\mathrm{a})", position = (0, 2.5))

    if type == "kde"
        # density!(data[!,"s"], color = (:darkorange, 0.8))
        # density!(data[!,"s"], color = (:grey30, 0.8), strokearound = true, strokewidth = 1) #(:grey20, 0.8))
        lines!(kde(data[!,"s"]), color = :grey30, linewidth = lw)
    elseif type == "hist"
        hist!(ax1, data[!,"s"], bins = 100, normalization = norm, color = (:darkorange, 0.8))
    end


    ax2 = Axis(f1[1,2], xlabel = L"\tilde{s} = \frac{s - s_\mathrm{pwp}}{s_\mathrm{sfc}-s_\mathrm{pwp}}", xlabelsize = lfs, xgridvisible = false, ygridvisible = false)
    data.sresc = (data.s .- data.spwp) ./ (data.sfc .- data.spwp)
    text!(L"\mathrm{b})", position = (-1, 4))

    if type == "kde"
        # density!(data[!,"sresc"], color = (:white, 0.8), strokearound = true, strokewidth = 1) # color = (:grey20, 0.8))
        lines!(kde(data[!,"sresc"]), color = :grey30, linewidth = lw)
        vlines!(ax2, 0.0, linewidth = 2.0, linestyle = :dash, color = :grey40)
        vlines!(ax2, 1.0, linewidth = 2.0, linestyle = :dash, color = :grey40)
    elseif type == "hist"
        hist!(ax2, data[!,"sresc"], bins = 100, normalization = norm, color = (:darkorange, 0.8))
    end

    ax3 = Axis(f1[1, 3], xlabel = L"Mean water vapor pass $w$", xlabelsize = lfs, xgridvisible = false, ygridvisible = false)
    text!(L"\mathrm{c})", position = (10, 0.05))

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

    ax4 = Axis(f1[1,4], xlabel = L"Precipitation ratio $\chi$", xlabelsize = lfs, xgridvisible = false, ygridvisible = false)
    text!(L"\mathrm{d})", position = (0.2, 17))

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
    gc = :white
    c = [:grey20, :grey60]
    s = collect(0.0:0.01:1.0)
    ep = 5.0
    spwp = [0.3, 0.4] 
    labls = [L"$s_\mathrm{pwp} = 0.3$", L"$s_\mathrm{pwp} = 0.4$"] 
    fig = Figure(resolution = (500, 400))
    ax = Axis(fig[1,1], xlabel = l["s"], ylabel = l["El"], ylabelsize = lfs, xlabelsize = lfs, xgridcolor = gc, ygridcolor = gc)

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