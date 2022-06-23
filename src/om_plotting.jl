#dom = CSV.read(datadir("sims", "open model pmscan/final/om_v2_sym_fixedpoints_runs100000_updated_ranges_final_all_quantities" * ".csv"), DataFrame)
#domold = CSV.read(datadir("sims", "open model pmscan/om_v2_fixedpoints_runs50000_sym_updatedparams_all_quantities" * ".csv"), DataFrame)
#domasl = CSV.read(datadir("sims", "open model pmscan/asym/om_v2_MC_fixedpoints_runs10000_asym3L1=L3_all_quantities" * ".csv"), DataFrame)
#domasr = CSV.read(datadir("sims", "open model pmscan/asym/om_v2_MC_fixedpoints_runs10000_asymL1=3L3_all_quantities" * ".csv"), DataFrame)

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

function eight_scatter_plots(df::DataFrame, nb_bins, land_dist)
    x = ["w0", "τ", "α"] #["spwp", "ep", "eo", "ϵ", "α", "nZr", "a", "b", "wsat", "w0", "u_ms", "τ"]
    x_units = [" [mm]", " [1/day]", ""] #["", " [mm/day]", " [mm/day]", "", "", " [mm]", "", "", " [mm]", " [mm]", " [m/s]", " [1/day]"]
    x_ranges = [(0, 70), (0, 1.8), (0, 1)]
    y = ["w1", "w2", "w3", "s", "Δw1", "Δw2", "Δw3", "Δw_net"]#["P1", "P2", "P3", "El", "R", "Φ", "Ptot", "PR"]
    y_units = [" [mm]", " [mm]", " [mm]", "", " [mm]", " [mm]", " [mm]", " [mm]"] #[" [mm/day]", " [mm/day]", " [mm/day]", " [mm/day]", " [mm/day]", " [mm]", " [mm]", ""]
    ms = 5.0
    nb_rows = nrow(df)

    for i = 1:length(x)
        fig = Figure(resolution = (1000, 1500))
        ax1 = Axis(fig[1, 1])
        ax2 = Axis(fig[1, 2])
        ax3 = Axis(fig[2, 1])
        ax4 = Axis(fig[2, 2])
        ax5 = Axis(fig[3, 1])
        ax6 = Axis(fig[3, 2])
        ax7 = Axis(fig[4, 1])
        ax8 = Axis(fig[4, 2])
        scatter!(ax1, df[!,x[i]], df[!,y[1]], markersize = ms)
        lines!(ax1, sort(df, x[i])[!,x[i]], mean_of_bins!(df,x[i],y[1],nb_bins), color = :red)
        hlines!(ax1, mean(df[:,y[1]]), color = :orange)
        scatter!(ax2, df[!,x[i]], df[!,y[2]], markersize = ms)
        lines!(ax2, sort(df, x[i])[!,x[i]], mean_of_bins!(df,x[i],y[2],nb_bins), color = :red)
        hlines!(ax2, mean(df[:,y[2]]), color = :orange)
        scatter!(ax3, df[!,x[i]], df[!,y[3]], markersize = ms)
        lines!(ax3, sort(df, x[i])[!,x[i]], mean_of_bins!(df,x[i],y[3],nb_bins), color = :red)
        hlines!(ax3, mean(df[:,y[3]]), color = :orange)
        scatter!(ax4, df[!,x[i]], df[!,y[4]], markersize = ms)
        lines!(ax4, sort(df, x[i])[!,x[i]], mean_of_bins!(df,x[i],y[4],nb_bins), color = :red)
        hlines!(ax4, mean(df[:,y[4]]), color = :orange)
        scatter!(ax5, df[!,x[i]], df[!,y[5]], markersize = ms)
        lines!(ax5, sort(df, x[i])[!,x[i]], mean_of_bins!(df,x[i],y[5],nb_bins), color = :red)
        hlines!(ax5, mean(df[:,y[5]]), color = :orange)
        scatter!(ax6, df[!,x[i]], df[!,y[6]], markersize = ms)
        lines!(ax6, sort(df, x[i])[!,x[i]], mean_of_bins!(df,x[i],y[6],nb_bins), color = :red)
        hlines!(ax6, mean(df[:,y[6]]), color = :orange)
        scatter!(ax7, df[!,x[i]], df[!,y[7]], markersize = ms)
        lines!(ax7, sort(df, x[i])[!,x[i]], mean_of_bins!(df,x[i],y[7],nb_bins), color = :red)
        hlines!(ax7, mean(df[:,y[7]]), color = :orange)
        scatter!(ax8, df[!,x[i]], df[!,y[8]], markersize = ms)
        lines!(ax8, sort(df, x[i])[!,x[i]], mean_of_bins!(df,x[i],y[8],nb_bins), color = :red)
        hlines!(ax8, mean(df[:,y[8]]), color = :orange)
      
        ax7.xlabel = string(x[i], x_units[i])
        ax8.xlabel = string(x[i], x_units[i])
        ax1.ylabel = string(y[1], y_units[1])
        ax2.ylabel = string(y[2], y_units[2])
        ax3.ylabel = string(y[3], y_units[3])
        ax4.ylabel = string(y[4], y_units[4])
        ax5.ylabel = string(y[5], y_units[5])
        ax6.ylabel = string(y[6], y_units[6])
        ax7.ylabel = string(y[7], y_units[7])
        ax8.ylabel = string(y[8], y_units[8])

        # ylims!(ax1, (0, 320))
        # ylims!(ax2, (0, 43))
        # ylims!(ax3, (0, 21))
        # ylims!(ax4, (0, 5))
        # ylims!(ax5, (0, 40))
        # ylims!(ax6, (-30, 7.1))
        # ylims!(ax7, (0, 26))
        # ylims!(ax8, (0, 2.1))


        xlims!(ax1, x_ranges[i])
        xlims!(ax2, x_ranges[i])
        xlims!(ax3, x_ranges[i])
        xlims!(ax4, x_ranges[i])
        xlims!(ax5, x_ranges[i])
        xlims!(ax6, x_ranges[i])
        xlims!(ax7, x_ranges[i])
        xlims!(ax8, x_ranges[i])


        hidespines!(ax1, :t, :r)
        hidespines!(ax2, :t, :r)
        hidespines!(ax3, :t, :r)
        hidespines!(ax4, :t, :r)
        hidespines!(ax5, :t, :r)
        hidespines!(ax6, :t, :r)
        hidespines!(ax7, :t, :r)
        hidespines!(ax8, :t, :r)

        #supertitle = fig[0, :] = Label(fig, "L = $(round(mean(df.L))) km", textsize = 22)

        save(plotsdir("Open model","v2","Parameter scatter plots",land_dist,"om_" * land_dist * "_runs$(nb_rows)_$(x[i])_morequantities.png"),fig)
        #return fig
    end
end


function eight_scatter_plots_totflux(df::DataFrame, nb_bins, land_dist)
    x = ["w0", "τ", "α"] #["spwp", "ep", "eo", "ϵ", "α", "nZr", "a", "b", "wsat", "w0", "u_ms", "τ"]
    x_units = [" [mm]", " [1/day]", ""] #["", " [mm/day]", " [mm/day]", "", "", " [mm]", "", "", " [mm]", " [mm]", " [m/s]", " [1/day]"]
    #x_ranges = [(0, 70), (0, 1.8), (0, 1)]
    y = ["P1", "P2", "P3", "E2", "R", "Atot", "PRtot", "PR"] #quantset1
    #y = ["w1", "w2", "w3", "s", "A1", "A2", "A3", "Atot"] #quantset2
    y_units = [" [mm^2/day]", " [mm^2/day]", " [mm^2/day]", " [mm^2/day]", " [mm^2/day]", " [mm^2/day]", "", ""] #quantset1
    #y_units = [" [mm^2]", " [mm^2]", " [mm^2]", "[mm^2]", " [mm]", " [mm]", " [mm]", " [mm]"] #quantset2
    ms = 5.0
    nb_rows = nrow(df)

    for i = 1:length(x)
        fig = Figure(resolution = (1000, 1500))
        ax1 = Axis(fig[1, 1])
        ax2 = Axis(fig[1, 2])
        ax3 = Axis(fig[2, 1])
        ax4 = Axis(fig[2, 2])
        ax5 = Axis(fig[3, 1])
        ax6 = Axis(fig[3, 2])
        ax7 = Axis(fig[4, 1])
        ax8 = Axis(fig[4, 2])
        scatter!(ax1, df[!,x[i]], df[!,y[1]], markersize = ms)
        lines!(ax1, sort(df, x[i])[!,x[i]], mean_of_bins!(df,x[i],y[1],nb_bins), color = :red)
        hlines!(ax1, mean(df[:,y[1]]), color = :orange)
        scatter!(ax2, df[!,x[i]], df[!,y[2]], markersize = ms)
        lines!(ax2, sort(df, x[i])[!,x[i]], mean_of_bins!(df,x[i],y[2],nb_bins), color = :red)
        hlines!(ax2, mean(df[:,y[2]]), color = :orange)
        scatter!(ax3, df[!,x[i]], df[!,y[3]], markersize = ms)
        lines!(ax3, sort(df, x[i])[!,x[i]], mean_of_bins!(df,x[i],y[3],nb_bins), color = :red)
        hlines!(ax3, mean(df[:,y[3]]), color = :orange)
        scatter!(ax4, df[!,x[i]], df[!,y[4]], markersize = ms)
        lines!(ax4, sort(df, x[i])[!,x[i]], mean_of_bins!(df,x[i],y[4],nb_bins), color = :red)
        hlines!(ax4, mean(df[:,y[4]]), color = :orange)
        scatter!(ax5, df[!,x[i]], df[!,y[5]], markersize = ms)
        lines!(ax5, sort(df, x[i])[!,x[i]], mean_of_bins!(df,x[i],y[5],nb_bins), color = :red)
        hlines!(ax5, mean(df[:,y[5]]), color = :orange)
        scatter!(ax6, df[!,x[i]], df[!,y[6]], markersize = ms)
        lines!(ax6, sort(df, x[i])[!,x[i]], mean_of_bins!(df,x[i],y[6],nb_bins), color = :red)
        hlines!(ax6, mean(df[:,y[6]]), color = :orange)
        scatter!(ax7, df[!,x[i]], df[!,y[7]], markersize = ms)
        lines!(ax7, sort(df, x[i])[!,x[i]], mean_of_bins!(df,x[i],y[7],nb_bins), color = :red)
        hlines!(ax7, mean(df[:,y[7]]), color = :orange)
        scatter!(ax8, df[!,x[i]], df[!,y[8]], markersize = ms)
        lines!(ax8, sort(df, x[i])[!,x[i]], mean_of_bins!(df,x[i],y[8],nb_bins), color = :red)
        hlines!(ax8, mean(df[:,y[8]]), color = :orange)
      
        ax7.xlabel = string(x[i], x_units[i])
        ax8.xlabel = string(x[i], x_units[i])
        ax1.ylabel = string(y[1], y_units[1])
        ax2.ylabel = string(y[2], y_units[2])
        ax3.ylabel = string(y[3], y_units[3])
        ax4.ylabel = string(y[4], y_units[4])
        ax5.ylabel = string(y[5], y_units[5])
        ax6.ylabel = string(y[6], y_units[6])
        ax7.ylabel = string(y[7], y_units[7])
        ax8.ylabel = string(y[8], y_units[8])

        # ylims!(ax1, (0, 320))
        # ylims!(ax2, (0, 43))
        # ylims!(ax3, (0, 21))
        # ylims!(ax4, (0, 5))
        # ylims!(ax5, (0, 40))
        # ylims!(ax6, (-30, 7.1))
        # ylims!(ax7, (0, 26))
        # ylims!(ax8, (0, 2.1))


        # xlims!(ax1, x_ranges[i])
        # xlims!(ax2, x_ranges[i])
        # xlims!(ax3, x_ranges[i])
        # xlims!(ax4, x_ranges[i])
        # xlims!(ax5, x_ranges[i])
        # xlims!(ax6, x_ranges[i])
        # xlims!(ax7, x_ranges[i])
        # xlims!(ax8, x_ranges[i])


        hidespines!(ax1, :t, :r)
        hidespines!(ax2, :t, :r)
        hidespines!(ax3, :t, :r)
        hidespines!(ax4, :t, :r)
        hidespines!(ax5, :t, :r)
        hidespines!(ax6, :t, :r)
        hidespines!(ax7, :t, :r)
        hidespines!(ax8, :t, :r)

        #supertitle = fig[0, :] = Label(fig, "L = $(round(mean(df.L))) km", textsize = 22)

        save(plotsdir("Open model","v2","Parameter scatter plots",land_dist,"om_" * land_dist * "_runs$(nb_rows)_$(x[i])_quantset1_total.png"),fig)
        #return fig
    end
end


function om_param_vars_scatter_plots(df::DataFrame, system, nb_bins = 100)
    ms = 5.0
    domain_size = Int(round(df[1, "L"], digits=0))
    nb_runs = nrow(df)

    x = ["spwp", "ep", "eo", "ϵ", "r", "α", "nZr", "a", "b", "wsat", "w0", "u"]
    x_units = ["", " [mm/day]", " [mm/day]", "", "", "", " [mm]", "", "", " [mm]", " [mm]", " [m/s]"]
    y = ["s", "w1", "w2", "w3"] 
    y_units = ["", " [mm]", " [mm]", " [mm]"]
    
    for n = 1:length(x)
        
        fig = Figure(resolution = (800, 1200))
        ax1 = Axis(fig[1, :], title = "Domain size = $(domain_size) km", textsize = 22)
        ax2 = Axis(fig[2, :])
        ax3 = Axis(fig[3, :])
        ax4 = Axis(fig[4, :])
        
        scatter!(ax1, df[!,x[n]], df[!,y[1]], markersize = ms)
        lines!(ax1, sort(df, x[n])[!,x[n]], mean_of_bins!(df,x[n],y[1],nb_bins), color = :red)
        hlines!(ax1, mean(df[:,y[1]]), color = :orange)
        scatter!(ax2, df[!,x[n]], df[!,y[2]], markersize = ms)
        lines!(ax2, sort(df, x[n])[!,x[n]], mean_of_bins!(df,x[n],y[2],nb_bins), color = :red)
        hlines!(ax2, mean(df[:,y[2]]), color = :orange)
        scatter!(ax3, df[!,x[n]], df[!,y[3]], markersize = ms)
        lines!(ax3, sort(df, x[n])[!,x[n]], mean_of_bins!(df,x[n],y[3],nb_bins), color = :red)
        hlines!(ax3, mean(df[:,y[3]]), color = :orange)
        scatter!(ax4, df[!,x[n]], df[!,y[4]], markersize = ms)
        lines!(ax4, sort(df, x[n])[!,x[n]], mean_of_bins!(df,x[n],y[4],nb_bins), color = :red)
        hlines!(ax4, mean(df[:,y[4]]), color = :orange)
        
        ax4.xlabel = string(x[n], x_units[n])
        
        ax1.ylabel = string(y[1],y_units[1])
        ax2.ylabel = string(y[2],y_units[2])
        ax3.ylabel = string(y[3],y_units[3])
        ax4.ylabel = string(y[3],y_units[4])

        hidespines!(ax1, :t, :r)
        hidespines!(ax2, :t, :r)
        hidespines!(ax3, :t, :r)
        hidespines!(ax4, :t, :r)

        #supertitle = fig[0, :] = Label(fig, , textsize = 22)
        save(plotsdir("Open model", system, "State variables scatter plots", "om_runs$(nb_runs)_domain$(domain_size)_$(x[n]).png"),fig)
        #return fig
    end
end


function one_scatter_plot(df::DataFrame, 
    x::String, 
    y::String, 
    xlabel::String, 
    ylabel::String, 
    domain::String,
    c = :darkorange2, ms = 5.0)
    l = nrow(df)
    d = round(mean(df.L), digits = 0)
    fig = Figure()
    ax = Axis(fig[1,1], title = domain)
    scatter!(ax, df[!,x], df[!,y], color = c, markersize = ms)
    #hlines!(ax, 0.0, linestyle = :dash, color = :gray22)
    #vlines!(ax, 0.0, linestyle = :dash, color = :gray22)
    ax.xlabel = xlabel
    ax.ylabel = ylabel
    hidespines!(ax, :t, :r)
    #save(plotsdir("Open model","v2","Parameter scatter plots", "runs$(l)_domain$(d)_$(x)_vs_$(y).png"),fig)
    return fig
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


### Presentation plots ###

function chi_pdf(data::DataFrame)
    lfs = 22
    fig = Figure(resolution = (500, 400))
    ax = Axis(fig[1, 1], xlabel = L"Precipitation ratio $\chi$", ylabel = L"\mathrm{PDF}", ylabelsize = lfs, xlabelsize = lfs, xgridcolor = :white, ygridcolor = :white)
    #lines!(kde(data[!,"PR"]), color = :grey30, linewidth = 2.5, label = L"\mathrm{data}\, \, \mathrm{from}\, \, 50000\, \, \mathrm{runs}")
    hist!(ax, data[!,"PR"], bins = 500, normalization = :pdf, label = L"\mathrm{data}\, \, \mathrm{from}\, \, 100000\, \, \mathrm{runs}", color = :grey30)
    vlines!(ax, mean(dcm.PR), color = :red, linewidth = 1.5, label = L"\chi_\mathrm{mean} \, =\, 0.94")
    #vlines!(ax, 0.9, color = :dodgerblue4, linestyle = :dash, linewidth = 1.5, label = L"\mathrm{range}\, \, \mathrm{from}\, \, \mathrm{obs}")
    #vlines!(ax, 1.04, color = :dodgerblue4, linestyle = :dash, linewidth = 1.5)
    axislegend(ax, position = :lt, framevisible = false, labelsize = 20)
    hidespines!(ax, :t, :r)
    return fig
end

function chi_param(data::DataFrame, xquant::String, yquant::String = "PR")
    l = full_labels_dict()
    lfs = 22
    ms = 5
    c1 = (:grey35, 0.8)
    fig = Figure(resolution = (500, 400))
    ax = Axis(fig[1, 1], xlabel = l[xquant], ylabel = l[yquant], ylabelsize = lfs, xlabelsize = lfs, xgridcolor = :white, ygridcolor = :white)
    scatter!(ax, data[!,xquant], data[!,yquant], markersize = ms, color = c1, label = L"\mathrm{data}\, \, \mathrm{from}\, \, 50000\, \, \mathrm{runs}")
    lines!(ax, sort!(data, xquant)[!,xquant], movingaverage(data[!,yquant], 20000), color = :orange, label = L"\mathrm{mean}")
    hidespines!(ax, :t, :r)
    #axislegend(ax, position = :rb, framevisible = false, labelsize = 20)
    return fig
end


function cm_statevar_pdfs(data::DataFrame)
    l = full_labels_dict()
    lfs = 24
    tls = 18
    lw = 2.5
    f1 = Figure(resolution = (1200, 500))
    ax1 = Axis(f1[1, 1], xlabel = l["s"], ylabel = L"\mathrm{PDF}", ylabelsize = lfs, xlabelsize = lfs, xgridvisible = false, ygridvisible = false, yticklabelsize = tls, xticklabelsize = tls)
    lines!(kde(data[!,"s"]), color = :grey30, linewidth = lw)

    ax2 = Axis(f1[1,2], xlabel = L"\tilde{s} = \frac{s - s_\mathrm{pwp}}{s_\mathrm{sfc}-s_\mathrm{pwp}}", xlabelsize = lfs, xgridvisible = false, ygridvisible = false, yticklabelsize = tls, xticklabelsize = tls)
    data.sresc = (data.s .- data.spwp) ./ (data.sfc .- data.spwp)
    lines!(kde(data[!,"sresc"]), color = :grey30, linewidth = lw)
    vlines!(ax2, 0.0, linewidth = 2.0, linestyle = :dash, color = :grey40)
    vlines!(ax2, 1.0, linewidth = 2.0, linestyle = :dash, color = :grey40)

    ax3 = Axis(f1[1, 3], xlabel = L"Water vapor pass $w$", xlabelsize = lfs, xgridvisible = false, ygridvisible = false, yticklabelsize = tls, xticklabelsize = tls)
    lines!(kde(data[!,"wo"]), color = :grey30, linewidth = lw, label = L"\mathrm{ocean}")
    lines!(kde(data[!,"wl"]), color = :grey30, linestyle = :dot, linewidth = lw + 0.5, label = L"\mathrm{land}")

    axislegend(ax3, position = :lt, framevisible = false, labelsize = 24)
    hidespines!(ax1, :t, :r)
    hidespines!(ax2, :t, :r)
    hidespines!(ax3, :t, :r)
    #xlims!(ax3, (8,60))
    vlines!(ax3, 48, color = :red)
    colsize!(f1.layout, 1, Relative(7/24))
    colsize!(f1.layout, 2, Relative(5/24))
    colsize!(f1.layout, 3, Relative(3/6))

    return f1
end

function P_alpha(data::DataFrame)
    l = full_labels_dict()
    lfs = 22
    ms = 5
    c1 = (:grey35, 0.8)
    fig = Figure(resolution = (500, 400))
    ax1 = Axis(fig[1, 1], xlabel = l["α"], ylabel = L"\mathrm{Precipitation}\,\,P\,\,\mathrm{[mm/day]}", xlabelsize = lfs, ylabelsize = lfs, xgridcolor = :white, ygridcolor = :white)
    #ax2 = Axis(fig[2, 1], xlabel = l["α"], ylabel = l["Po"], xlabelsize = lfs, ylabelsize = lfs, xgridcolor = :white, ygridcolor = :white)
    #scatter!(ax1, data[!,"α"], data[!,"Pl"], markersize = ms, color = c1)
    lines!(ax1, sort!(data, "α")[!,"α"], movingaverage(data[!,"Pl"], 20000), color = :darkorange3, label = L"\mathrm{land}")
    #scatter!(ax2, data[!,"α"], data[!,"Po"], markersize = ms, color = c1)
    lines!(ax1, sort!(data, "α")[!,"α"], movingaverage(data[!,"Po"], 20000), color = :dodgerblue4, label = L"\mathrm{ocean}")
    hidespines!(ax1, :t, :r)
    #hidespines!(ax2, :t, :r)
    axislegend(ax1, position = :lb, framevisible = false, labelsize = 20)
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
        "spwp" => L"s_\mathrm{pwp}",
        "sfc"  => L"s_\mathrm{fc}",
        "ep"   => L"e_\mathrm{p}", 
        "eo"   => L"e_\mathrm{o}",
        "ϵ"    => L"\epsilon",
        "r"    => L"r",  
        "α"    => L"\alpha", 
        "nZr"  => L"nz_\mathrm{r}",   
        "a"    => L"a",  
        "b"    => L"b",
        "wsat" => L"w_\mathrm{sat}", 
        "u"    => L"u",
        "L"    => L"L",
        "Pl"   => L"P_\mathrm{l}",
        "Po"   => L"$P_\mathrm{o}$",
        "El"   => L"$E_\mathrm{l}$",
        "R"    => L"$R$",
        "PR"   => L"$\chi$",
        "Ptot" => L"$P_\mathrm{mean}$",
        "A"    => L"$A_\mathrm{l}$",
        "B"    => L"$-A_\mathrm{o}$", 
        "τ"    => L"\tau",
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