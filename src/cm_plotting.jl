
function cm_t_evolution_plot(s, wl, wo, t, 
                            fs = 14.0,
                            c1 = "dodgerblue", 
                            c2 = "darkblue",
                            c3 = "darkgreen")
    fig = Figure()
    ax1 = Axis(fig[1,1])
    ax2 = Axis(fig[2,1])
    ax2.xlabel = "time [s]"
    ax1.ylabel = "water vapour pass [mm]"
    ax2.ylabel = "relative soil moisture"
    ax2.ylabelsize = fs
    ax2.xlabelsize = fs
    ax2.ylabelsize = fs
    hidespines!(ax1, :t, :r)
    hidespines!(ax2, :t, :r)
    lines!(ax1, t, wl, label="w_l", color=c1)
    lines!(ax1, t, wo, label="w_o", color=c2)
    lines!(ax2, t, s, color=c3)
    axislegend(ax1, framevisible = false)
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

    x = ["spwp", "sfc", "ϵ", "u", "wsat", "a", "b", "ep", "nZr", "α", "eo"]
    x_units = ["", "", "", " [m/s]", " [mm]", "", "", " [mm/day]", " [mm]", "", " [mm/day]"]
    #x = ["u"]
    #x_units = [" [m/s]"]
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
        
        scatter!(ax1, df[!,x[n]], df[!,y[1]], markersize = ms)
        lines!(ax1, sort(df, x[n])[!,x[n]], cm_mean_of_bins!(df,x[n],y[1],nb_bins), color = :red)
        #lines!(ax1, sort(df, x[n])[!,x[n]], cm_mean_of_bins!(df,x[n],y[2],nb_bins) ./ cm_mean_of_bins!(df, x[n], y[3], nb_bins), color = :red)
        hlines!(ax1, mean(df[:,y[1]]), color = :orange)
        scatter!(ax2, df[!,x[n]], df[!,y[2]], markersize = ms)
        lines!(ax2, sort(df, x[n])[!,x[n]], cm_mean_of_bins!(df,x[n],y[2],nb_bins), color = :red)
        hlines!(ax2, mean(df[:,y[2]]), color = :orange)
        scatter!(ax3, df[!,x[n]], df[!,y[3]], markersize = ms)
        lines!(ax3, sort(df, x[n])[!,x[n]], cm_mean_of_bins!(df,x[n],y[3],nb_bins), color = :red)
        hlines!(ax3, mean(df[:,y[3]]), color = :orange)
        scatter!(ax4, df[!,x[n]], df[!,y[4]], markersize = ms)
        lines!(ax4, sort(df, x[n])[!,x[n]], cm_mean_of_bins!(df,x[n],y[4],nb_bins), color = :red)
        hlines!(ax4, mean(df[:,y[4]]), color = :orange)
        scatter!(ax5, df[!,x[n]], df[!,y[5]], markersize = ms)
        lines!(ax5, sort(df, x[n])[!,x[n]], cm_mean_of_bins!(df,x[n],y[5],nb_bins), color = :red)
        hlines!(ax5, mean(df[:,y[5]]), color = :orange)
        scatter!(ax6, df[!,x[n]], df[!,y[6]], markersize = ms)
        lines!(ax6, sort(df, x[n])[!,x[n]], cm_mean_of_bins!(df,x[n],y[6],nb_bins), color = :red)
        hlines!(ax6, mean(df[:,y[6]]), color = :orange)
        
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

        supertitle = fig[0, :] = Label(fig, "Domain size = $(round(df[1, "L"], digits = 0)) km", textsize = 22)

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
    #scatter!(ax, df2[!,x], df2[!,y], color = (:chartreuse4, 0.6), markersize = 5.0)
    ax.xlabel = xlabel
    ax.ylabel = ylabel
    #save(plotsdir("Closed model", "Parameter scatter plots/", system * "/cm_10000runs_domain" * domain * "_" * x * "_" * y * ".png"),fig)
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
