function eight_scatter_plots(df::DataFrame, nb_bins)
    x = ["spwp", "ep", "eo", "ϵ", "r", "α", "nZr", "a", "b", "wsat", "w0", "u"]
    x_units = ["", " [mm/day]", " [mm/day]", "", "", "", " [mm]", "", "", " [mm]", " [mm]", " [m/s]"]
    y = ["P1", "P2", "P3", "El", "R", "Φ", "Ptot", "PR"]
    y_units = [" [mm/day]", " [mm/day]", " [mm/day]", " [mm/day]", " [mm/day]", "", " [mm/day]", ""]
    

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

        hidespines!(ax1, :t, :r)
        hidespines!(ax2, :t, :r)
        hidespines!(ax3, :t, :r)
        hidespines!(ax4, :t, :r)
        hidespines!(ax5, :t, :r)
        hidespines!(ax6, :t, :r)
        hidespines!(ax7, :t, :r)
        hidespines!(ax8, :t, :r)

        supertitle = fig[0, :] = Label(fig, "L = $(round(mean(df.L))) km", textsize = 22)

        save(plotsdir("Open model","v2","Parameter scatter plots","runs9950_$(x[i])_params_scatter_domain10000.png"),fig)
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
