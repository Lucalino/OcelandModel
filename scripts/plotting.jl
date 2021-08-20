using DrWatson
@quickactivate "Oceland Model"
using DataFrames
using CSV
using CairoMakie
using Statistics
using Colors
include(srcdir("parametrisations.jl"))
include(srcdir("utils.jl"))
include(srcdir("analysis.jl"))

fs= 14.0 #fontsize
lw= 2.0  #linewidth
ms= 5.0  #markersize

#read in raw solution data (parameters and solutions)
df_raw = CSV.read(datadir("sims", "om_eq_MonteCarlo_scan_100000_runs.csv"), DataFrame)

#add derived quantities to dataframe
df = derived_quantities!(df_raw)

nb_bins = 1000

function eight_scatter_plots(df::DataFrame, nb_bins = nb_bins)
    x = ["spwp", "sfc", "ϵ", "r", "w0", "Li", "u", "w_sat"]
    y = ["P1", "P2", "P3", "El", "infilt", "runoff", "PR"]
    y_units = [" [mm/day]", " [mm/day]", " [mm/day]", " mm/day", "", " [mm/day]", ""]
    x_units = ["", "", "", "", " [mm]", " [km]", " [m/s]", " [mm]"]

    for i = 1:length(y)
        fig = Figure(resolution = (1000, 1500))
        ax1 = Axis(fig[1, 1])
        ax2 = Axis(fig[1, 2])
        ax3 = Axis(fig[2, 1])
        ax4 = Axis(fig[2, 2])
        ax5 = Axis(fig[3, 1])
        ax6 = Axis(fig[3, 2])
        ax7 = Axis(fig[4, 1])
        ax8 = Axis(fig[4, 2])
        scatter!(ax1, df[!,x[1]], df[!,y[i]], markersize = ms)
        lines!(ax1, sort(df, x[1])[!,x[1]], mean_of_bins!(df,x[1],y[i],nb_bins), color = :red)
        hlines!(ax1, mean(df[:,y[i]]), color = :orange)
        scatter!(ax2, df[!,x[2]], df[!,y[i]], markersize = ms)
        lines!(ax2, sort(df, x[2])[!,x[2]], mean_of_bins!(df,x[2],y[i],nb_bins), color = :red)
        hlines!(ax2, mean(df[:,y[i]]), color = :orange)
        scatter!(ax3, df[!,x[3]], df[!,y[i]], markersize = ms)
        lines!(ax3, sort(df, x[3])[!,x[3]], mean_of_bins!(df,x[3],y[i],nb_bins), color = :red)
        hlines!(ax3, mean(df[:,y[i]]), color = :orange)
        scatter!(ax4, df[!,x[4]], df[!,y[i]], markersize = ms)
        lines!(ax4, sort(df, x[4])[!,x[4]], mean_of_bins!(df,x[4],y[i],nb_bins), color = :red)
        hlines!(ax4, mean(df[:,y[i]]), color = :orange)
        scatter!(ax5, df[!,x[5]], df[!,y[i]], markersize = ms)
        lines!(ax5, sort(df, x[5])[!,x[5]], mean_of_bins!(df,x[5],y[i],nb_bins), color = :red)
        hlines!(ax5, mean(df[:,y[i]]), color = :orange)
        scatter!(ax6, df[!,x[6]], df[!,y[i]], markersize = ms)
        lines!(ax6, sort(df, x[6])[!,x[6]], mean_of_bins!(df,x[6],y[i],nb_bins), color = :red)
        hlines!(ax6, mean(df[:,y[i]]), color = :orange)
        scatter!(ax7, df[!,x[7]], df[!,y[i]], markersize = ms)
        lines!(ax7, sort(df, x[7])[!,x[7]], mean_of_bins!(df,x[7],y[i],nb_bins), color = :red)
        hlines!(ax7, mean(df[:,y[i]]), color = :orange)
        scatter!(ax8, df[!,x[8]], df[!,y[i]], markersize = ms)
        lines!(ax8, sort(df, x[8])[!,x[8]], mean_of_bins!(df,x[8],y[i],nb_bins), color = :red)
        hlines!(ax8, mean(df[:,y[i]]), color = :orange)
        ax1.xlabel = string(x[1], x_units[1])
        ax1.ylabel = string(y[i],y_units[i])
        ax2.xlabel = string(x[2], x_units[2])
        ax3.xlabel = string(x[3], x_units[3])
        ax3.ylabel = string(y[i],y_units[i])
        ax4.xlabel = string(x[4], x_units[4])
        ax5.xlabel = string(x[5], x_units[5])
        ax5.ylabel = string(y[i],y_units[i])
        ax6.xlabel = string(x[6], x_units[6])
        ax7.xlabel = string(x[7], x_units[7])
        ax7.ylabel = string(y[i],y_units[i])
        ax8.xlabel = string(x[8], x_units[8])

        hidespines!(ax1, :t, :r)
        hidespines!(ax2, :t, :r)
        hidespines!(ax3, :t, :r)
        hidespines!(ax4, :t, :r)
        hidespines!(ax5, :t, :r)
        hidespines!(ax6, :t, :r)
        hidespines!(ax7, :t, :r)
        hidespines!(ax8, :t, :r)

        supertitle = fig[0, :] = Label(fig, "w0 < w_sat", textsize = 22)

        save(plotsdir("Open model","Parameter scatter plots","100000_$(y[i])_params_scatter.png"),fig)
        #return fig
    end
end

function one_scatter_plot(df::DataFrame, 
                           x::String, 
                           y::String, 
                      xlabel::String, 
                      ylabel::String, 
                      c = :darkorange2, ms = 5.0)
    fig = Figure()
    ax  = Axis(fig[1,1], title = "w0 < w_sat")
    scatter!(ax, df[!,x], df[!,y], color = c, markersize = ms)
    hlines!(ax, 0.0, linestyle = :dash, color = :gray22)
    vlines!(ax, 0.0, linestyle = :dash, color = :gray22)
    ax.xlabel = xlabel
    ax.ylabel = ylabel
    hidespines!(ax, :t, :r)
    #supertitle = fig[0, :] = Label(fig, "w0 < w_sat", textsize = 22)
    #return fig
    save(plotsdir("Open model","Water vapour changes","$(x)_vs_$(y).png"),fig)
    return fig
end

function colormap_scatterplot(df::DataFrame,
                       x::String, 
                       y::String, 
                       z::String,
                       xlabel::String, 
                       ylabel::String,
                       zlabel::String)
    fig = Figure()
    ax  = Axis(fig[1,1])
    pnts = scatter!(ax,df[!,x], df[!,y], color = df[!,z], colormap = :darkrainbow, markersize = 3)
    cbar = Colorbar(fig, pnts, height = Relative(0.75), label = zlabel)
    fig[1, 2] = cbar
    ax.xlabel = xlabel
    ax.ylabel = ylabel
    ax.title  = "N = 50 000"
    save(plotsdir("Open model", "Parameter scatter plots", "$(x)_vs_$(y)_vs_$(z)_N=50000.png"), fig)
    return fig
end

function state_variable_plot(df::DataFrame)
    c = :dimgrey
    fig = Figure(resolution = (1000, 1200))
    ax1, ax2, ax3, ax4, ax5, ax6  = Axis(fig[1,1]), Axis(fig[1,2]), Axis(fig[2,1]), Axis(fig[2,2]), Axis(fig[3,1]), Axis(fig[3,2])
    scatter!(ax1, df.w1, df.w2, markersize = ms, color = c)
    scatter!(ax2, df.w1, df.s, markersize = ms, color = c)
    scatter!(ax3, df.w1, df.w3, markersize = ms, color = c)
    scatter!(ax4, df.w2, df.s, markersize = ms, color = c)
    scatter!(ax5, df.w2, df.w3, markersize = ms, color = c)
    scatter!(ax6, df.s, df.w3, markersize = ms, color = c)
    ax1.xlabel = "w1 [mm]"
    ax1.ylabel = "w2 [mm]"
    ax2.xlabel = "w1 [mm]"
    ax2.ylabel = "s"
    ax3.xlabel = "w1 [mm]"
    ax3.ylabel = "w3 [mm]"
    ax4.xlabel = "w2 [mm]"
    ax4.ylabel = "s"
    ax5.xlabel = "w2 [mm]"
    ax5.ylabel = "w3 [mm]"
    ax6.xlabel = "s"
    ax6.ylabel = "w3 [mm]"
    hidespines!(ax1, :t, :r)
    hidespines!(ax2, :t, :r)
    hidespines!(ax3, :t, :r)
    hidespines!(ax4, :t, :r)
    hidespines!(ax5, :t, :r)
    hidespines!(ax6, :t, :r)
    fig[0, :] = Label(fig, "State variables, w0 < w_sat", textsize = 22)
    save(plotsdir("Open model", "State variable scatter plots", "state_variables_N=100000.png"), fig) 
end