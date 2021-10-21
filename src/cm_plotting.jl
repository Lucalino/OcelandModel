
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