using DrWatson
@quickactivate "Oceland Model"
using CairoMakie
using DataFrames
using Colors

pw1000 = CSV.read(datadir("sims", "Closed model pmscan/cm_piecewise_eq_MC_fixedpoints_10000_runs_domain1000_all_quantities.csv"), DataFrame)
pw10000 = CSV.read(datadir("sims", "Closed model pmscan/cm_piecewise_eq_MC_fixedpoints_10000_runs_domain10000_all_quantities.csv"), DataFrame)

fig1 = Figure()
ax1  = Axis(fig1[1,1], xlabel = "u [m/s]", ylabel = "PR", title = "Full dataset, L = 1000 km")
scatter!(ax1, pw1000.u, pw1000.PR)

fig2 = Figure()
ax2  = Axis(fig2[1,1], xlabel = "u [m/s]", ylabel = "PR", title = "Full dataset, L = 10000 km")
scatter!(ax2, pw10000.u, pw10000.PR)

pw1000la = pw1000[(pw1000.α .> 0.1) .& (pw1000.α .< 0.15), :]
pw10000la = pw10000[(pw10000.α .> 0.1) .& (pw10000.α .< 0.15), :]
pw1000ma = pw1000[(pw1000.α .> 0.4) .& (pw1000.α .< 0.45), :]
pw10000ma = pw10000[(pw10000.α .> 0.4) .& (pw10000.α .< 0.45), :]
pw1000ha = pw1000[(pw1000.α .> 0.7) .& (pw1000.α .< 0.75), :]
pw10000ha = pw10000[(pw10000.α .> 0.7) .& (pw10000.α .< 0.75), :]

fig3 = Figure()
ax3  = Axis(fig3[1,1], xlabel = "u [m/s]", ylabel = "PR", title = "land fraction 0.1 < α < 0.15")
scatter!(ax3, pw1000la.u, pw1000la.PR, color = "dodgerblue", label = "L=1000 km")
scatter!(ax3, pw10000la.u, pw10000la.PR, color = "dodgerblue4", label = "L=10000 km")
axislegend(ax3, position = :rb, framevisible = false)

fig4 = Figure()
ax4  = Axis(fig4[1,1], xlabel = "u [m/s]", ylabel = "PR", title = "land fraction 0.4 < α < 0.45")
scatter!(ax4, pw1000ma.u, pw1000ma.PR, color = "dodgerblue", label = "L=1000 km")
scatter!(ax4, pw10000ma.u, pw10000ma.PR, color = "dodgerblue4", label = "L=10000 km")
axislegend(ax4, position = :rb, framevisible = false)

fig5 = Figure()
ax5  = Axis(fig5[1,1], xlabel = "u [m/s]", ylabel = "PR", title = "land fraction 0.7 < α < 0.75")
scatter!(ax5, pw1000ha.u, pw1000ha.PR, color = "dodgerblue", label = "L=1000 km")
scatter!(ax5, pw10000ha.u, pw10000ha.PR, color = "dodgerblue4", label = "L=10000 km")
axislegend(ax5, position = :rb, framevisible = false)

fig6 = Figure()
ax6  = Axis(fig6[1,1], xlabel = "u [m/s]", ylabel = "P_l [mm/day]", title = "land fraction 0.1 < α < 0.15")
scatter!(ax6, pw1000la.u, pw1000la.Pl, color = "dodgerblue", label = "L=1000 km")
scatter!(ax6, pw10000la.u, pw10000la.Pl, color = "dodgerblue4", label = "L=10000 km")
axislegend(ax6, position = :rb, framevisible = false)

fig7 = Figure()
ax7  = Axis(fig7[1,1], xlabel = "u [m/s]", ylabel = "P_l [mm/day]", title = "land fraction 0.4 < α < 0.45")
scatter!(ax7, pw1000ma.u, pw1000ma.Pl, color = "dodgerblue", label = "L=1000 km")
scatter!(ax7, pw10000ma.u, pw10000ma.Pl, color = "dodgerblue4", label = "L=10000 km")
axislegend(ax7, position = :rb, framevisible = false)

fig8 = Figure()
ax8  = Axis(fig8[1,1], xlabel = "u [m/s]", ylabel = "P_l [mm/day]", title = "land fraction 0.7 < α < 0.75")
scatter!(ax8, pw1000ha.u, pw1000ha.Pl, color = "dodgerblue", label = "L=1000 km")
scatter!(ax8, pw10000ha.u, pw10000ha.Pl, color = "dodgerblue4", label = "L=10000 km")
axislegend(ax8, position = :rb, framevisible = false)