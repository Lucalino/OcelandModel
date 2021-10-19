using DrWatson
@quickactivate "Oceland Model"
include(srcdir("parametrisations.jl"))
include(srcdir("utils.jl"))
using CairoMakie
using Colors
using Random
using Distributions

fs= 20.0 #fontsize
lw= 2.0  #linewidth
ms= 5.0  #markersize

"""
    evap_tanh_plot(mode::String)

Takes two different modes of computation: 
"pt" to produce a plot of smooth evapotranspiration curves with various tuning parameter values and
"pwp" to produce a plot of smooth evapotranspiration curves with various spwp and sfc and constant pt = 10.

"""
function evap_tanh_plot(mode::String)

    fig = Figure()
    ax  = Axis(fig[1,1])
    ax.xlabel = "soil moisture saturation"
    ax.ylabel = "evaporation rate [mm/day]"
    ax.xlabelsize = fs
    ax.ylabelsize = fs
    hidespines!(ax, :t, :r)
    
    c = [:dodgerblue, :chocolate1, :green, :orangered2]
    s    = collect(0.0:0.01:1.0)
    Ep   = 4.38
    
    if mode == "pt"
        for n = 1:4
            spwp = 0.4
            sfc  = spwp + 0.3
            A    = Ep/2
            B    = 16 - 2 * n
            C    = (spwp + sfc)/2
            parameters = @dict A B C      
            
            if n == 1
                piecewise = lines!(ax, s, [land_evap(el, spwp, sfc, Ep) for el in s], linestyle = :dot, color = :gray22, linewidth = 2.5, label = "piecewise parametrisation")
            end

            new = lines!(ax, s, [evap_tanh(el, parameters) for el in s], linewidth = lw, color = c[n], label = "tuning parameter = $(B)")
            ax.title = "s_pwp = 0.4, s_fc = 0.7"
        
            save(plotsdir("Sketches", "Evap_tanh_p_tuning_varied.png"), fig)

        end
    
    elseif mode == "pwp"
        for n = 1:4
            spwp = 0.1 + n * 0.1
            sfc  = spwp + 0.3
            A    = Ep/2
            B    = 10
            C    = (spwp + sfc)/2
            parameters = @dict A B C      
      
            new = lines!(ax, s, [evap_tanh(el, parameters) for el in s], linewidth = lw, color = c[n], label = "s_pwp = $(round(spwp, digits=1)), s_fc = $(round(sfc, digits=1))")
            piecewise = lines!(ax, s, [land_evap(el, spwp, sfc, Ep) for el in s], linestyle = :dash, color = c[n], linewidth = lw)
            ax.title = "Tuning parameter = 10"
            axislegend(position = :lt, labelsize = 16, framevisible = false)
            #save(plotsdir("Sketches", "Evap_tanh_pwp_varied.png"), fig)
        end

    else
        println("Which parameter do you want to vary?")
    
    end

    return fig
end