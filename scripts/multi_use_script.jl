using DrWatson
@quickactivate "Oceland Model"

using CairoMakie
using GLMakie
#using Colors
#using Random
using DataFrames
#using Distributions
using DifferentialEquations
using DynamicalSystems
using InteractiveDynamics
include(srcdir("parametrisations.jl"))
include(srcdir("utils.jl"))
include(srcdir("model_versions.jl"))
include(srcdir("cm_analysis.jl"))


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



function diffeq_issue()
    #x0 = [0.3, 50.0, 40.0]
    x0 = [0.3, 0.5, 0.4]
    #p = [2.0, 70.0, 4.0e8]
    p = Dict(
        :a => 2.0,
        :b => 70.0,
        :c => 4.0e8,
    )

    function my_system(x, p, t)
        @unpack a, b, c = p
        dx = exp(x[2] - 0.5) * (1.0 - x[1]) - a * tanh(10.0 * (x[1] - 0.5)) - a
        dy = a * tanh(10.0 * (x[1] - 0.5)) + a - exp(x[2] - 0.5) + (x[3] - x[2]) * 2.0
        dz = 3.0 - exp(x[3] - 0.5) - (x[3] - x[2]) * 0.5
        # dx = (15.0 * (exp(x[2]/b - 0.6)) * (1.0 - x[1]^2) - a * tanh(10.0 * (x[1] - 0.5)) - a)/100.0
        # dy = a * tanh(10.0 * (x[1] - 0.5)) + a - exp(15.0 * (x[2]/b - 0.6)) + (x[3] - x[2]) * 2.0
        # dz = 3.0 - exp(15.0 * (x[3]/b - 0.6)) - (x[3] - x[2]) * 0.5
        return SVector(dx, dy, dz)
    end

    #p_fixed = cm_fixed_params()
    ds = ContinuousDynamicalSystem(my_system, x0, p)
    u0s =  [ds.u0]
    diffeq = (alg = Rodas5(), adaptive = false, dt = 0.001, reltol = 1e-8, abstol = 1e-8)
    tr = trajectory(ds, 100; diffeq...)
    ps = Dict(
        :a => 1.9:0.01:2.1
        #:spwp  => 0.2:0.01:0.54,   
        #:eo => 2.5:0.1:3.5,
    )

    fig, obs = interactive_evolution_timeseries(
        ds, u0s, ps; tail = 100000, diffeq, idxs = (1, 2, 3), 
        lims =((0.0,1.0), (0.0, 2.0), (0.0, 2.0))
    )
    #prob = ODEProblem(my_system, x0, tspan, p)
    #sol = solve(prob, alg = Rodas5(), reltol=1e-9, abstol=1e-9, saveat=1.0)
    return fig
end

td = cm_rand_params()