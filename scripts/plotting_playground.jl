using DrWatson
@quickactivate "Oceland Model"

using CairoMakie
using ColorSchemes
using DataFrames
using LaTeXStrings
using DynamicalSystems
include(srcdir("parametrizations.jl"))
include(srcdir("figure_labels.jl"))
include(srcdir("create_model_output.jl"))
include(srcdir("utils.jl"))

lfs = 20
lw = 2.5


function time_series(data::DataFrame, y::String, x::String = "t")
    l = full_labels_dict()
    lfs = 20
    lw = 3.0
    fig = Figure(resolution = (500,400))
    ax = Axis(fig[1,1], xlabel = l[x], ylabel = l[y], ylabelsize = lfs, xlabelsize = lfs, xgridcolor = :white, ygridcolor = :white)
    lines!(ax, data[!, x], data[!, y], linewidth = lw)
    #xlims!(ax, -10, 10)
    hidespines!(ax, :t, :r)
    return fig
end

function time_series_exp(data::DataFrame, y::String, x::String = "t")
    l = full_labels_dict()
    lfs = 20
    lw = 3.0
    fig = Figure(resolution = (500,400))
    ax = Axis(fig[1,1], xlabel = l[x], ylabel = l[y], ylabelsize = lfs, xlabelsize = lfs, xgridcolor = :white, ygridcolor = :white)
    lines!(ax, data[!, x], data[!, y], linewidth = lw)
    # vlines!(ax, 2, color=:black, linestyle = :dot)
    # vlines!(ax, 60, color=:black, linestyle = :dot)
    # vlines!(ax, 130, color=:black, linestyle = :dot)
    # vlines!(ax, 240, color=:black, linestyle = :dot)
    # vlines!(ax, 450, color=:black, linestyle = :dot)
    #ylims!(ax, (0, 10^(-10)))
    xlims!(ax, 900, 910)
    hidespines!(ax, :t, :r)
    return fig
end

# function multiple_time_series(data_var::Vector{Vector{Float64}}, data_time::Vector{Vector{Float64}}, lab::Vector{Float64})
#     l = full_labels_dict()
#     lfs = 20
#     lw = 3.0
#     fig = Figure(resolution = (600,400))
#     ax = Axis(fig[1,1], xlabel = l["t"], ylabel = l["PRmean"], ylabelsize = lfs, xlabelsize = lfs, xgridcolor = :white, ygridcolor = :white)
#     for n = 1:length(data_var)
#         lines!(ax, data_time[n][1:end-100], data_var[n][1:end-100], linewidth = lw, label = "α = $(lab[n])")
#     end
#     lines!(ax, data_time[1][1:end-100], [1.0 for el in data_time[1][1:end-100]], linewidth = lw, color = :black)
#     hidespines!(ax, :t, :r)
#     axislegend(ax, framevisible = false, labelsize = 16)
#     #xlims!(0,10)
#     #ylims!(0,5)
#     return fig
# end

function multiple_time_series(data::DataFrame, y1::String, y2::String, y3::String, x::String = "t")
    l = full_labels_dict()
    lbn = labels_norm_dict()
    lfs = 20
    lw = 3.0
    fig = Figure(resolution = (500,400))
    ax = Axis(fig[1,1], xticks = [0.0, 0.25, 0.5, 0.75, 1.0], xlabel = l[x], ylabel = L"\mathrm{Fluxes\, \, \, [mm/day]}", ylabelsize = lfs, xlabelsize = lfs, xgridcolor = :white, ygridcolor = :white)
    lines!(ax, data[!,x] .- 800, data[!,y1], linewidth = lw, label = lbn[y1])
    lines!(ax, data[!,x] .- 800, data[!,y2], linewidth = lw, label = lbn[y2])
    lines!(ax, data[!,x] .- 800, data[!,y3], linewidth = lw, label = lbn[y3])
    # lines!(ax, data[!,x], data[!,y1], linewidth = lw, label = L"s")
    # lines!(ax, data[!,x], data[!,y2]./45, linewidth = lw, label = L"w_\ell /45.0")
    # lines!(ax, data[!,x], data[!,y3]./45, linewidth = lw, label = L"w_\mathrm{o} /45.0")
    # vlines!(ax, 2, color=:black, linestyle = :dot)
    # vlines!(ax, 20, color=:black, linestyle = :dot)
    # vlines!(ax, 60, color=:black, linestyle = :dot)
    # vlines!(ax, 115, color=:black, linestyle = :dot)
    # vlines!(ax, 220, color=:black, linestyle = :dot)
    hidespines!(ax, :t, :r)
    axislegend(ax, framevisible = false, labelsize = 16, position = :rc)
    #xlims!(-10,400)
    #ylims!(0,5)
    return fig
end

function PR_param_plot(param::String, data, p)
    l = full_labels_dict()
    lfs = 20
    lw = 3.0

    x = data[:,1]
    PR= [precip(elm,p) for elm in data[:,3]] ./ [precip(elm,p) for elm in data[:,4]]
    
    fig = Figure(resolution = (500,400))
    ax = Axis(fig[1,1], xlabel = l[param], ylabel = l["PRmean"], ylabelsize = lfs, xlabelsize = lfs, xgridcolor = :white, ygridcolor = :white)
    lines!(ax, x, PR, linewidth = lw)
    hlines!(ax, 1.0, linestyle = :dot, color = :black, linewidth = 2.0)
    hidespines!(ax, :t, :r)
    return fig
end

function four_time_series(data::DataFrame, y1::String, y2::String, y3::String, y4::String, y5::String)
    l = labels_dict()
    lfs = 20
    lw = 2.5

    f1 = Figure(resolution = (1200, 400))
    ax1 = Axis(f1[1, 1], xlabel = l["t"], ylabel = L"$w$ [mm]", ylabelsize = lfs, xlabelsize = lfs, xgridvisible = false, ygridvisible = false)
    ax2 = Axis(f1[1,2], xlabel = l["t"], xlabelsize = lfs,  ylabel = l[y3], ylabelsize = lfs,xgridvisible = false, ygridvisible = false)
    ax3 = Axis(f1[1, 3], xlabel = l["t"], xlabelsize = lfs,  ylabel = l[y4], ylabelsize = lfs,xgridvisible = false, ygridvisible = false)
    ax4 = Axis(f1[1,4], xlabel = l["t"], xlabelsize = lfs,  ylabel = l[y5], ylabelsize = lfs,xgridvisible = false, ygridvisible = false)

    lines!(ax1, data[!,"t"], data[!,y1], color = :ForestGreen, linewidth = lw, label = L"w_\ell")  
    lines!(ax1, data[!,"t"], data[!,y2], linewidth = lw, label = L"w_\mathrm{o}")  
    lines!(ax2, data[!,"t"], data[!,y3], color = :ForestGreen, linewidth = lw)
    lines!(ax3, data[!,"t"], data[!,y4], linewidth = lw)
    lines!(ax4, data[!,"t"], data[!,y5], color = :black, linewidth = lw)

    hidespines!(ax1, :t, :r)
    hidespines!(ax2, :t, :r)
    hidespines!(ax3, :t, :r)
    hidespines!(ax4, :t, :r)
    axislegend(ax1, framevisible = false, labelsize = 16, position = :rb)
    #xlims!(ax3, (0,60))
    #xlims!(ax1, (-0.02, 0.82))
    ylims!(ax1, (0,70))
    #ylims!(ax2, (0, 150))
    ylims!(ax3, (2.5, 3.0))
    #colsize!(f1.layout, 1, Relative(5/24))
    #colsize!(f1.layout, 2, Relative(5/24))
    #colsize!(f1.layout, 3, Relative(1/3))
    #colsize!(f1.layout, 4, Relative(1/4))
    return f1
end


function wind_profile_xt(lfs = lfs, lw = lw)
    t = [0.0, 0.1, 0.25, 0.5, 0.6, 0.8]
    L = 100.0 #km
    boxlen = 1.0
    boxnb = Int(L / boxlen)
    α = 0.3 
    boxnb_land = Int(round(α * boxnb))
    boxnb_ocean = boxnb - boxnb_land
    lsmask = vcat(zeros(boxnb_ocean+1), ones(boxnb_land))
    x = collect(0.0:boxlen:L)

    fig = Figure(resolution = (600,400))
    ax = Axis(fig[1,1], xlabel = L"Distance $x$ [km]", ylabel = L"Wind speed $u$ [m/s]", ylabelsize = lfs, xlabelsize = lfs, xgridcolor = :white, ygridcolor = :white)
    for i in t
        uo = cos.(pi .* x ./ ((1-α)*L)) .* cos(2*pi*i)
        ul = -cos.(pi .* (x .- (1-α)*L) ./ (α * L)) .* cos(2*pi*i) 
        u = vcat(uo[lsmask .== 0.0], ul[lsmask .> 0.0])
        lines!(ax, x, u, label = "$(i)")
    end
    vlines!(ax, (1-α) * L, linestyle = :dash, color = :black, linewidth = 2.0)
    #vlines!(ax, 0, linestyle = :dash, color = :black, linewidth = 2.0)
    #vlines!(ax, L, linestyle = :dash, color = :black, linewidth = 2.0)
    text!((1-α)*L/3, -1.0, text = "ocean")
    text!((1-α)*L+α*L/3, -1.0, text = "land")
    hidespines!(ax, :t, :r)
    fig[1,2] = Legend(fig, ax, "Time t [days]", framevisible = false)
    return fig
end


function animation(s, p, lfs = lfs, lw = lw)
    @unpack x, α, L = p
    time = Observable(1)
    vapor = @lift(s[:, $time])
    wind = @lift((wind_DC_xt.(x, s.t[$time], Ref(p)) .* mm2m(1.0)./day2s(1.0)))
    fig = Figure()
    ax1 = Axis(fig[1,1], xlabel = "x [km]", ylabel = "Water vapor path [mm]", ylabelsize = lfs, xlabelsize = lfs, xgridcolor = :white, ygridcolor = :white, title = @lift("t = $(round(s.t[$time], digits = 1)) days"))
    ax2 = Axis(fig[2,1], xlabel = "x [km]", ylabel = "Wind [m/s]", ylabelsize = lfs, xlabelsize = lfs, xgridcolor = :white, ygridcolor = :white)
    lines!(ax1, mm2km.(x), vapor)
    vlines!(ax1, (1-α) * mm2km(L), linestyle = :dash, linewidth = 2.0)
    lines!(ax2, mm2km.(x), wind)
    vlines!(ax2, (1-α) * mm2km(L), linestyle = :dash, linewidth = 2.0)
    ylims!(ax1, (0, 80))
    hidespines!(ax1, :t, :r)
    hidespines!(ax2, :t, :r)
    # record(fig, "/Users/lucaschmidt/Documents/Oceland Model/plots/Project 2/animations/vapor_animation_pset1_10days_x0hom_AandC.mp4", 1:50:length(s.t); framerate = 100) do t
    #     time[] = t
    record(fig, "/Users/lucaschmidt/Documents/Oceland Model/plots/Project 2/animations/test.mp4", 1:50:length(s.t); framerate = 100) do t
        time[] = t
    end
end

function A_C_animation(s, p, lfs = lfs, lw = lw)
    A, C, AnC = AandC_matrices(s, p)
    @unpack x, α, L = p
    time = Observable(1)
    advection = @lift(A[$time])
    convergence = @lift(C[$time])
    adv_and_conv = @lift(AnC[$time])
    w = @lift(s.u[$time])
    wind = @lift((wind_DC_xt.(x, s.t[$time], Ref(p)) .* mm2m(1.0)./day2s(1.0)))
    dwind = @lift((dwind_dx_DC_xt.(x, s.t[$time], Ref(p)) .* mm2m(1.0)./(day2s(1.0)^2)))
    fig = Figure(resolution = (800, 1200))
    ax1 = Axis(fig[1,1], xlabel = "x [km]", ylabel = "Advection [mm/day]", ylabelsize = lfs, xlabelsize = lfs, xgridcolor = :white, ygridcolor = :white, title = @lift("t = $(round(s.t[$time], digits = 1)) days"))
    ax2 = Axis(fig[2,1], xlabel = "x [km]", ylabel = "Convergence [mm/day]", ylabelsize = lfs, xlabelsize = lfs, xgridcolor = :white, ygridcolor = :white)
    #ax3 = Axis(fig[3,1], xlabel = "x [km]", ylabel = "Adv. + Conv. [mm/day]", ylabelsize = lfs, xlabelsize = lfs, xgridcolor = :white, ygridcolor = :white)
    ax3 = Axis(fig[3,1], xlabel = "x [km]", ylabel = "Water vapor path [mm]", ylabelsize = lfs, xlabelsize = lfs, xgridcolor = :white, ygridcolor = :white)
    #ax3 = Axis(fig[3,1], xlabel = "x [km]", ylabel = "Change of wind [mm/day^2]", ylabelsize = lfs, xlabelsize = lfs, xgridcolor = :white, ygridcolor = :white)
    ax4 = Axis(fig[4,1], xlabel = "x [km]", ylabel = "Wind [m/s]", ylabelsize = lfs, xlabelsize = lfs, xgridcolor = :white, ygridcolor = :white)
    lines!(ax1, mm2km.(x), advection)
    vlines!(ax1, (1-α) * mm2km(L), linestyle = :dash, linewidth = 2.0)
    lines!(ax2, mm2km.(x), convergence)
    vlines!(ax2, (1-α) * mm2km(L), linestyle = :dash, linewidth = 2.0)
    #lines!(ax3, mm2km.(x), adv_and_conv)
    lines!(ax3, mm2km.(x), w)
    #lines!(ax3, mm2km.(x), dwind)
    vlines!(ax3, (1-α) * mm2km(L), linestyle = :dash, linewidth = 2.0)
    lines!(ax4, mm2km.(x), wind)
    vlines!(ax4, (1-α) * mm2km(L), linestyle = :dash, linewidth = 2.0)
    ylims!(ax1, (-20, 20))
    hidespines!(ax1, :t, :r)
    hidespines!(ax2, :t, :r)
    hidespines!(ax3, :t, :r)
    hidespines!(ax4, :t, :r)
    # record(fig, "/Users/lucaschmidt/Documents/Oceland Model/plots/Project 2/animations/vapor_animation_pset1_10days_x0hom_AandC.mp4", 1:50:length(s.t); framerate = 100) do t
    #     time[] = t
    record(fig, "/Users/lucaschmidt/Documents/Oceland Model/plots/Project 2/animations/Removing the Zigzag/AnC_α03_dt=00001_dx=10km_w.mp4", 1:50:length(s.t); framerate = 1000) do t
        time[] = t
    end
end



function AandC_matrices(s, p)
    nb_rows = size(s)[2] #number of time steps
    nb_cols = size(s)[1] #number of boxes in x direction
    A   = Array{Array{Float64}}(undef, nb_rows)
    C   = Array{Array{Float64}}(undef, nb_rows)
    AnC = Array{Array{Float64}}(undef, nb_rows)
    dwdx = zeros(nb_cols)
    #dwinddx = zeros(nb_cols)
    @unpack x, dx = p
    
    for i=1:nb_rows # i indexes time
        w = s.u[i,:][1]
        wind = wind_DC_xt.(x, s.t[i], Ref(p))
        dwinddx = dwind_dx_DC_xt.(x, s.t[i], Ref(p))
        for j=1:nb_cols-1 # j indexes location in x
            dwdx[j] = (w[j+1] - w[j]) / dx
            #dwinddx[j] = (wind[j+1] - wind[j]) / dx
        end
        dwdx[nb_cols] = (w[1] - w[nb_cols]) / dx
        #dwinddx[nb_cols] = (wind[1] - w[nb_cols]) / dx

        A[i] = wind .* dwdx
        C[i] = w .* dwinddx
        AnC[i] = A[i] .+ C[i]
    end
    return A, C, AnC
end