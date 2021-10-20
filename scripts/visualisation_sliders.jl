
using DrWatson
@quickactivate "Oceland Model"
include(srcdir("parametrisations.jl"))
include(srcdir("utils.jl"))
include(srcdir("model_versions.jl"))
using InteractiveDynamics
using DynamicalSystems, OrdinaryDiffEq
import GLMakie

spwp = 0.3
sfc  = spwp + 0.3
ep   = 4.38
pt   = 10
eo   = 3.0
ϵ    = 1.0
r    = 2.0
α    = 0.2
nZr  = 100.0
a    = 15.6
b    = 0.603
wsat = 72.0
u    = 5.0 * km2mm(1)/s2day(1)
L    = 1000 * km2mm(1) 

p = @dict spwp sfc ep pt eo ϵ r α nZr a b wsat u L

diffeq = (alg = Vern9(), adaptive = false, dt = 0.001, reltol = 1e-8, abstol = 1e-8)

ds = ContinuousDynamicalSystem(closed_model_smooth, [0.2, 50.0, 50.0], p)

ps = Dict(
    :spwp  => 0.2:0.01:0.54,
    :eo => 2.5:0.1:3.5,
    :ϵ  => 0.9:0.01:1.1,
    :r  => 1.9:0.01:2.1,
    :α  => 0.1:0.01:0.5,
    :nZr  => 90.0:1.0:110.0,
    :a  => 11.4:0.1:15.6,
    :b  => 0.52:0.01:0.61,
    :wsat  => 65.0:1.0:80.0,
    :u => 5.0* km2mm(1)/s2day(1):1.0* km2mm(1)/s2day(1):15.0* km2mm(1)/s2day(1),
    :L => 1000.0* km2mm(1):1000.0* km2mm(1):10000* km2mm(1), 
)


u0s =  [ds.u0, [0.3, 60.0, 50.0], [0.5, 50.0, 60.0], [0.2, 35.0, 30.0]]

fig, obs = interactive_evolution_timeseries(
    ds, u0s, ps; tail = 10000, idxs = (1, 2, 3), diffeq,
    lims =((0.0,1.0), (0.0, 100.0), (0.0, 100.0))
)

ax = GLMakie.content(fig[1,1])
ax.xlabel = "s"
ax.ylabel = "wl"
ax.zlabel = "wo"