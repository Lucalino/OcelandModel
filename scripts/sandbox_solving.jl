#using NLsolve
#using Roots
using GLMakie
using OrdinaryDiffEq
using DynamicalSystems
#using Colors

#solw1 = roots(W1 -> 2 * u * (w0 - W1) + Lo1 * eo - Lo1 * exp(a * (W1/Wsat - b)), W1_intv)

# #system of equations
# function f!(fvec, W)
#     fvec[1] = 2 * (w0 - W[1]/Lo1) * u + Lo1 * (eo - precip_from_water_content(W[1], Lo1, Wsat, a, b))
#     fvec[2] = 2 * (2 * W[1]/Lo1 - w0 - W[2]/Li) * u + Li * (land_evap(W[4], spwp, sfc, Ep) - precip_from_water_content(W[2], Li, Wsat, a, b))
#     fvec[3] = 2 * (2 * W[2]/Li - 2 * W[1]/Lo1 - W[3]/Lo2 + w0) * u + Lo2 * (eo - precip_from_water_content(W[3], Lo2, Wsat, a, b))
#     fvec[4] = precip_from_water_content(W[2], Li, Wsat, a, b) * infiltration(W[4], ϵ, r) - land_evap_(W[4], spwp, sfc, Ep)
# end

# #Jacobian
# function j!(J, W)
#     J[1, 1] = -2 * u / Lo1 - Lo1 * precip_diff_from_water_content(W[1], Lo1, Wsat, a, b)
#     J[1, 2] = 0.
#     J[1, 3] = 0.
#     J[1, 4] = 0.
#     J[2, 1] = 4 * u / Lo1
#     J[2, 2] = -2 * u / Li - Li * precip_diff_from_water_content(W[2], Li, Wsat, a, b)
#     J[2, 3] = 0. 
#     J[2, 4] = Li * land_evap_diff(W[4], spwp, sfc, Ep) 
#     J[3, 1] = - 4 * u / Lo1 
#     J[3, 2] = 4 * u / Li 
#     J[3, 3] = - 2 * u / Lo2 - Lo2 * precip_diff_from_water_content(W[3], Lo2, Wsat, a, b)
#     J[3, 4] = 0.
#     J[4, 1] = 0. 
#     J[4, 2] = infiltration(W[4], ϵ, r) * precip_diff_from_water_content(W[2], Li, Wsat, a, b) 
#     J[4, 3] = 0. 
#     J[4, 4] = precip_from_water_content(W[2], Li, Wsat, a, b) * infiltration_diff(W[4], ϵ, r) - land_evap_diff(W[4], spwp, sfc, Ep)
# end

# #W[1], W[2], W[3] in [mm^2] and W[4] = s dimensionless

# #sol = nlsolve(f!, j!, [50. *Lo1, 50. * Li, 50. * Lo2, 0.8], xtol = 1e-5)

# #Visualising Eq. 1
# #W1 = range(10. * Lo1, length= 100, stop = 80. * Lo1)
# #plot(W1, 2 .* (w0 .- W1./Lo1) .* u + Lo1 .* (eo .- exp.(a.*(W1 ./ (Wsat .* Lo1) .- b))))


# #original definition of functions to be solved, problem: I somehow need to hand over the
# #dictionary that contains the variables
# #I found an inconsistency in these equations which arrised because I wasn't clear about
# #whether W1 should be water vapour pass or water content. The below form with Wsat = saturated
# #water vapour pass in [mm], is not correct for either definition. (9.8.2021)

# function f_one( (W1, W2, W3, s) )
#     return SVector(2 * u * (w0 - W1) + Lo1 * eo - Lo1 * exp(a * (W1/Wsat - b)), 
#                     2 * u * (2 * W1 - W2 - w0) - Li * exp(a * (W2/Wsat - b)),
#                     2 * u * (2 * W2 - 2 * W1 - W3 + w0) + Lo2 * eo - Lo2 * exp(a * (W3/Wsat - b)),
#                     exp(a * (W2/Wsat - b)) * infiltration(s, ϵ, r)
#                     )
# end

# function f_two( (W1, W2, W3, s) )
#         return SVector(2 * u * (w0 - W1) + Lo1 * eo - Lo1 * exp(a * (W1/Wsat - b)), 
#                     2 * u * (2 * W1 - W2 - w0) + Li * Ep/(sfc - spwp) * (s - spwp) - Li * exp(a * (W2/Wsat - b)),
#                     2 * u * (2 * W2/Li - 2 * W1/Lo1 - W3/Lo2 + w0) + Lo2 * eo - Lo2 * exp(a * (W3/Wsat - b)),
#                     exp(a * (W2/Wsat - b)) * infiltration(s, ϵ, r) - Ep/(sfc - spwp) * (s - spwp)
#                     )
# end

# function f_three( (W1, W2, W3, s) )
#     return SVector(2 * u * (w0 - W1) + Lo1 * eo - Lo1 * exp(a * (W1/Wsat - b)), 
#                    2 * u * (2 * W1 - W2 - w0) + Li * Ep - Li * exp(a * (W2/Wsat - b)),
#                    2 * u * (2 * W2 - 2 * W1 - W3 + w0) + Lo2 * eo - Lo2 * exp(a * (W3/Wsat - b)),
#                    exp(a * (W2/Wsat - b)) * infiltration(s, ϵ, r) - Ep
#                     )
# end



# #om_eq_solution() function with one s-interval and subsequent case checking

# function om_eq_solution(f1, f2, f3, d, W1_lower=0.0, W1_upper=100.0, W2_lower=0.0, W2_upper=100.0, W3_lower=0.0, W3_upper=100.0, s_lower=0.0, s_upper=1.0)

#     @unpack spwp, sfc, Ep, eo, ϵ, r, nZr, w0, L, Li, Lo1, Lo2, u, a, b, Wsat = d
    
#     #search intervals
#     W1_intv = W1_lower..W1_upper
#     W2_intv = W2_lower..W2_upper
#     W3_intv = W3_lower..W3_upper
#     s_intv  = s_lower..s_upper

#     #find roots with IntervalRootFinding.jl in the three soil moisture regimes
    
#     input_functions = (f1, f2, f3)

#     #anonymous function that reduces f1, f2, f3 to functions with one argument 
#     #as required by roots() by giving d and t as constants
#     #The first element is equivalent to: closure_one(u) = f1(u, d, 0.0)
#     closures = [variables -> f(variables,d,0.0) for f in input_functions] 
#     sol_dry = roots(closures[1], W1_intv × W2_intv × W3_intv × s_intv)
#     sol_trans = roots(closures[2], W1_intv × W2_intv × W3_intv × s_intv)
#     sol_wet = roots(closures[3], W1_intv × W2_intv × W3_intv × s_intv)


#     #initialise matrix of solutions
#     sol_matrix = Array{Float64}(undef, 0, output_columns)

#     #validate solutions

#     found_solution = false
#     more_than_one_solution_dry = false
#     more_than_one_solution_trans = false
#     more_than_one_solution_wet = false

#     if (length(sol_dry) == 1 && isunique(sol_dry[1]) == true)
#         sol_dry = mid.(interval.(sol_dry))[1]
#         if sol_dry[4] < spwp
#             sol_dry_W1, sol_dry_W2, sol_dry_W3, sol_dry_s = rou(sol_dry[1],2), rou(sol_dry[2],2), rou(sol_dry[3],2), rou(sol_dry[4],2)
#             found_solution = true
#             sol_matrix = [sol_matrix; sol_dry_W1 sol_dry_W2 sol_dry_W3 sol_dry_s]
#         end
#     elseif length(sol_dry) > 1
#         more_than_one_solution_dry = true
#     end

#     if (length(sol_trans) == 1 && isunique(sol_trans[1]) == true)
#         sol_trans = mid.(interval.(sol_trans))[1]
#         if sol_trans[4] >= spwp && sol_trans[4] <= sfc
#             sol_trans_W1, sol_trans_W2, sol_trans_W3, sol_trans_s = rou(sol_trans[1],2), rou(sol_trans[2],2), rou(sol_trans[3],2), rou(sol_trans[4],2)
#             found_solution = true
#             sol_matrix = [sol_matrix; sol_trans_W1 sol_trans_W2 sol_trans_W3 sol_trans_s]            
#         end
#     elseif length(sol_trans) > 1
#         more_than_one_solution_trans = true
#     end

#     if (length(sol_wet) == 1 && isunique(sol_wet[1]) == true)
#         sol_wet = mid.(interval.(sol_wet))[1]
#         if sol_wet[4] > sfc
#             sol_wet_W1, sol_wet_W2, sol_wet_W3, sol_wet_s = rou(sol_wet[1],2), rou(sol_wet[2],2), rou(sol_wet[3],2), rou(sol_wet[4],2)
#             found_solution = true
#             sol_matrix = [sol_matrix; sol_wet_W1 sol_wet_W2 sol_wet_W3 sol_wet_s]     
#         end
#     elseif length(sol_wet) > 1
#         more_than_one_solution_wet = true
#     end

#     return sol_matrix

#     # P_o_mean = (Lo1 * precip(sol_trans_W1, Wsat, a, b) + Lo2 * precip(sol_trans_W3, Wsat, a, b)) / (Lo1 + Lo2)
#     # precip_ratio =  precip(sol_trans_W2, Wsat, a, b) / P_o_mean
#     # println("Precipitation ratio: ", rou(precip_ratio,4))
#     # println("P_1: ", rou(precip(sol_trans_W1, Wsat, a, b),2))
#     # println("P_2: ", rou(precip(sol_trans_W2, Wsat, a, b),2))
#     # println("P_3: ", rou(precip(sol_trans_W3, Wsat, a, b),2))
#     # println("E_l: ", rou(land_evap(sol_trans_s, spwp, sfc, Ep),2))
# end


################################## BASINS ISSUE ###################################

function my_system(x, p, t)
    dx = (exp(0.22 * x[2] - 9.4) * (1.0 - x[1]^2) - 2.2 * tanh(10.0 * x[1] - 4.5) - 2.2) / 100.0
    dy = 2.2 * tanh(10.0 * x[1] - 4.5) + 2.2 - exp(0.22 * x[2] - 9.4) + (x[3] - x[2]) * 4.32
    dz = 3.0 - exp(0.22 * x[3] - 9.4) - (x[3] - x[2]) * 0.48
    return SVector(dx, dy, dz)
end

p = Dict(
)

xg = range(-0.1, 1.0; length = 100)
yg = zg = range(0.0, 72.0; length = 100)
grid = (xg, yg, zg)

# Initial condition that is predicted to diverge but doesn't
x0 = xg[60]
y0 = yg[90]
z0 = zg[90] 

u0 = [x0, y0, z0]
ds = ContinuousDynamicalSystem(my_system, u0, p)

diffeq = (alg = Vern9(), adaptive = false, dt = 0.001, reltol = 1e-8, abstol = 1e-8)

# xg = range(-0.1, 1.0; length = 30)
# yg = zg = range(0.0, 72.0; length = 50)
# grid = (xg, yg, zg)
# x0 = xg[1]; y0 = yg[1]; z0 = zg[1]
# u0 = [x0, y0, z0]
# ds = ContinuousDynamicalSystem(my_system, u0, p)
# diffeq = (alg = Vern9(), reltol = 1e-8, abstol = 1e-8)
# basins, attractors = basins_of_attraction(grid, ds; diffeq = diffeq)



tr = trajectory(ds, 50)

x, y, z = columns(tr)
t = collect(1:1:length(x))

fig = Figure()
ax1 = Axis(fig[1,1])
ax2 = Axis(fig[2,1])
ax3 = Axis(fig[3,1])
lines!(ax1, t, x)
lines!(ax2, t, y)
lines!(ax3, t, z)
ax1.ylabel = "x"
ax2.ylabel = "y"
ax3.xlabel = "t"
ax3.ylabel = "z"

basins, attractors = basins_of_attraction(grid, ds)

function cm_basins_plot(b, g1, g2)
    #cmap = colormap("Blues"; logscale=false)
    fig = Figure()
    ax = Axis(fig[1,1])
    hm = heatmap!(ax, g1, g2, b)
    Colorbar(fig[1, 2], hm)
    return fig
end