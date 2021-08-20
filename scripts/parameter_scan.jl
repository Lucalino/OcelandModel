using DrWatson
@quickactivate "Oceland Model"
include(srcdir("parametrisations.jl"))
include(srcdir("list_of_combinations.jl"))

using DataFrames
using PyPlot
using CSV
pygui(true)

comptype = "load data"
filepath = "sims/closed_model_pmscan_highres.csv" #inside datadir

if comptype == "create data"
    #parameters
    spwp = collect(range(0.1, length=3, stop=0.3))
    sfc  = collect(range(0.7, length=3, stop=0.9))
    Ep   = collect(range(4.0, length=3, stop=4.5))
    s    = collect(range(0.0, length=100, stop=1.0))
    Ptot = collect(range(2.5, length=3, stop=3.5))
    α    = collect(range(0.3, length=3, stop=0.7))
    ϵ    = collect(range(0.9, length=3, stop=1.1))
    r    = collect(range(1.9, length=3, stop=2.1))

    #Create array of all points in the parameter space
    param = [spwp, sfc, Ep, s, Ptot, α, ϵ, r]
    pmpoints = listofpoints(param)

    #Number of points
    nbpoints = length(param)

    #Function that evaluates PR at one point in the parameter space
    function precipratio(spwp, sfc, Ep, s, Ptot, α, ϵ, r)
        El = land_evapo(spwp,sfc,Ep,s)
        phi  = infiltration(ϵ,r,s)
        #precipitation ratio
        PR = El * (1-α) / (Ptot * phi - α * El)
        return PR
    end

    #Function that evaluates dPR at one point in the parameter space
    function dPR(spwp, sfc, Ep, s, Ptot, α, ϵ, r)
        El = land_evapo(spwp,sfc,Ep,s)
        phi  = infiltration(ϵ,r,s)

        #partial derivatives of the precipitation ratio
        denom = (El * α - Ptot * phi)^2
        dPR_El = Ptot * phi * (1-α) / denom
        dPR_Ptot = -El * phi * (1-α) / denom
        dPR_Φ = -El * Ptot * (1-α) / denom
        dPR_α = El * (El - Ptot * phi) / denom

        dPR_ϵ = dPR_Φ * (-1) * s^r
        dPR_r = dPR_Φ * (-1) * ϵ * s^r * log(s)
        
        if s <= spwp
            dPR_Ep   = 0.0
            dPR_sfc  = 0.0
            dPR_spwp = 0.0
            dPR_s    = 0.0
        elseif s <= sfc
            dPR_Ep   = dPR_El*(s-spwp)/(sfc-spwp)
            dPR_sfc  = dPR_El*Ep*(spwp-s)/(sfc-spwp)^2
            dPR_spwp = dPR_El*Ep*(s-sfc)/(sfc-spwp)^2
            dPR_s    = dPR_El*Ep/(sfc-spwp)+dPR_Φ*(-1)ϵ*r*s^(r-1)
        else
            dPR_Ep   = dPR_El
            dPR_sfc  = 0.0 
            dPR_spwp = 0.0
            dPR_s    = dPR_Φ*(-1)ϵ*r*s^(r-1)
        end 
        dPR_tot = dPR_spwp + dPR_sfc + dPR_Ep + dPR_s + dPR_Ptot + dPR_α + dPR_ϵ + dPR_r
        der    = [dPR_spwp dPR_sfc dPR_Ep dPR_s dPR_Ptot dPR_α dPR_ϵ dPR_r dPR_tot]   
        return der
    end


    data = Array{Float64}(undef, 0, 18)

    for i=1:length(pmpoints)
        pr = precipratio(pmpoints[i]...)
        der = dPR(pmpoints[i]...)
        datarow = [transpose(pmpoints[i]) pr der]
        data = [data;datarow]
    end

    df = convert(DataFrame, data)
    datanames = [:spwp, :sfc, :Ep, :s, :Ptot, :α, :ϵ, :r, :PR, :dPR_spwp, :dPR_sfc, :dPR_Ep, :dPR_s, :dPR_Ptot, :dPR_α, :dPR_ϵ, :dPR_r, :dPR_tot]
    rename!(df, datanames)
    CSV.write(datadir("sims/closed_model_pmscan_highres.csv"),df,delim='\t')

elseif comptype == "load data"
    
    #reading the file 
    data = DataFrame(CSV.File(datadir(filepath)))
    
    #filter out the points after divergence for which PR is negative
    data = filter(row -> row.PR >= 0.0,data)

    #filter out the rows in which s and alpha give the largest contribution to the total derivative
    #I verified that none of the other parameters is dominant at any point in the parameter space.

    sdom = filter(row -> abs(row.dPR_s) > max(abs(row.dPR_spwp), abs(row.dPR_sfc), abs(row.dPR_Ep), abs(row.dPR_Ptot), abs(row.dPR_α), abs(row.dPR_ϵ), abs(row.dPR_r)),data)
    αdom = filter(row -> abs(row.dPR_α) > max(abs(row.dPR_spwp), abs(row.dPR_sfc), abs(row.dPR_Ep), abs(row.dPR_Ptot), abs(row.dPR_s), abs(row.dPR_ϵ), abs(row.dPR_r)),data)

end


