# This file contains unit conversion functions and other utilities for the Oceland Model.
# Created on 5th July 2021 by Luca Schmidt
# Last edited on 5th July 2021 by Luca Schmidt

function m2mm(val::Number)
    return val * 10^3
end

function mm2m(val::Number)
    return val * 10^(-3)
end

function km2mm(val::Number)
    return val * 10^6
end

function mm2km(val::Number)
    return val * 10^(-6)
end

function s2day(val::Number)
    return val / (60 * 60 * 24)
end

function day2s(val::Number)
    return val * 60 * 60 * 24
end

function rou(val::Number,d::Integer)
    return round(val, digits=d)
end

function mean_of_bins!(df::DataFrame, param::String, var::String, nb_bins::Int, name::String = "μ")
    df = sort!(df, param)
    col_name = string(name, "_", param)
    df[!,col_name] .= 1.0
    nb_rows = nrow(df)
    mod(nb_rows, nb_bins) == 0 || error("Number of bins must divide the data set's number of rows")
    bin_size = Int(nb_rows / nb_bins)
    bin_start = 1  

    for bin_end in Iterators.countfrom(bin_size, bin_size)
        bin_end > nb_rows && break
        df[bin_start:bin_end, col_name] .= DataFrames.mean(df[bin_start:bin_end, var])
        bin_start = bin_end+1
    end
    return nothing    #df[!,col_name]
end

function movingaverage(X::Vector, numofele::Int)
    BackDelta = div(numofele,2) 
    ForwardDelta = isodd(numofele) ? div(numofele,2) : div(numofele,2) - 1
    len = length(X)
    Y = similar(X)
    for n = 1:len
        if BackDelta > n-1
            lo = 1
            hi = n + (n-1)
        elseif n + ForwardDelta > len
            lo = n - (len-n)
            hi = len 
        else
            lo = n - BackDelta
            hi = n + ForwardDelta
        end
        Y[n] = DataFrames.mean(X[lo:hi])
    end
    return Y
end


function sensitivity_v1!(df::DataFrame, parameter::String, quantity::String = "PR", nb_bins::Int = 100)
    
    #compute μ for the quantity of interest in each bin of parameter values
    mean_of_bins!(df, parameter, quantity, nb_bins)

    #compute the square of the deviation from the mean for each data point
    df[:, string("dev^2_μ_", parameter)] = (df[:,quantity] .- df[:,string("μ_", parameter)]).^2
    
    #compute variance for each bin
    mean_of_bins!(df, parameter, string("dev^2_μ_", parameter), nb_bins, "σ^2")

    #compute standard deviation for each bin
    df[:, string("σ_", parameter)] = sqrt.(df[:, string("σ^2_", parameter)])

    σ_mean = DataFrames.mean(df[:,string("σ_", parameter)])

    return σ_mean
end


function sensitivity_v2!(df::DataFrame, parameter::String, quantity::String = "PR", nb_bins::Int = 100)

    #compute total mean of all data points
    μtot = DataFrames.mean(df[:, quantity])

    #compute difference between bin means and total mean
    mean_of_bins!(df, parameter, quantity, nb_bins)
    df[:, string("Δμ_", parameter)] = abs.(df[:, string("μ_", parameter)] .- μtot)
    #df[:, string("Δμ_", parameter)] = df[:, string("μ_", parameter)] .- μtot
    meandiff = DataFrames.mean(df[:, string("Δμ_", parameter)])

    return meandiff
end

function information_entropy(x::Vector, y::Vector, nb_bins::Int)
    ds = Dataset(x,y)
    H = genentropy(ds, VisitationFrequency(RectangularBinning(nb_bins)); q = 1.0, base = 2)
    return H
end

function mutual_information(x::Vector, y::Vector, bin_length = 0.1)
    xnorm = normalise(x)
    ynorm = normalise(y)
    Hx  = genentropy(Dataset(xnorm), bin_length)
    Hy  = genentropy(Dataset(ynorm), bin_length)
    Hxy = genentropy(Dataset(xnorm,ynorm), bin_length)
    m = Hx + Hy - Hxy
    return m
end

function normalise(x::Vector)
    n = length(x)
    xnorm = zeros(n)
    for i = 1:n
        xnorm[i] = (x[i] - minimum(x)) / (maximum(x) - minimum(x))
    end
    return xnorm
end

function bootstrapping(x::Vector, y::Vector, bin_length = 0.1, N::Int = 10000)
    xnorm = normalise(x)
    ynorm = normalise(y)
    null  = zeros(N)
    Hx  = genentropy(Dataset(xnorm), bin_length)
    Hy  = genentropy(Dataset(ynorm), bin_length)
    for i in 1:N
        shuffle!(xnorm)
        shuffle!(ynorm)
        Hxy = genentropy(Dataset(xnorm, ynorm), bin_length)
        null[i] = Hx + Hy - Hxy
    end
    μ = mean(null)
    std = Statistics.std(null)
    three_sigma = 3 * std
    return null, μ, three_sigma
end

function relative_mi(x::Vector, y::Vector, bin_length = 0.1, N = 10000)
    mi = mutual_information(x, y)
    null, μ, three_sigma = bootstrapping(x,y)
    mi_rel = mi / (μ + three_sigma)
    return mi_rel
end

function all_parameter_sensitivities(data::DataFrame, yquant::String, filename::String, om = false, τ = true)
    if τ == true
        if om == false
            p = ["α", "b", "ep", "spwp", "nZr", "sfc", "ϵ", "a", "wsat", "τ", "eo", "r"]
        else
            p = ["α", "b", "ep", "spwp", "nZr", "w0", "ϵ", "a", "wsat", "τ", "u", "L", "eo", "r"]
        end
        rmi = zeros(length(p))
        for i = 1:length(p)
            rmi[i] = relative_mi(data[:,p[i]], data[:,yquant])
            #println("$(p[i])" * " and " * yquant * " have relative mutual information mi_rel = $(rmi[i]).")
        end
    else
        println("Only implemented for τ dataset so far.")
    end
    #size = size(data)[1]
    df = DataFrame(pnames = p, MI_rel = rmi)
    CSV.write(datadir("sims", "mutual information", filename * ".csv"), df)
    return df
end




function cm_load_data()
    s1000 = CSV.read(datadir("sims", "closed model pmscan/cm_smooth_eq_MC_fixedpoints_10000_runs_domain1000_all_quantities" * ".csv"), DataFrame)
    s10000 = CSV.read(datadir("sims", "closed model pmscan/cm_smooth_eq_MC_fixedpoints_10000_runs_domain10000_all_quantities" * ".csv"), DataFrame)
    return s1000, s10000
end

function cm_load_alphavar_data()
    df = CSV.read(datadir("sims", "closed model pmscan/cm_smooth_variable-α_bunt_eq_MC_fixedpoints_5000_runs_all_quantities.csv"), DataFrame)
    return df
end

function om_v2_load_data()
    s1000 = CSV.read(datadir("sims", "open model pmscan/om_v2_MC_fixedpoints_runs10000_domain1000_all_quantities" * ".csv"), DataFrame)
    s10000 = CSV.read(datadir("sims", "open model pmscan/om_v2_MC_fixedpoints_runs10000_domain10000_all_quantities" * ".csv"), DataFrame)
    return s1000, s10000
end

function om_v2_closed_load_data()
    #s1000 = CSV.read(datadir("sims", "open model pmscan/om_v2_closed_MC_fixedpoints_runs10000_domain1000_all_quantities" * ".csv"), DataFrame)
    s10000 = CSV.read(datadir("sims", "open model pmscan/om_v2_closed_MC_fixedpoints_runs10000_domain10000_all_quantities" * ".csv"), DataFrame)
    return s10000
end