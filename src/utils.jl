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

function mean_of_bins!(df::DataFrame, param::String, var::String, nb_bins::Int)
    df = sort(df, param)
    col_name = string("mean_",var)
    df[!,col_name] .= 1.0
    nb_rows = nrow(df)
    mod(nb_rows, nb_bins) == 0 || error("Number of bins must divide the data set's number of rows")
    bin_size = Int(nb_rows / nb_bins)
    bin_start = 1  

    for bin_end in Iterators.countfrom(bin_size, bin_size)
        bin_end > nb_rows && break
        df[bin_start:bin_end, col_name] .= mean(df[bin_start:bin_end, var])
        bin_start = bin_end+1
    end
    return df[!,col_name]
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