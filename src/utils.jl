#length conversion
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


#time conversion
function s2day(val::Number)
    return val / (60 * 60 * 24)
end

function day2s(val::Number)
    return val * 60 * 60 * 24
end


#pressure conversion
function hPa2Pa(val::Number)
    return val * 10^2
end

function Pa2hPa(val::Number)
    return val / 10^2
end


#temperature conversion
function C2K(val::Number)
    return val + 273.25
end

function K2C(val::Number)
    return val - 273.25
end


#miscellaneous
function rou(val::Number,d::Integer)
    return round(val, digits=d)
end


function normalise(x::Vector)
    n = length(x)
    xnorm = zeros(n)
    for i = 1:n
        xnorm[i] = (x[i] - minimum(x)) / (maximum(x) - minimum(x))
    end
    return xnorm
end


function mean_of_bins!(df::DataFrame, param::String, var::String, nb_bins::Int)
    df = sort!(df, param)
    col_name = string("mean_",var)
    df[!,col_name] .= 1.0
    nb_rows = nrow(df)
    mod(nb_rows, nb_bins) == 0 || error("Number of bins must divide the data set's number of rows.")
    bin_size = Int(nb_rows / nb_bins)
    bin_start = 1  

    for bin_end in Iterators.countfrom(bin_size, bin_size)
        bin_end > nb_rows && break
        df[bin_start:bin_end, col_name] .= DataFrames.mean(df[bin_start:bin_end, var])
        bin_start = bin_end+1
    end
    return df[!,col_name]
end


function movingaverage(X::Vector, nb_of_elm::Int)
    BackDelta = div(nb_of_elm,2) 
    ForwardDelta = isodd(nb_of_elm) ? div(nb_of_elm,2) : div(nb_of_elm,2) - 1
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


