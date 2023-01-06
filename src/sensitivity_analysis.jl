function information_entropy(x::Vector, y::Vector, nb_bins::Int)
    ds = Dataset(x,y)
    H = genentropy(ds, VisitationFrequency(RectangularBinning(nb_bins)); q = 1.0, base = 2)
    return H
end


function mutual_information_binlength(x::Vector, y::Vector, bin_length = 0.1)
    xnorm = normalise(x)
    ynorm = normalise(y)
    Hx  = genentropy(probabilities(Dataset(xnorm), bin_length))
    Hy  = genentropy(probabilities(Dataset(ynorm), bin_length))
    Hxy = genentropy(probabilities(Dataset(xnorm,ynorm), bin_length))
    m = Hx + Hy - Hxy
    return m
end

function mutual_information_nbbins(x::Vector, y::Vector, nb_bins::Int = 10)
    Hx = genentropy(probabilities(Dataset(x), nb_bins))
    Hy = genentropy(probabilities(Dataset(y), nb_bins))
    Hxy = genentropy(probabilities(Dataset(x,y), nb_bins))
    m = Hx + Hy - Hxy
    return m
end

function bootstrapping_binlength(x::Vector, y::Vector, bin_length = 0.1, N::Int = 1000)
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

function bootstrapping_nbbins(x::Vector, y::Vector, nb_bins = 10, N::Int = 1000)
    null  = zeros(N)
    Hx  = genentropy(probabilities(Dataset(x), nb_bins))
    Hy  = genentropy(probabilities(Dataset(y), nb_bins))
    for i in 1:N
        shuffle!(x)
        shuffle!(y)
        Hxy = genentropy(probabilities(Dataset(x, y), nb_bins))
        null[i] = Hx + Hy - Hxy
    end
    μ = mean(null)
    std = Statistics.std(null)
    three_sigma = 3 * std
    return null, μ, three_sigma
end


function relative_mi_binlength(x::Vector, y::Vector, bin_length = 0.1, N = 1000)
    mi = mutual_information_binlength(x, y)
    null, μ, three_sigma = bootstrapping_binlength(x,y)
    mi_rel = mi / (μ + three_sigma)
    return mi_rel
end

function relative_mi_nbbins(x::Vector, y::Vector, nb_bins = 10, N = 1000)
    mi = mutual_information_nbbins(x, y, nb_bins)
    null, μ, three_sigma = bootstrapping_nbbins(x,y, nb_bins)
    mi_rel = mi / (μ + three_sigma)
    return mi_rel
end


function all_parameter_sensitivities(data::DataFrame, yquant::String, filename::String, nbbins = 10, om = false, τ = true, lin = false)
    if lin == false   
        if τ == true

            if om == false
                p = ["α", "b", "ep", "spwp", "nZr", "ϵ", "a", "wsat", "τ", "eo", "r"]
            else
                p = ["α", "b", "ep", "spwp", "nZr", "w0", "ϵ", "a", "wsat", "τ", "u", "L", "eo", "r"]
            end

            rmi = zeros(length(p))
            for i = 1:length(p)
                rmi[i] = relative_mi_nbbins(data[:,p[i]], data[:,yquant], nbbins)
            end

        else
            println("Only implemented for τ dataset so far.")
        end
    else
        if om == false
            p = ["α", "τ", "p", "e", "r", "eo"]
            rmi = zeros(length(p))
            for i = 1:length(p)
                rmi[i] = relative_mi_nbbins(data[:,p[i]], data[:,yquant])
            end
        else
           @warn "Closed model version not yet implemented."
        end
    end
    df = DataFrame(pnames = p, MI_rel = rmi)
    CSV.write(datadir("mutual information", filename * ".csv"), df)
    return df
end