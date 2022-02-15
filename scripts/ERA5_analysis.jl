using DrWatson
@quickactivate "Oceland Model"
using ClimateBase


data1 = datadir("era", "era5_monthly_av_sl_2020.nc")
data2 = datadir("era", "era5_monthly_av_pl_2020.nc")
p = ncread(data1,"tp") #total precip 
wvp = ncread(data1,"tcwv") #total column water vapour
rwp = ncread(data1,"tcrw") #total column rain water
u = ncread(data2, "u")
v = ncread(data2, "v")
wind = sqrt.(u.^2 .+ v.^2)
meanwind = ClimateBase.mean(wind, dims = Pre)


function rolling_average(data::ClimArray, len::Int, month::Int)
    bin_mean_data = deepcopy(data)

    if length(size(data)) == 4

        for i in 1:size(bin_mean_data)[1]
            i_start = i < len ? 1 : i-len+1 
            i_end   = (i+len) < size(bin_mean_data)[1] ? i + len : size(bin_mean_data)[1]
            for n in 1:size(bin_mean_data)[2]
                n_start = n < len ? 1 : n-len+1
                n_end   = (n+len) < size(bin_mean_data)[2] ? n + len : size(bin_mean_data)[2]
                bin_mean_data[i, n, 1, month] = ClimateBase.mean(bin_mean_data[i_start:i_end, n_start:n_end, 1, month])
            end
        end

    elseif length(size(data)) == 3

        for i in 1:size(bin_mean_data)[1]
            i_start = i < len ? 1 : i-len+1 
            i_end   = (i+len) < size(bin_mean_data)[1] ? i + len : size(bin_mean_data)[1]
            for n in 1:size(bin_mean_data)[2]
                n_start = n < len ? 1 : n-len+1
                n_end   = (n+len) < size(bin_mean_data)[2] ? n + len : size(bin_mean_data)[2]
                bin_mean_data[i, n, month] = ClimateBase.mean(bin_mean_data[i_start:i_end, n_start:n_end, month])
            end
        end
    end

    return bin_mean_data
end