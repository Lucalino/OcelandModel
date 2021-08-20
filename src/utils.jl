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