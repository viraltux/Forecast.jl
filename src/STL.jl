"""
Package: Forecast

Store results from the function stl

# Arguments
    `ta::TimeArray`    A time array with three time series from a fitted STL model
    `call::String`     method called to generate ta
""" 
mutable struct STL{T<:TimeArray}
    ta::T         # A time array with three time series from a fitted STL model
    call::String  # method called to generate ta
end

