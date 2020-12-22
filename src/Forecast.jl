module Forecast

using LinearAlgebra, TimeSeries

# types
export STL
# methods
export loess, sma, stl

# source files
include("utils.jl")
include("sma.jl")
include("loess.jl")

"""
Collection of methods for Time Series analysis

Methods implemented:

    loess:      Locally weighted smoothed series.
    sma:        Simple moving average.
    stl:        Seasonal and Trend decomposition using loess.
"""
Forecast

end
