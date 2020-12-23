module Forecast

using CSV, DataFrames, LinearAlgebra, Plots, RecipesBase, TimeSeries

# types
export STL
# methods
export loess, sma, stl

# source files
include("loess.jl")
include("sma.jl")
include("stl.jl")
include("utils.jl") 

# recipes
include("plotrecipes.jl")

"""
Collection of methods for Time Series analysis

Methods implemented:

    loess:      Locally weighted smoothed series.
    sma:        Simple moving average.
    stl:        Seasonal and Trend decomposition using loess.
"""
Forecast

end
