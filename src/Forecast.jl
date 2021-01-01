module Forecast

using CSV, Distributions, DataFrames, LinearAlgebra, Plots, RecipesBase, TimeSeries, Statistics

# types
export STL
# methods
export acf, ccf, d, loess, sma, stl
# datasets
export co2

# source files
include("acf.jl")
include("ccf.jl")
include("d.jl")
include("loess.jl")
include("sma.jl")
include("stl.jl")

include("utils.jl")

include("datasets.jl")

## recipes
include("plotrecipes.jl")

"""
Collection of methods for Time Series analysis

Methods implemented:

    acf:        Auto-correlation or auto-covariance of a univariate serie. 
    ccf:        Cross-correlation or cros-covariance of two univariate series.
    d:          Lagged differences of a given order for Vector, Array and TimeSeries.
    loess:      Locally weighted smoothed series.
    sma:        Simple moving average.
    stl:        Seasonal and Trend decomposition using loess.
"""
Forecast

end
