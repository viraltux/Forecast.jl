module Forecast

using CSV, Distributions, DataFrames, GZip, LinearAlgebra, Plots,
      RecipesBase, TimeSeries, Statistics

# types
export CCF, STL
# methods
export acf, ccf, d, loess, sma, stl, hma, hmaSymmetricWeights
# datasets
export co2

# source files

## ccf
include("ccf.jl")
include("acf.jl")

include("d.jl")

## stl 
include("loess.jl")
include("sma.jl")
include("stl.jl")

## hma
include("hma.jl")

include("utils.jl")

include("datasets.jl")

## recipes
include("plotrecipes.jl")

"""
Collection of methods for Time Series analysis

Methods implemented:

    acf:        Auto-correlation or auto-covariance of univariate series. 
    ccf:        Cross-correlation or cros-covariance of two univariate series.
    d:          Lagged differences of a given order for Vector, Array and TimeSeries.
    hma:        Henderson moving average filters.
    loess:      Locally estimated scatterplot smoothing.
    sma:        Simple moving average.
    stl:        Seasonal and Trend decomposition using loess.
"""
Forecast

end
