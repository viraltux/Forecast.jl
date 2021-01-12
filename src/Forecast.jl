module Forecast

using CSV, Distributions, DataFrames, GZip, LinearAlgebra, Plots,
      RecipesBase, TimeSeries, Statistics

# types
export CCF, STL

# methods
export acf, ccf, d, hma, loess, pacf, sma, stl

# datasets
export co2

# types
include("CCF.jl")
include("STL.jl")

# source files
include("acf.jl") 
include("ccf.jl")
include("d.jl")
include("datasets.jl")
include("hma.jl") 
include("loess.jl")
include("pacf.jl")
include("sma.jl")
include("stl.jl") 
include("utils.jl")

# recipes
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
