module Forecast

using CSV, Distributions, DataFrames, GZip, LinearAlgebra, Plots,
      RecipesBase, TimeSeries, Statistics

# types
export AR, CCF, STL

# methods
export acf, ar, arsim, ccf, d, hma, loess, p, pacf, predict, sma, stl

# datasets
export co2

# types
include("AR.jl")
include("CCF.jl")
include("STL.jl")

# source files
include("acf.jl")
include("ar.jl")
include("arsim.jl") 
include("ccf.jl")
include("d.jl")
include("datasets.jl")
include("hma.jl") 
include("loess.jl")
include("p.jl")
include("pacf.jl")
include("predict_ar.jl")
include("sma.jl")
include("stl.jl") 
include("utils.jl")

# plot recipes
##  type
include("CCFplot.jl")
include("STLplot.jl")

## user
include("splot.jl")

"""
Collection of methods for Time Series analysis

Methods implemented:

    acf:        Auto-correlation or auto-covariance of univariate series. 
    ar:         Multivariate Autoregressive Model.
    arsim:      Simulated Multivariate Autoregressive Model.
    ccf:        Cross-correlation or cros-covariance of two univariate series.
    d:          Lagged differences of a given order for Vector, Array and TimeSeries.
    hma:        Henderson moving average filter.
    loess:      Locally estimated scatterplot smoothing.
    p:          Reverse lagged differences of a given order for types Vector, Array and TimeArray.
    pacf:       Partial Auto-correlation function.
    sma:        Simple moving average.
    splot:      Plot a seasonal plot for types Vector and TimeArray.
    stl:        Seasonal and Trend decomposition using loess.
"""
Forecast

end
