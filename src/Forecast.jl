module Forecast

using CSV, Distributions, ColorSchemes, DataFrames, DataFramesMeta, Dates, GZip,
    HypothesisTests, LinearAlgebra, Plots, PrettyTables, RecipesBase,
    Statistics, StatsBase

# types
export AR, CCF, STL

# methods
export acf, ar, arsim, ccf, d, hma, loess, p, pacf, forecast, sma, stl, summarize

# datasets
export co2, seaborne, quakes

# types
include("AR.jl")
include("CCF.jl")
include("FORECAST.jl")
include("STL.jl")
include("SUMMARIZE.jl")

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
include("forecast_ar.jl")
include("sma.jl")
include("stl.jl")
include("summarize.jl") 
include("utils.jl")
include("utils_datetime.jl")

# plot recipes
##  type
include("plot_CCF.jl")
include("plot_DataFrame.jl")
include("plot_FORECAST.jl")
include("plot_STL.jl")

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
    forecast:   Forecast values of fitted time series models.
    hma:        Henderson moving average filter.
    loess:      Locally estimated scatterplot smoothing.
    p:          Reverse lagged differences of a given order for types Vector, Array and TimeArray.
    pacf:       Partial Auto-correlation function.
    sma:        Simple moving average.
    splot:      Plot a seasonal plot for types Vector and TimeArray.
    stl:        Seasonal and Trend decomposition using loess.
    summarize:  Statistical summary.
"""
Forecast

end
