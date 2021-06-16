module Forecast

using CSV, Distributions, ColorSchemes, DataFrames, DataFramesMeta, Dates, GZip,
    HypothesisTests, LinearAlgebra, Optim, Plots, PrettyTables, RecipesBase, Statistics, StatsBase

# types
export AR, CCF, STL

# methods
export acf, ar, arsim, boxcox, ccf, d, hma, iboxcox, loess, p, pacf, forecast, 
       setnames!, sma, stl, transform, summarize

# datasets
export air, co2, london, quakes, seaborne

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
include("boxcox.jl") 
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

Featured Methods:

    acf:        Auto-correlation or auto-covariance of univariate series. 
    ar:         Multivariate Autoregressive Model.
    arsim:      Simulated Multivariate Autoregressive Model.
    boxcox:     Box-Cox Transformations.
    ccf:        Cross-correlation or cros-covariance of two univariate series.
    d:          Lagged differences of a given order for Vector and Array.
    forecast:   Forecast values of fitted time series models.
    hma:        Henderson moving average filter.
    iboxcox:    Inverse Box-Cox Transformations.
    loess:      Locally estimated scatterplot smoothing.
    p:          Reverse lagged differences of a given order for types Vector and Array.
    pacf:       Partial Auto-correlation function.
    sma:        Simple moving average.
    splot:      Plot a seasonal plot for types Vector and TimeArray.
    stl:        Seasonal and Trend decomposition using loess.
    summarize:  Statistical summary.
"""
Forecast

end
