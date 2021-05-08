```@meta
CurrentModule = Forecast
```

# What is this about?

*It's about time*

Julia package containing utilities intended for Time Series analysis.

**Warning**: This package is under development and its functionality has not been thoroughly tested. Please, consider to report issues if you find any.

## Methods Available


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
	

## Manual Outline

```@contents
Pages = [
    "man/quickstart.md",
    "man/methods.md",
    "man/datasets.md",
    "man/docstrings.md",		
]
Depth = 1
```
