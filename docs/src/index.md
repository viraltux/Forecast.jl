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
	boxcox:     Box-Cox transformation.
    d:          Lagged differences of a given order.
    forecast:   Forecast values of fitted time series models.
    hma:        Henderson moving average filter.
    loess:      Locally estimated scatterplot smoothing.
    p:          Reverse lagged differences of a given order.
    pacf:       Partial Auto-correlation function.
	rename!:    Renames FORECAST objects' variables.
    sma:        Simple moving average.
    splot:      Plot a seasonal plot.
    stl:        Seasonal and Trend decomposition using loess.
    summarize:  Statistical summary.
