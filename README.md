# Forecast [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://viraltux.github.io/Forecast.jl/dev)

Julia package containing utilities intended for Time Series analysis.

> :warning: This package is in an early development stage and its functionality has not been thoroughly tested. Please, consider to report issues if you find any.

## Featured Methods

* Autocorrelated models for univariate and multivariate data (ar)
* Autocorrelated model simulations for univariate and multivariate data (arsim)
* Autocorrelated model forecasting with custom plots (forecast)
* Autocorrelation & Autocovariance function and custom plots (acf)
* Crosscorrelation & Crosscovariance function and custom plots (ccf)
* Henderson Moving Average Filter (hma)
* Lagged differences and Reverse lagged differences of a given order (d | p)
* Locally Estimated Scatterplot Smoothing (loess)
* Seasonal and Trend decomposition based on Loess and custom plot (stl) 
* Seasonal Plot for series and Time Series (splot)
* Simple Moving Average (sma)

<img src="./docs/src/images/stl_readme.png">

## Datasets

* Atmospheric Carbon Dioxide Dry Air Mole Fractions from quasi-continuous measurements at Mauna Loa, Hawaii.

K.W. Thoning, A.M. Crotwell, and J.W. Mund (2020), Atmospheric Carbon Dioxide Dry Air Mole Fractions from continuous measurements at Mauna Loa, Hawaii, Barrow, Alaska, American Samoa and South Pole. 1973-2019, Version 2020-08 National Oceanic and Atmospheric Administration (NOAA), Global Monitoring Laboratory (GML), Boulder, Colorado, USA https://doi.org/10.15138/yaf1-bk21 FTP path: ftp://aftp.cmdl.noaa.gov/data/greenhouse_gases/co2/in-situ/surface/


* Estimates of world seaborne trade from AIS data collected by MarineTraffic; available at UN COMTRADE Monitor.

Cerdeiro, Komaromi, Liu and Saeed (2020). Subset with imports and exports data for France, Germany and the United Kingdom from 2015-04-01 to 2021-05-02.  https://comtrade.un.org/data/ais

## References

* [Cleveland et al. 1990]  Cleveland,  R.  B.;  Cleveland,  W.  S.;McRae, J. E.; and Terpenning, I.  1990.  STL: A seasonal-trend decomposition procedure based on loess. Journal of Official Statistics 6(1):3–73.


[![Coverage](https://codecov.io/gh/viraltux/Forecast.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/viraltux/Forecast.jl)


* 


