```@setup examples
using Forecast
using Plots; gr()
Plots.reset_defaults()
default(size=(800,500))
```
# Quakes (Univariate)

## Introduction
In this example we will predict the of number of earthquakes per year with a magnitude higher or equal to six. The data for the analysis has been collected from [USGS](https://earthquake.usgs.gov/) and aggregated from 1950 to 2020.

```@example examples
plot(quakes()) #hide
```
## Descriptive Analysis
First let's use a few utilities contained in the Forecast package to have a first impression on the data:

### Numerical Summary
```@example examples
qk = quakes()
sqk = summarize(qk)
show(sqk) # hide
```
Data's behavior seems to follow a Normal distribution with no strong indications of seasonal patterns in its plot.

### Autoregressive Behavior
```@example examples
plot(pacf(qk.quakes))
```
The partial autoregression function shows us that there seems to be a significant correlation to the number of earthquakes taking place last year. If we were looking for seasonality we could check on periods of 11 or 15 years since they show a nearly significant correlations but since they're most likely spureous (...or are they?) we will ignore them in this analysis.

## Fitting an AR Model
```@example examples
ar_qk = ar(qk)
show(ar_qk) # hide
```
In the AR model of order one we have highly significant coefficients and increasing its order does not provide important changes in the Information Criteria, however, the residuals show a barely significant normality behavior and we may consider to transform our data to improve on that. 
Given tha large noise in the model, tranformations to improve results will not be dramatic and therefore we will continue with a simple AR model of order one for our forecasting.


## Forecasting Earthquakes
```@example examples
fc_qk = forecast(ar_qk,10);
plot(fc_qk)
```
The plot shows us the forecast for the next ten years and, as we see, the large noise in the model does not allow us to be very accurate in our forecasting, but at least we can confidently say that there is a resonable chance to have a larger number of big earthquakes in 2021 and 2022 than the number we had in 2020.

```@example examples
show(fc_qk) # hide
```
