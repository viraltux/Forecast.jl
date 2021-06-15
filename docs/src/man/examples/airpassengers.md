```@setup examples
using Forecast
using DataFrames
using Plots; gr()
Plots.reset_defaults()
default(size=(800,500))
```
# Air Passengers (Seasonality)

## Introduction
In this example we will predict the monthly totals in thounsads of international airline passengers with data from 1949 to 1960.

```@example examples
plot(air()) #hide
```
## Descriptive Analysis & Transformations

### Heteroskedasticity
A cople of  obvious feature in the plot of the data is that it shows an non-stationary heteroskedastic time series. However, to facilitate the fitting of an AR model we want a Stationary Homoskedastic time series. With that goal in mind we could consider first a boxcox transformation to deal with the heteroskedasticity but in this case a simple log transformation works just fine.

```@example examples
log_ap = air()
log_ap.Passengers = log.(log_ap.Passengers)
rename!(log_ap, :Passengers => :LogPassengers)
plot(log_ap) #hide
```

### Stationary Behavior
We still have a clear non-stationary time series displaying a marked trend and montnly seasonality.
```@example examples
stl_ap = stl(log_ap.LogPassengers,12,robust=true)
plot(stl_ap)
```
To faciliate the fitting of an AR model we want an stationary time series, and in order to have one we will be using a differentation of order 12 for the trend and seasonality.
```@example examples
d12_log_ap = d(log_ap,1,12)
plot(d12_log_ap) # hide
```
There is a case for second differentation of oder one since it furthers reduces the variance of the resulting time series indicating a potential remaining trend, however this potential improvement might not be enough to justify this new transformation since the existing one offers now a good enough stationary time series to consider an AR model with a constant coefficient.

### Evaluating Seasonal Differentation

```@example examples
splot(d12_log_ap)
```
After a differentation of order 12 we can also see that the seasonality has mostly disappeared and we can continue our analysis with a reasonalble stationary dataset.

### Autoregressive Behavior

```@example examples
plot(pacf(d12_log_ap[:,2]))
```
We can see clear auto-correlations with values of the previous two months and with values of the same months in the previous year. Beyond this there seems to be not enough correlation to justify a more complex model, but we can always check this hypothesis lookig at the Information Criteria.

## Fitting an AR model

Let's then fit a AR with 13 parameters with a constant despite being differentiated and see how it looks:

```@example examples
ar_tap = ar(d12_log_ap,13)
show(ar_tap) # hide
```
### Fixing Coefficients
As expected we see that ``\Phi1,\Phi2,\Phi12,`` and ``\Phi13`` have highly significant coefficients, also we can see some significance in ``\Phi0``. This confirms the case for a further differentation of order one however, since doing so decreases the normality profile of the residuals we will keep it as it is.

Since we want to know the falues of these four parameters without the influcence of the rest we will now fit again the model fixing all coefficients except those five.

```@example examples
Φ = ar_tap.Φ
fΦ = copy(ar_tap.Φ)
fΦ[3:11] .= 0 # fixing coefficients from 3 to 11 to zero.
dΦ = (Φ,fΦ)   # Tuple informing AR which coefficients to fix.
arf_tap = ar(d12_log_ap,13;dΦ)
show(arf_tap) # hide
```

## Forecast of Transformed Series

```@example examples
fct = forecast(arf_tap,3*12)
plot(fct)
```

## Forecast Original Data

```@example examples
x0 = reshape(log_ap[:,2][1:12]',1,1,12)
fct2 = p(fct,x0)
fc = Forecast.transform(fct2,exp)
new_names = ["#Passengers"]
setnames!(fc,new_names)
plot(fc, title = "Forecast Air Passengers")
```
And the forecasted values in their prediction intervals for the next three years are:

```@example examples
show(fc) # hide
```

	

