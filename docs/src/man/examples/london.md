```@setup examples
using Forecast
using Plots; gr()
Plots.reset_defaults()
default(size=(800,500))
using PrettyTables
```
# London Weather/Crime (Multivariate)

## Introduction
In this example we will study the relationship between weather and crime in in Greater London. In particular we will be considering ten years of monthly data about weather and crime in Greater London from 2008 to 2018.


```@example examples
x = london()
plot(x)
```

The variables at play will be:

* **MaxTemp**:  Monthly data for the mean daily maximum temperature in C° for each month.
* **Violence**: Violent crimes comprising; Assault with Injury, Common Assault, Grievous Bodily Harm, Harassment, Murder, Offensive Weapon, Other violence, Wounding/GBH
* **Sexual**: Rape, Other sexual crimes.

## Descriptive Analysis
### Seasonality
```@example examples
splot(x[:,["Date","MaxTemp"]], title="Seasonal Maximum Temperature")
```
```@example examples
splot(x[:,["Date","Violence"]], title="Seasonal Violent Crimes")
```
```@example examples
splot(x[:,["Date","Sexual"]], title="Seasonal Sexual Crimes")
```
Sesonality patterns for `MaxTemp` are obvious and noticeable for `Violence`, however there's barely any to be seen in `Sexual`.

### Normalization
To compare time series with different positive magnitudes and to facilitate an stationary forecasting we will normalize data by dividing each variable by its maximum value.
```@example examples
tx = copy(x)
tx[!,2:end] = Array(x[:,2:end]) ./   (1. * mapslices(maximum, Array(x[:,2:end]), dims=[1]))
tx = tx[:,["Date","MaxTemp","Violence","Sexual"]]
# hide
```
To achieve an stationary time series to start our analysis we need to differentiate `MaxTemp` with a lag of 12, however this might not be enough for `Violence` and `Sexual` since they also show a marked trend from 2013 onwards.

```@example examples
dtx = d(tx,1,12) 
plot(dtx)
```

An extra differentation might be required for `Violence` and `Sexual` but this seems to be one to many for `MaxTemp` indicated by an increase in its variance when doing so. Also, the apparent need for a differentation in `Violence` and `Sexual` seems to be a temporary effect rather than a permanent trend or seasonality therefore we will continue with just one differentation of order 12.

## Fitting an AR(1) model
Given the scarcity of data for a multivariate model we cannot fit too many parameters while keeping its significance. We will therefore fit an AR of order one.

```@example examples
xar = ar(dtx,1,false)
show(xar) # hide
```
The first thing we can do is to use Φ1 to infer directional causality, in this case we see how `MaxTemp` is not influenced at all (with statistical significance) by `Violence` and `Sexual`, which could not possibly be otherwise in this case. On the other hand we have `Violence` and `Sexual` being affected significantly by `MaxTemp`.

The relathionship between `Violence` and `Sexual` is both ways, however, the effect `Sexual` has on `Violence` is small and barely significant which might prompt us to consider to remove it altogether, the opposite though is not true, `Violence` has a strong and significant influence in `Sexual`.


## Fixing Coefficients
Let's now fix to 0 the influcence of `Violence` and `Sexual` on `MaxTemp`, remove the influence from `MaxTemp` on `Sexual` since it is accounted for via `Violence` and the non-significant influence of `Sexual` on `Violence`. If now we fit again the model we have

```@example examples
dϕ1 = xar.Φ
dϕ2 = copy(dϕ1)
dϕ2[1,2:3,1] .= 0.0
dϕ2[3,1,1] = 0.0
dϕ2[2,3,1] = 0.0
dΦ = (dϕ1, dϕ2)
fxar = ar(dtx,1,false; dΦ)
show(fxar) # hide
```

## Transformed Forecasting
Next we can see the two years forecast for the transformed dataset.

```@example examples
dtfc = forecast(fxar,12*2)
plot(dtfc)
```
In order to obtain a forecast for the original data we need to invert the transformations carried out previously.

### Scaled Forecasting

Next we can see the two years forecast for the scaled dataset.

```@example examples
x0 = Array(tx[1:12,2:end])
x0 = reshape(x0',1,3,12)
tfc = p(dtfc, x0)
setnames!(tfc,["MaxTemp / $(maximum(x.MaxTemp))",
               "Violence / $(maximum(x.Violence))",
               "Sexual / $(maximum(x.Sexual))"])
plot(tfc,title="Greater London Crime/Weather Scaled Forecast")
```
The two years forecast follows:

### Final Forecast
```@example examples
fc = transform(tfc,(v->v*maximum(x.MaxTemp)),1);
fc = transform(fc,(v->round(v*maximum(x.Violence))),2)
fc = transform(fc,(v->round(v*maximum(x.Sexual))),3)
setnames!(fc,["MaxTemp","Violence","Sexual"])

fc.mean
```


