```@setup examples
using Forecast
using Plots; gr()
Plots.reset_defaults()
default(size=(800,500))
```
# WHAT-IF Scenarios

## Introduction
There are situations in multivariate analysis in which we want to know things like:

* How robust is our model to impacts
* How does our model evolves given known events
* How does overall accuracy improves if we improve forecast accuracy for specific variables
* How does our model reacts to different scenarios

To answer this questions and similar ones we need a way to instruct our forecast to test this scenarios.

## London Crime/Weather Scenarios

Let's show how to use the Great London Crime/Weather example. In this example we see how weather temperature affects/explains criminality in Greater London. However, the prediction of temperature is not as accurate as it would be if we use forecast coming from (metoffice)[metoffice.gov.uk] instead those predicted by an AR(1) model. Therefore in order to have more accurate predictions about crime we want to overrule the temperature forecast of the AR(1) model with the forecast comming from `metoffice`.

These metoffice predictions would not be much different than those from the AR(1) model therefore, in order to better visualize impacts in the forecast, we will play with an hypothetical scenario in which temperatures increase one degree per month for two years straight.

Also, let's imagine that the accuracy of the one degree per year scenario is much better than the one coming for the AR(1) model making its variance 25 times smaller.

## Fitting and Forecast

Let's first fit our model and retrieve its forecast. We will be using the scaled forecast to better see how our changes affect the model.

```@example examples
x = london()
tx = copy(x)
tx[!,2:end] = Array(x[:,2:end]) ./   (1. * mapslices(maximum, Array(x[:,2:end]), dims=[1]))
tx = tx[:,["Date","MaxTemp","Violence","Sexual"]]

dtx = d(tx,1,12) 

xar = ar(dtx,1,false)

dϕ1 = xar.Φ
dϕ2 = copy(dϕ1)
dϕ2[1,2:3,1] .= 0.0
dϕ2[3,1,1] = 0.0
dϕ2[2,3,1] = 0.0
dΦ = (dϕ1, dϕ2)
fxar = ar(dtx,1,false; dΦ)

dtfc = forecast(fxar,12*2)

x0 = Array(tx[1:12,2:end])
x0 = reshape(x0',1,3,12)
tfc = p(dtfc, x0)
setnames!(tfc,["MaxTemp / $(maximum(x.MaxTemp))",
               "Violence / $(maximum(x.Violence))",
               "Sexual / $(maximum(x.Sexual))"])

plot(tfc,title="Greater London Crime/Weather Scaled Forecast")
```

## Scenario Building

The last forecast we ran is for `tfc`, therefore we need to transform as well our what-if scenario to be able to use it.

```@example examples
# Scenario: starting at 10 C° and increasing one degree per month for two years
witemp = vcat(x.MaxTemp, 10 .+ collect(1:2*12)) 
witemp = witemp ./ maximum(x.MaxTemp) 
witemp = d(witemp,1,12)[end-12*2+1:end]
# hide
```
Now we create a DataFrame with the values we want to fix in our scenario

```@example examples
fixMean = copy(tfc.mean)
fixMean[!,3:4] .= missing
fixMean[!,2] = witemp
# hide
```
We also fix the Variance matrix with a 25 times smaller variance for `MaxTemp` and update accordingly its covariance interactions.

```@example examples
fixΣ2  = copy(dtfc.model.Σ2)
fixΣ2[1,:] = fixΣ2[1,:]./5
fixΣ2[:,1] = fixΣ2[:,1]./5
# hide
```

## Forecasting Scenario

Finally we forecast the `fxar` model with the scenario built previously.

```@example examples
dtfc = forecast(fxar,12*2; fixMean, fixΣ2)
plot(dtfc,title="Greater London Crime/Weather Scaled Forecast & Fixed MaxTemp")
            
x0 = Array(tx[1:12,2:end])
x0 = reshape(x0',1,3,12)
tfc0 = p(dtfc)
tfc = p(dtfc, x0)


setnames!(tfc,["Fixed MaxTemp / $(maximum(x.MaxTemp))",
               "Violence / $(maximum(x.Violence))",
               "Sexual / $(maximum(x.Sexual))"])

plot(tfc,title="Greater London Crime/Weather Scaled Forecast & Fixed MaxTemp")
```
We can see how our scenario displays a much narrower forecast interval, however the new scenario in the already small covariance does not affect much the intervals for crime.
	
Another noticeable effect would be the inverse relationship between crime and temperature, this is something we can see directly into the AR(1) coefficients but now we can clearly visualize it.
