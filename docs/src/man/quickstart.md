```@setup quickstart
using Plots; gr()
Plots.reset_defaults()
```
# Quick Start

## Installing Forecast.jl

To begin exploring the package functionality type the lines below in
the Julia REPL:

    Pkg> add https://github.com/viraltux/Forecast.jl

    julia> using Forecast

Alternatively you can also type:

    julia> using Pkg

    julia> Pkg.add(url="https://github.com/viraltux/Forecast.jl")

    julia> using Forecast


## LOESS extrapolated

The example below compares LOESS result with extrapolated predictions compared
to a simple moving average result with a window size of 100.

```@example quickstart
using Plots
using Forecast

n = 1000
axb = LinRange(-1/2,pi+1/2,n)
x = LinRange(0,pi,n)
y = sin.(x) .+ rand(n)
scatter(x, y, xlims=(-1/2,pi+1/2), ma=.5, label = "Data", color = :grey,size=(800,500))
plot!(axb,loess(x,y,predict=axb), linewidth = 4, label = "LOESS", color = :blue)
plot!(x,sma(y,100), linewidth = 2, label= "MA 100", color = :orange)
```

## STL on CO2 dataset

For this example we will be using the co2 data used by the creators of STL to
demostrate its funcitonality, below we can see such time series.

```@example quickstart
using Plots
using Forecast

plot(co2(), legend=:bottomright,size=(800,500))
```

The parameters used for STL they're also from the orginal paper, a period of
365 days is used (removing leap years extra day), a robust fit is required and
seasonality post-smoothing is applied.

```@example quickstart
using Plots
using Forecast

stl_co2 = stl(co2(),365; robust=true, spm=true)
plot(stl_co2,size=(800,1000))
```
The image below comes from the original paper for comparison purposes.

```@raw html
<img src="../../images/stl.png" width="800px"/>
```

## Cross-Correlation on shifted dasaset

Here we cross-correlate two identical series shifted by six positions, the plot
shows how the peak correlation is at position six.

```@example quickstart
using Plots
using Random
using Forecast

Random.seed!(36)
x1 = rand(100)
x2 = circshift(x1,6)
res = ccf(x1, x2; type="cor")
plot(res,size=(800,500))
```

## PACF on random dataset

The `pacf` function is useful to identify significant parameters in ARIMA models. For instance, in R the default `pacf` function estimates partial auto-correlation in a stepwise fashion, however in cases where the model is highly correlated with many previous steps this approach identifies the first lag as highly correlated and the rest as near zeroes when, in reality, all partial auto-correlations should be around zero since that's the information left once taken away the linear influence from the all other lags. Below is an example of such effect where the `stepwise` (in blue) and `real` (in red) partial auto-correlations are compared for a series where all lags highly correlate to each other.

```@example quickstart
using Plots
using Random
using Forecast

Random.seed!(36)
x = collect(1:100) + rand(100)
res = pacf(x)
plot(res,size=(800,500))
```

## Seasonal Plot on Seaborne dataset 

To compare seasonal behaviour we can use splot to display it side by side, in thi case it seems on Thursdays there is on average an higher trade activity in Germany's ports.

```@example quickstart
using Plots
using Forecast
sb = seaborne()
splot(sb.Germany)
```

## Multivariate Autoregressive Model

Random Walk simulated, fitted and forecast plotted.
```@example quickstart
using Plots
using Forecast
plot(forecast(ar(arsim( 1,0, 0,1,100),1,false),100), size=(800,500))
```

Random Zigzag Walk simulated, fitted and forecast plotted.
```@example quickstart
using Plots
using Forecast
plot(forecast(ar(arsim(-1,0,0,1,100),1),10), size=(800,500))
```

Bivariate dataset simulated, fitted and forecast plotted.
```@example quickstart
using Plots
using Forecast

Φ = [1 .1;
    -.2 1]
Φ0 = [2, 3]
x0 = [.1, .1]
Σ = [.2 .05;
     .05 .2]
ar_model = ar(arsim(Φ,Φ0,x0,Σ,100),1)
fc = forecast(ar_model,50)
plot(fc, size=(800,500))
```

Output from fitting a three dimensional seaborne trade dataset with two parameters.
```@example quickstart
using Plots
using Forecast

sb = seaborne()
ar_model = ar(sb,2)
show(ar_model)
```
