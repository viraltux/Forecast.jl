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

## LOESS example

```@example tutorial
using Plots
using Forecast
n = 1000
axb = LinRange(-1/2,pi+1/2,n)
x = LinRange(0,pi,n)
y = sin.(x) .+ rand(n)
scatter(x, y, xlims=(-1/2,pi+1/2), ma=.5, label = "Data", color = :grey)
plot!(axb,loess(x,y,predict=axb), linewidth = 4, label = "LOESS", color = :blue)
plot!(x,sma(y,100), linewidth = 2, label= "MA 100", color = :orange)
```

## STL example

### CO2 data

```@example tutorial
using Plots
using Forecast
plot(co2())
```

### CO2 Decomposition
```@example tutorial
using Plots
using Forecast
stl_co2 = stl(co2(),365; robust=true, spm=true)
plot(stl_co2)
```

## Cros-Correlation example
```@example tutorial
x1 = rand(100);
x2 = circshift(x1,6);
res = ccf(x1, x2; type="cor");
plot(res)
```

