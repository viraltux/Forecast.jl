```@setup examples
using Forecast
using Plots; gr()
Plots.reset_defaults()
```
# Seaborne (Multivariate)

## Introduction
In this example we will predict seaborne deadweight imports for France, Germany and the United Kingdom with daily data from 2015-04-01 to 2021-05-02. Below we can see the last 100 days of data.

```@example examples
plot(seaborne()[end-100:end,:]) #hide
```
## Descriptive Analysis
```@example examples
using StatsPlots
sb = seaborne()
@df sb corrplot(cols(1:3), grid = false)
```
The data is quite random in apperance showing a much larger seaborne average of imports for the UK than for France and Germany, whoever there seems to be a small positive correlation among the threea countries which gives hope for a multivariate analysis.


	


