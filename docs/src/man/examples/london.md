```@setup examples
using Forecast
using Plots; gr()
Plots.reset_defaults()
```
# Weather/Crime (Multivariate)

## Introduction
In this example we will study the relationship between weather and crime in in Greater London. In particular we will be considering ten years of monthly data about weather and crime in Greater London from 2008 to 2018.

The variables at play will be:

**MaxTemp**:  Monthly data for the mean daily maximum temperature in C° for each month.
**Violence**: Violent crimes comprising; Assault with Injury, Common Assault, Grievous Bodily Harm, Harassment, Murder, Offensive Weapon, Other violence, Wounding/GBH
**Sexual**: Rape, Other sexual crimes.


```@example examples
plot(london()) #hide
```
## Descriptive Analysis
```@example examples
using StatsPlots
sb = seaborne()
@df sb corrplot(cols(1:3), grid = false)
```
The data is quite random in apperance showing a much larger seaborne average of imports for the UK than for France and Germany, whoever there seems to be a small positive correlation among the threea countries which gives hope for a multivariate analysis.
