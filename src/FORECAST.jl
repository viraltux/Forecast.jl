"""
Package: Forecast

Store results from the function `forecast`

# Arguments
   model: Model object containing information about the fitted model.
       x: Original time series.
   level: The confidence values associated with the prediction intervals.
    mean: Point forecasts.
   lower: Lower limits for prediction intervals.
   upper: Upper limits for prediction intervals.
"""
mutable struct FORECAST
    model
    levels::Tuple
    mean
    upper
    lower
    call
end

function Base.show(io::IO, f::FORECAST)
    println("Prediction Intervals levels at: ", string(f.levels))
    println()
    
    # display(DataFrame((mu=f.mean, lower1=f.lower[:,1], upper1=f.upper[:,1],
    #                                lower2=f.lower[:,2], upper2=f.upper[:,2])))
end
