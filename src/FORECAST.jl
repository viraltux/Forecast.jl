"""
Package: Forecast

Store results from the function `forecast`

# Arguments
   model:  Model object containing information about the fitted model.
   x:      Original time series.
   alpha:  The confidence levels associated with the prediction intervals.
   mean:   Point forecasts.
   lower:  Lower limits for prediction intervals.
   upper:  Upper limits for prediction intervals.
"""
mutable struct FORECAST
    model
    alpha::Tuple
    mean
    upper
    lower
    call
end

function Base.show(io::IO, fc::FORECAST)
    ts = eltype(fc.mean[:,1]) in [Date, DateTime, Time,
                                  Year, Month, Quarter, Week, Day,
                                  Hour, Second, Millisecond,
                                  Microsecond, Nanosecond]
    printstyled("Forecast Information\n",bold=true,color=:underline)
    print("\n    ",fc.call, "\n")
    printstyled("\nMean Forecasting\n",bold=true,color=:underline)
    pretty_table(fc.mean, noheader = false, nosubheader = true, show_row_number=false)
    printstyled("\nPrediction Intervals alpha at: "*string(fc.alpha)*"\n",bold=true,color=:underline)
    printstyled("\nUpper:\n",color=:underline)
    pretty_table(fc.upper, noheader = false, nosubheader = true, show_row_number=false,
                 vlines = [0,1,div(size(fc.upper,2),2)+1,size(fc.upper,2)])
    printstyled("\nLower:\n",color=:underline)
    pretty_table(fc.lower, noheader = false, nosubheader = true, show_row_number=false,
                 vlines = [0,1,div(size(fc.lower,2),2)+1,size(fc.lower,2)])
end


