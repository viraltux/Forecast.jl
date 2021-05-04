"""
Package: Forecast

Store results from the function `forecast`

# Arguments
   model:  Model object containing information about the fitted model.
   x:      Original time series.
   levels: The confidence values associated with the prediction intervals.
   mean:   Point forecasts.
   lower:  Lower limits for prediction intervals.
   upper:  Upper limits for prediction intervals.
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
    printstyled("Forecast Information\n",bold=true,color=:underline)
    print("\n    ",f.call, "\n")
    printstyled("\nMean Forecasting\n",bold=true,color=:underline)
    pretty_table(f.mean, noheader = true, nosubheader = true, show_row_number=false)
    printstyled("\nPrediction Intervals levels at: "*string(f.levels)*"\n",bold=true,color=:underline)
    printstyled("\nUpper:\n",color=:underline)
    pretty_table(f.upper, noheader = true, nosubheader = true, show_row_number=false,
                 vlines =0:2:size(f.mean,2)+size(f.mean,2))
    printstyled("\nLower:\n",color=:underline)
    pretty_table(f.lower, noheader = true, nosubheader = true, show_row_number=false,
                 vlines =0:2:size(f.mean,2)+size(f.mean,2))
end
