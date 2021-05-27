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


function rename!(fc::FORECAST,new_names)
    n_mean = names(fc.mean)
    n_mean[2:end] = new_names
    rename!(fc.mean,n_mean)

    nn_upper = []
    for n in new_names
        push!(nn_upper,"upper1_"*n)
        push!(nn_upper,"upper2_"*n)
    end
    n_upper = names(fc.upper)
    n_upper[2:end] = nn_upper
    rename!(fc.upper,n_upper)
    
    nn_lower = []
    for n in new_names
        push!(nn_lower,"lower1_"*n)
        push!(nn_lower,"lower2_"*n)
    end
    n_lower = names(fc.lower)
    n_lower[2:end] = nn_lower
    rename!(fc.lower,n_lower)

    if fc.model.x isa DataFrame
        n_model = names(fc.model.x)
        n_model[2:end] = new_names
        rename!(fc.model.x,n_model)
        fc.model.varnames=new_names
    end
end
