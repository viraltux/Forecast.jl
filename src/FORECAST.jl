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

"""
Transform a FORECAST object value with given function
"""
function transform(fc::FORECAST, f::Function)

    ts_x = tots(fc.model.x)
    names_x = names(ts_x)
    names_x[2:end] = "$(string(f))_" .* names_x[2:end]

    ts_mean = tots(fc.mean)
    ts_lower = tots(fc.lower)
    ts_upper = tots(fc.upper)
    names_mean = names(ts_mean)
    names_mean[2:end] = "$(string(f))_" .* names_mean[2:end]
    names_lower = names(ts_lower)
    names_lower[2:end] = "$(string(f))_" .* names_lower[2:end]
    names_upper = names(ts_upper)
    names_upper[2:end] = "$(string(f))_" .* names_upper[2:end]

    xts = ts_x[:,1]
    x = f.(Array(ts_x[:,2:end]))

    fmean = f.(Array(ts_mean[:,2:end]))
    flower = f.(Array(ts_lower[:,2:end]))
    fupper = f.(Array(ts_upper[:,2:end]))

    n = size(x,1)
    m = size(x,2)

    pfc = deepcopy(fc)
    pfc.model.varnames = names_x[2:end]

    pfc.model.x = hcat(ts_x[:,1:1], DataFrame(x,names_x[2:end]))
    pfc.mean    = hcat(ts_mean[:,1:1], DataFrame(fmean,names_mean[2:end]))
                        

    if size(fmean,2) > 1
        z = fupper .- repeat(fmean,1,2)
        fmean = repeat(pfmean,1,2)
    else
        z = fupper .- fmean
    end
    
    pfc.upper = hcat(ts_mean[:,1:1],
                     DataFrame(fmean .+ z,names_upper[2:end]))
    pfc.lower = hcat(ts_mean[:,1:1],
                     DataFrame(fmean .- z,names_lower[2:end]))
    
    pfc.call = fc.call * "\nData transformed with function: $(f)"

    return(pfc)
end
