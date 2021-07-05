@userplot SPlot

"""
Package: Forecast

    splot(x, labels)

Plot a seasonal plot of x considering the parameter `labels`

# Arguments
- `x`: regular timed observations
- `labels`: This parameter accepts Integer, String and Vector values. When an Integer the labels are 1:labels, when a Vector the labels are specified within and when a String it accepts values "month", "day" and "quarter" expecting the first value of x to fall in "Jan", "Mon" or "Q1" unless x is a DataFrame in which case it is treated as a Time Series where the first Date typed column and value columns ares considered, observations are then automatically ordered either by "month", "day" or "quarter" and labels may be use to rename the default values.

# Returns
Sesonal plot

# Examples
```julia-repl
julia> splot(rand(120),"month")
julia> splot(rand(120),"quarter")
julia> splot(rand(120),"day")
```
"""
splot()

@recipe function f(sp::SPlot)

    day_labels = ["Mon","Tue","Wed","Thu","Fri","Sat","Sun"]
    month_labels = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]
    quarter_labels = ["Q1","Q2","Q3","Q4"]

    na = min(length(sp.args),2)

    x = sp.args[1]

    # DataFrame
    if isa(x, DataFrame)
        x = tots(x)
        t = x[:,1]
        x = x[:,2]

        dt = Δt(t)
        
        if dt == Day(1)
            labels = (na == 2) ? sp.args[2] : day_labels
            season = 7
            rp = dayofweek(t[1]) - 1
        elseif dt == Month(1)
            labels = (na == 2) ? sp.args[2] : month_labels
            season = 12
            rp = month(t[1]) - 1
        else
            labels = (na == 2) ? sp.args[2] : quarter_labels
            season = 4
            rp = quarterofyear(t[1]) - 1
        end

        # Shift values to assign correct labels
        x = vcat(repeat([missing],rp),x)
    end

    # Numerical Series
    if (na == 1) & !isa(sp.args[1], DataFrame)
        labels = 1:12
        season = 12
    end

    if (na == 2) & !isa(sp.args[1], DataFrame)

        if isa(sp.args[2], String)
            if sp.args[2] == "month"
                labels = month_labels
                season = 12

            elseif sp.args[2] == "day"
                labels = day_labels
                season = 7
                
            elseif sp.args[2] == "quarter"
                labels = quarter_labels
                season = 4
            else
                @error "Only 'month', 'day' and 'quarter' are valid string values for `labels`"
            end
        end

        if isa(sp.args[2], Vector)
            labels = sp.args[2]
            season = length(sp.args[2])
        end

        if isa(sp.args[2], Integer)
            season = sp.args[2]
            labels = 1:season
        end
        
    end
    
    # set up the subplots
    legend --> false
    grid --> false

    ls = ceil(length(x)/season)

    # padding to assure correct size
    x = vcat(x,repeat([missing],Int(season*ls-length(x))))

    # labels for seasons 
    xticks := (ls/2:ls+1:length(x)+ls, labels)
    
    x = Array(reshape(x,:,Int(length(x)/season))')

    # season plots
    xs = vcat(x, Array(repeat([missing],inner=season)') )
    xs = reshape(xs,:,1)

    # season means
    m = mapslices(mean ∘ skipmissing, x, dims=1)
    m = repeat(m, size(x)[1] )
    m = vcat(m, Array(repeat([missing],inner=season)') )
    m = reshape(m,:,1)
    
    @series begin
        seriestype := :line
        xs
    end

    @series begin
        seriestype := :line
        m
    end
    
end
