function Δt(ts::AbstractVector{<:Date}) 

    @assert length(ts) >= 3 "time series must have a length >= 3"
    
    ts1 = ts[1]
    ts2 = ts[2]
    ts3 = ts[3]

    dy = year(ts[2])-year(ts[1])
    dm = month(ts[2])-month(ts[1])
    dd = day(ts[2])-day(ts[1])

    same_day = day(ts1) == day(ts2) == day(ts3)
    same_month = month(ts1) == month(ts2) == month(ts3)
    same_year = year(ts1) == year(ts2) == year(ts3)
    
    !same_day && dd % 7 == 0 && return Week(dd÷7)
    !same_day && return Day(dd)
    same_day && !same_month && dm % 4 == 0 && return Quarter(dm÷4)
    same_day && !same_month && return Month(dm)
    same_day && !same_year && return Year(dy)
    return ts[2]-ts[1]
    
end

function Δt(ts::AbstractVector{<:DateTime}) 

    @assert length(ts) >= 3 "time series must have a length >= 3"
    
    ts1 = ts[1]
    ts2 = ts[2]
    ts3 = ts[3]
    
    dt = Dates.value(ts2-ts1) # Milliseconds
    
    dt % (24*3_600_000) == 0 && return Δt(convert(AbstractVector{Date},ts))
    dt % 3_600_000 == 0 && return Hour(dt/3_600_000)
    dt % 60_000 == 0 && return Minute(dt/60_000)
    dt % 1_000 == 0 && return Second(dt/1_000)
    return ts[2]-ts[1]

end

function Δt(ts::AbstractVector{<:Time}) 

    @assert length(ts) >= 3 "time series must have a length >= 3"
    
    ts1 = ts[1]
    ts2 = ts[2]
    ts3 = ts[3]
    
    dt = Dates.value(ts2-ts1) # Nanoseconds

    dt % 3_600_000_000_000 == 0 && return Hour(dt/3_600_000_000_000)
    dt % 60_000_000_000 == 0 && return Minute(dt/60_000_000_000)
    dt % 1_000_000_000 == 0 && return Second(dt/1_000_000_000)
    dt % 1_000_000 == 0 && return Millisecond(dt/1_000_000)
    dt % 1_000 == 0  && return Microsecond(dt/1_000)
    return Nanosecond(dt)
    
end

Δt(ts::AbstractVector) = ts[2]-ts[1]

# Next n increments for a timestamp vector
function nΔt(ts::AbstractVector, n::Integer)
    dt = Δt(ts)
    ts1 = n > 0 ? ts[end] + dt : ts[1] - dt 
    n > 0 && return reshape(collect(ts1 : dt : ts1 + (n-1)*dt),:,1)
    n < 0 && return reshape(collect(ts1 + (n+1)*dt: dt : ts1),:,1)
    return reshape(ts,:,1)
end


# DataFrame to canon Time Series Dataframe
function tots(df::DataFrame, interval = Day; start::Integer = 1)::DataFrame

    x = df[:,eltype.(eachcol(df)) .<: Union{Missing,Real}]
    names_x = names(x)
    ts = df[:,eltype.(eachcol(df)) .<: Union{Date, DateTime, Time,
                                             Year, Quarter, Month, Day,
                                             Hour, Second, Millisecond,
                                             Microsecond, Nanosecond}]
    if size(ts) != (0,0)
        names_x = vcat(names(ts)[:,1],names_x)
        return df[!,names_x]
    else
        ts = reshape(collect(interval(start):interval(1):interval(start-1+size(x,1))),:,1)
        tdf = DataFrame(ts,["Time"])
        return hcat(tdf,df[!,names_x])
    end
    
end

function tots(ar::Matrix, interval = Day; start::Integer = 1)::DataFrame
    df = DataFrame(ar,:auto)
    tots(df,interval; start)
end

function tots(ar::Vector, interval = Day; start::Integer = 1)::DataFrame
    df = DataFrame(reshape(ar,:,1),:auto)
    tots(df,interval; start)
end

