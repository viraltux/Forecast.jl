function Δt(ts::AbstractVector{<:Date}) 
    ts1 = ts[1]
    ts2 = ts[2]
    
    dt = Dates.value(ts2-ts1) # Days
    
    dt in 365:366 && return Year(1)
    dt in 90:92 && return Quarter(1)
    dt in 28:31 && return Month(1)
    dt % 7 == 0 && return Week(dt)
    return Day(dt)
end

function Δt(ts::AbstractVector{<:DateTime}) 
    ts1 = ts[1]
    ts2 = ts[2]
    
    dt = Dates.value(ts2-ts1) # Milliseconds
    
    dt % (24*3_600_000) == 0 && return Δt(convert(AbstractVector{Date},ts))
    dt % 3_600_000 == 0 && return Hour(dt/3_600_000)
    dt % 60_000 == 0 && return Minute(dt/60_000)
    dt % 1_000 == 0 && return Second(dt/1_000)
    return Millisecond(dt)
end

function Δt(ts::AbstractVector{<:Time}) 
    ts1 = ts[1]
    ts2 = ts[2]
    
    dt = Dates.value(ts2-ts1) # Nanoseconds

    dt % 3_600_000_000_000 == 0 && return Hour(dt/3_600_000_000_000)
    dt % 60_000_000_000 == 0 && return Minute(dt/60_000_000_000)
    dt % 1_000_000_000 == 0 && return Second(dt/1_000_000_000)
    dt % 1_000_000 == 0 && return Millisecond(dt/1_000_000)
    dt % 1_000 == 0  && return Microsecond(dt/1_000)
    return Nanosecond(dt)
end

function Δt(ts::AbstractVector)
    ts[2]-ts[1]
end


# Next n increments for a timestamp vector
function nΔt(ts::AbstractVector, n::Integer)
    dt = Δt(ts)
    ts1 = n > 0 ? ts[end] + dt : ts[1] - dt 
    n > 0 && return reshape(collect(ts1 : dt : ts1 + (n-1)*dt),:,1)
    n < 0 && return reshape(collect(ts1 + (n+1)*dt: dt : ts1),:,1)
    return ts
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

