"""
Package: Forecast

    co2(full = false)

Return dataset with atmospheric Carbon Dioxide Dry Air Mole Fractions from quasi-continuous measurements at Mauna Loa, Hawaii.

K.W. Thoning, A.M. Crotwell, and J.W. Mund (2020), Atmospheric Carbon Dioxide Dry Air Mole Fractions from continuous measurements at Mauna Loa, Hawaii, Barrow, Alaska, American Samoa and South Pole. 1973-2019, Version 2020-08 National Oceanic and Atmospheric Administration (NOAA), Global Monitoring Laboratory (GML), Boulder, Colorado, USA https://doi.org/10.15138/yaf1-bk21 FTP path: ftp://aftp.cmdl.noaa.gov/data/greenhouse_gases/co2/in-situ/surface/

# Arguments
- `full`: if `true` Returns the full original dataset from 1973 to 2020 in a DataFrame, otherwise returns the subset used in "STL: A Seasonal-Trend Decomposition Procedure Based on Loess" from Cleveland et. al. Its default value is `false`.

# Returns
 Dataframe or TimeArray containing the descrived dataset.

# Examples
```julia-repl
julia> co2()
[ Info: Dataset used in Cleveland et al. paper
4609×1 TimeArray{Union{Missing, Float64},1,Date,Array{Union{Missing, Float64},1}}
1974-05-17 to 1986-12-31
│            │ co2     │
├────────────┼─────────┤
│ 1974-05-17 │ 333.38  │
│ 1974-05-18 │ 333.11  │
│ 1974-05-19 │ 333.46  │
   [...]
```
"""
function co2(full = false)

    data = "data/co2.csv.gz"
    path = joinpath(splitdir(@__DIR__)[1], data)
    
    co2_df = GZip.open(path, "r") do io
        CSV.read(io,DataFrame)
    end

    if full
        @info "Full dataset from 1973 to 2020"
        return co2_df
    else
        @info "Dataset used in Cleveland et al. paper"
    end

    # recode missing values
    co2_df.value = replace(co2_df.value, -999.99 => missing)

    dates_co2 = Dates.Date(1973,1,1):Dates.Day(1):Dates.Date(2019,12,31)
    co2_ts = TimeSeries.TimeArray(dates_co2, co2_df.value)

    # paper dates
    # using 17 may instead 17 april (as the paper claims)
    # since there seems to be no data for april in 1974
    dates_co2_stl = Dates.Date(1974,5,17):Dates.Day(1):Dates.Date(1986,12,31)

    # revome leap year day
    dates_co2_stl = filter(dates_co2_stl) do x
        (Dates.isleapyear(x) & (Dates.month(x) == 2)) ? Dates.day(x) != 29 : true
    end

    co2_stl_ts = co2_ts[dates_co2_stl]
    TimeSeries.rename!(co2_stl_ts,:co2)
    co2_stl_ts
    
end


"""
Package: Forecast

    seaborne(full = false)

Cerdeiro, Komaromi, Liu and Saeed (2020). Estimates of world seaborne trade from AIS data collected by MarineTraffic; available at UN COMTRADE Monitor.

By default a TimeArray containing deadweight imports for France, Germany and the United Kingdom from 2015-04-01 to 2021-05-02 is returned, otherwise a DataFrame is returned with import and exports wiht the below fields:

num_pc:     number of port calls
mtc:        metric tons of cargo
dwt:        deadweight tonnage
suffix_ma:  30-day moving averages


# Returns
 Dataframe or TimeArray containing the seaborne dataset.

# Examples
```julia-repl
julia> seaborne()
[ Info: Seaborne deadweight trade imports from AIS
2204×3 TimeArray{Int64, 2, Date, Matrix{Int64}} 2015-04-01 to 2021-04-12
│            │ Germany │ France  │ UK      │
├────────────┼─────────┼─────────┼─────────┤
│ 2015-04-01 │ 173686  │ 309178  │ 873765  │
│ 2015-04-02 │ 571960  │ 750539  │ 929067  │
   [...]
```
"""
function seaborne(full = false)

    data = "data/seaborne.csv.gz"
    path = joinpath(splitdir(@__DIR__)[1], data)
    
    sb_df = GZip.open(path, "r") do io
        CSV.read(io,DataFrame, dateformat = "mm/dd/yyyy HH:MM:SS pp", types = Dict(:date => Date))
    end

    if full
        @info "Full dataset with estimates of world seaborne trade from AIS"
        return sb_df
    else
        @info "Seaborne deadweight trade imports from AIS"
    end

    # Group by Country and Flow to select Imports
    gbt = groupby(sb_df, [:country_name, :flow])
    data = (datetime = Dates.Date(2015,4,1):Dates.Day(1):Dates.Date(2021,4,12),
            Germany = gbt[1].dwt,
            France = gbt[3].dwt,
            UK = gbt[5].dwt)

    TimeArray(data; timestamp = :datetime, meta = "Seaborne Deadweight Imports")

end


"""
Package: Forecast

    quakes()

Return the number of earthquakes per year on earth with a magnitude higher or equal to six since 1950 up to 2020.

# Examples
```julia-repl
julia> quakes()
71×1 TimeSeries.TimeArray{Int64, 1, Dates.Date, Vector{Int64}} 1950-01-01 to 2020-01-01
│            │ quakes │
├────────────┼────────┤
│ 1950-01-01 │ 138    │
│ 1951-01-01 │ 151    │
│ 1952-01-01 │ 181    │
   [...]
```
"""
function quakes()

    data = "data/quakes.csv.gz"
    path = joinpath(splitdir(@__DIR__)[1], data)
    
    qk_df = GZip.open(path, "r") do io
        CSV.read(io,DataFrame, dateformat = "yyyy", types = Dict(:year => Date))
    end

    TimeArray(qk_df,timestamp = :year)
    
end


