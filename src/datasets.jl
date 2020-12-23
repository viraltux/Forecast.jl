"""
co2(full = false)

Return dataset with atmospheric Carbon Dioxide Dry Air Mole Fractions from quasi-continuous measurements at Mauna Loa, Hawaii.

K.W. Thoning, A.M. Crotwell, and J.W. Mund (2020), Atmospheric Carbon Dioxide Dry Air Mole Fractions from continuous measurements at Mauna Loa, Hawaii, Barrow, Alaska, American Samoa and South Pole. 1973-2019, Version 2020-08 National Oceanic and Atmospheric Administration (NOAA), Global Monitoring Laboratory (GML), Boulder, Colorado, USA https://doi.org/10.15138/yaf1-bk21 FTP path: ftp://aftp.cmdl.noaa.gov/data/greenhouse_gases/co2/in-situ/surface/

    Args:
        `full`: if `true` Returns full original dataset from 1973 to 2020 in a DataFrame, otherwise returns the subset used in "STL: A Seasonal-Trend Decomposition Procedure Based on Loess" from Cleveland et. al.
    Returns:
        Dataframe or TimeArray containing the descrived dataset.

# Examples
```julia-repl
julia> co2()
[ Info: Dataset used in Cleveland et al. paper
4609×1 TimeArray{Union{Missing, Float64},1,Date,Array{Union{Missing, Float64},1}} 1974-05-17 to 1986-12-31
│            │ co2     │
├────────────┼─────────┤
│ 1974-05-17 │ 333.38  │
│ 1974-05-18 │ 333.11  │
│ 1974-05-19 │ 333.46  │
   ⋮
```
"""
function co2(full = false)

    data = "data/co2.csv"
    path = joinpath(splitdir(@__DIR__)[1], data)
    co2_df = CSV.read(path, DataFrame)

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
