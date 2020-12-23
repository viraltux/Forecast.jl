function co2(full=false)

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
