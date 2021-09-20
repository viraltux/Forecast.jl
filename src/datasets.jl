"""
Package: Forecast

    co2(full = false)

Return dataset with atmospheric Carbon Dioxide Dry Air Mole Fractions from quasi-continuous measurements at Mauna Loa, Hawaii.

K.W. Thoning, A.M. Crotwell, and J.W. Mund (2020), Atmospheric Carbon Dioxide Dry Air Mole Fractions from continuous measurements at Mauna Loa, Hawaii, Barrow, Alaska, American Samoa and South Pole. 1973-2019, Version 2020-08 National Oceanic and Atmospheric Administration (NOAA), Global Monitoring Laboratory (GML), Boulder, Colorado, USA https://doi.org/10.15138/yaf1-bk21 FTP path: ftp://aftp.cmdl.noaa.gov/data/greenhouse_gases/co2/in-situ/surface/

# Arguments
- `full`: if `true` Returns the full original dataset from 1973 to 2020 in a DataFrame, otherwise returns the subset used in "STL: A Seasonal-Trend Decomposition Procedure Based on Loess" from Cleveland et. al. Its default value is `false`.

# Returns
Dataframe containing the descrived dataset.

# Examples
```julia-repl
julia> co2()
[ Info: Dataset used in Cleveland et al. paper
4612×2 DataFrame
  Row │ date        co2        
      │ Date        Float64?   
──────┼────────────────────────
    1 │ 1974-05-17      333.38
    2 │ 1974-05-18      333.11
    3 │ 1974-05-19      333.46
   [...]
```
"""
function co2(full::Bool = false)

    data = "data"
    path = joinpath(splitdir(@__DIR__)[1], data)
    
    co2_df = open(joinpath(path, "co2.csv.gz")) do file
        CSV.read(GzipDecompressorStream(file),DataFrame)
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
    co2_sdf = DataFrame([dates_co2, co2_df.value],[:date,:co2])

    # paper dates
    # using 17 may instead 17 april (as the paper claims)
    # since there seems to be no data for april in 1974
    dates_co2_stl = Dates.Date(1974,5,17):Dates.Day(1):Dates.Date(1986,12,31)

    # revome leap year day
    dates_co2_stl = filter(dates_co2_stl) do x
        (Dates.isleapyear(x) & (Dates.month(x) == 2)) ? Dates.day(x) != 29 : true
    end

    @where(co2_sdf, in.(:date,  Ref(dates_co2_stl)))
    
end


"""
Package: Forecast

    seaborne(full = false)

Return estimates of world seaborne trade from AIS data collected by Marine Traffic.

By default a DataFrame containing deadweight imports for France, Germany and the United Kingdom from 2015-04-01 to 2021-05-02 is returned, otherwise a DataFrame is returned for the same countries with import and exports for the below fields:

num_pc:     number of port calls
mtc:        metric tons of cargo
dwt:        deadweight tonnage
suffix_ma:  30-day moving averages

Data available at UN COMTRADE Monitor.Cerdeiro, Komaromi, Liu and Saeed (2020). 

# Returns
 Dataframe containing the seaborne dataset.

# Examples
```julia-repl
julia> seaborne()
[ Info: Seaborne deadweight trade imports from AIS
2199×4 DataFrame
  Row │ Date        France   Germany  UK     
      │ Date        Int64    Int64    Int64  
──────┼──────────────────────────────────────
    1 │ 2015-04-01   507946   878377  599573
    2 │ 2015-04-02   332043  1501614  772714
    3 │ 2015-04-03   810077   941663  262994
   [...]
```
"""
function seaborne(full::Bool = false)

    data = "data"
    path = joinpath(splitdir(@__DIR__)[1], data)
    
    sb_df = open(joinpath(path, "seaborne.csv.gz")) do file
        CSV.read(GzipDecompressorStream(file),DataFrame,
                 dateformat = "mm/dd/yyyy HH:MM:SS pp",
                 types = Dict(:year => Date))
    end

    if full
        @info "Full dataset with estimates of world seaborne trade from AIS"
        return sb_df
    else
        @info "Seaborne deadweight trade imports from AIS"
    end

    # Group by Country and Flow to select Imports
    gbt = groupby(sb_df, [:country_name, :flow], sort=true)
    DataFrame(Dict(:Date => gbt[1].date ,
                   :UK => gbt[1].dwt, :Germany => gbt[3].dwt, :France => gbt[5].dwt))
end


"""
Package: Forecast

    quakes()

Return the number of earthquakes per year on earth with a magnitude higher or equal to six from 1950 to 2020. The data has been collected from https://earthquake.usgs.gov/ and aggregated.

# Examples
```julia-repl
julia> quakes()
71×2 DataFrame
 Row │ year        quakes 
     │ Date        Int64  
─────┼────────────────────
   1 │ 1950-01-01     138
   2 │ 1951-01-01     151
   3 │ 1952-01-01     181
   [...]
```
"""
function quakes()

    data = "data"
    path = joinpath(splitdir(@__DIR__)[1], data)
    
    open(joinpath(path, "quakes.csv.gz")) do file
        CSV.read(GzipDecompressorStream(file),DataFrame,
                 dateformat = "yyyy",
                 types = Dict(:year => Date))
    end
    
end

"""
Package: Forecast

    air()

Return the classic Box & Jenkins airline data. Monthly totals of international airline passengers from 1949 to 1960.

Box, G. E. P., Jenkins, G. M. and Reinsel, G. C. (1976) _Time Series Analysis, Forecasting and Control._ Third Edition. Holden-Day. Series G.

# Returns
Dataframe containing the descrived dataset.

# Examples
```julia-repl
julia> air()
71×2 DataFrame
 Row │ year        quakes 
     │ Date        Int64  
─────┼────────────────────
   1 │ 1950-01-01     138
   2 │ 1951-01-01     151
   3 │ 1952-01-01     181
   [...]
```
"""
function air()

    data = "data"
    path = joinpath(splitdir(@__DIR__)[1], data)
    
    df = open(joinpath(path, "air.csv.gz")) do file
        CSV.read(GzipDecompressorStream(file),DataFrame)
    end
    
    t =  DataFrame([Date.(df.year,df.month)],[:Date])
    x = DataFrame([df.passengers],[:Passengers])

    [t x]
end


"""
Package: Forecast

    london()

Return ten years of monthly data about weather and crime in Greater London from 2008 to 2018.

Data has been collected and joined from london.gov.uk and metoffice.gov.uk (Heathrow Station).

# Weather Variables

- `MaxTemp`:  Mean daily maximum temperature in C°
- `MinTemp`:  Mean daily minimum temperature in C°
- `AirFrost`: Days of air frost
- `Rain`:     Total rainfall in mm
- `Sun`:      Total sunshine durationin hours

# Crime Variables and its aggregated categories
```
┌─────────────────────────────┬────────────────────────────────────────┐
│                       Crime │                               Category │
├─────────────────────────────┼────────────────────────────────────────┤
│                    Burglary │            Burglary in Other Buildings │
│                             │                 Burglary in a Dwelling │
│                      Damage │            Criminal Damage To Dwelling │
│                             │       Criminal Damage To Motor Vehicle │
│                             │      Criminal Damage To Other Building │
│                             │                  Other Criminal Damage │
│                       Drugs │                       Drug Trafficking │
│                             │                            Other Drugs │
│                             │                    Possession Of Drugs │
│                       Fraud │                     Counted per Victim │
│                             │                  Other Fraud & Forgery │
│                       Other │                         Going Equipped │
│                             │                       Other Notifiable │
│                     Robbery │                      Business Property │
│                             │                      Personal Property │
│                      Sexual │                           Other Sexual │
│                             │                                   Rape │
│                       Theft │                  Handling Stolen Goods │
│                             │ Motor Vehicle Interference & Tampering │
│                             │                            Other Theft │
│                             │                     Other Theft Person │
│                             │               Theft From Motor Vehicle │
│                             │                       Theft From Shops │
│                             │          Theft/Taking Of Motor Vehicle │
│                             │            Theft/Taking of Pedal Cycle │
│                    Violence │                    Assault with Injury │
│                             │                         Common Assault │
│                             │                   Grievous Bodily Harm │
│                             │                             Harassment │
│                             │                                 Murder │
│                             │                       Offensive Weapon │
│                             │                         Other violence │
│                             │                           Wounding/GBH │
└─────────────────────────────┴────────────────────────────────────────┘
```

# Returns
Dataframe containing the descrived dataset.

# Examples
```julia-repl
julia> london()
132×15 DataFrame
 Row │ Date        MaxTemp  MinTemp  AirFrost  
     │ Date        Float64  Float64  Int64     
─────┼─────────────────────────────────────────[...]
   1 │ 2008-01-01     10.4      4.7         0  
   2 │ 2008-02-01     11.0      2.0         7  
   3 │ 2008-03-01     10.6      3.7         2
   [...]
```
"""
function london()

    data = "data"
    path = joinpath(splitdir(@__DIR__)[1], data)
    
    df = open(joinpath(path, "london.csv.gz")) do file
        CSV.read(GzipDecompressorStream(file),DataFrame,
                 dateformat = "yyyy-mm-dd",
                 types = Dict(:Date => Date))
    end
    
end
