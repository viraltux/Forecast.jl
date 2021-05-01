"""
Package: Forecast

Store results from the functions `acf`, `ccf` and `pacf`

# Arguments
`ccf::AbstractArray`    An array with results from ccf, acf and pacf
`N::Integer`            Length of ccf
`type::String`          Type of CCF
`lag::Integer`          Maximum number of lags
`levels::Tuple`         CI thresholds
`ci::Tuple`             CI for the levels
`auto::Bool`            Auto-correlation
`call::String`          Method called to generate ccf
"""
mutable struct CCF
    ccf::AbstractArray
    N::Integer
    type::String
    lag::Integer
    levels::Tuple{Float64,Float64}
    ci::Tuple{Float64,Float64}
    auto::Bool
    call::String
end
