"""
Package: Forecast

Store results from the functions acf, ccf and pacf

# Arguments

    `ccf::T`          An array with results from ccf, acf and pacf
    `N::Integer`      Length of ccf
    `type::String`    Type of CCF
    `lag::Integer`    Maximum number of lags
    `alpha::Tuple`    CI thresholds
    `ci::Tuple`       CI for alpha
    `auto::Bool`      Auto-correlation
    `call::String`    Method called to generate ccf
"""
mutable struct CCF{T<:AbstractArray}
    ccf::T             # An array with results from ccf, acf and pacf
    N::Integer         # Length of ccf
    type::String       # Type of CCF
    lag::Integer       # Maximum number of lags
    alpha::Tuple       # CI thresholds
    ci::Tuple          # CI for alpha
    auto::Bool         # Auto-correlation
    call::String       # Method called to generate ccf
end
