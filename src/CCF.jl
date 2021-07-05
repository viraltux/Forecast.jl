"""
Package: Forecast

Store results from the functions `acf`, `ccf` and `pacf`

# Arguments
`ccf`    An array with results from ccf, acf and pacf
`N`      Length of ccf
`type    Type of CCF
`lag`    Maximum number of lags
`alph`   CI thresholds
`ci`     CI for the alpha
`auto`   Auto-correlation
`call`   Method called to generate ccf
"""
mutable struct CCF
    ccf::AbstractArray{<:Real}
    N::Integer
    type::String
    lag::Integer
    alpha::Tuple{Real,Real}
    ci::Tuple{Real,Real}
    auto::Bool
    call::String
end
