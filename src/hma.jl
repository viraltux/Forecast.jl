# Ported https://github.com/bpalmer4 Python implementation 
# Reference:
# - https://www.abs.gov.au/websitedbs/d3310114.nsf/4a256353001af3ed4b2562bb00121564/5fc845406def2c3dca256ce100188f8e!OpenDocument

# Henderson symmetric weights
function hmaSymmetricWeights(n::Integer)

    m = (n-1)÷2
    m1 = (m+1)^2
    m2 = (m+2)^2
    m3 = (m+3)^2
    d = 8*(m+2)*(m2-1)*(4*m2-1)*(4*m2-9)*(4*m2-25)

    w = map(x -> round(315*(m1-x^2)*(m2-x^2)*(m3-x^2)*(3*m2-11*x^2-16)/d, digits=3), 0:m+1)
    u = vcat(w[end-1:-1:1],w[2:end-1])
    
    return mod(n, 2) != 0 ? u : vcat(u, missing)

end

function hmaAsymmetricWeights(m::Integer, w::AbstractArray{<:Real})
    n = length(w)

    @assert m <= n "The m argument must be less than w"
    @assert m >= (n-1)÷2 "The m argument must be greater than (n-1)/2"

    sumResidual = sum(w[range(m + 1, n, step = 1)]) / m
    sumEnd =  sum(map(x -> (x-(m+1.0)/2.0)*w[x], m+1:n))

    ic = n < 13 ? 1.0 : (13 <= n < 15 ? 3.5 : 4.5)
    
    b2s2 = 4/pi/ic^2
    denominator = 1.0 + ((m*(m-1.0)*(m+1.0) / 12.0 ) * b2s2)
    map(r -> round(w[r] + sumResidual + (((r-1+1.0) - (m+1.0)/2.0) * b2s2) / denominator * sumEnd, digits = 3), 1:m)
end


"""
Package: Forecast

    hma(s, n)

Applies the Henderson moving average filter to dataset `s` with `n`-term.


"Henderson moving averages are filters which were derived by Robert Henderson in 1916 for use in actuarial applications. They are trend filters, commonly used in time series analysis to smooth seasonally adjusted estimates in order to generate a trend estimate. They are used in preference to simpler moving averages because they can reproduce polynomials of up to degree 3, thereby capturing trend turning points.

The ABS uses Henderson moving averages to produce trend estimates from a seasonally adjusted series. The trend estimates published by the ABS are typically derived using a 13 term Henderson filter for monthly series, and a 7 term Henderson filter for quarterly series.

Henderson filters can be either symmetric or asymmetric. Symmetric moving averages can be applied at points which are sufficiently far away from the ends of a time series. In this case, the smoothed value for a given point in the time series is calculated from an equal number of values on either side of the data point." - Australian Bureau of Statistics (www.abs.gov.au)

# Arguments
- `s`: Observations' support.
- `n`: Observation values. Tipically 13 or 7 for quarterly data, larger values may need a BigInt type.

# Returns
An array of Henderson filter smoothed values provided in `s`.

# Examples
```julia-repl
julia> hma(rand(1000), BigInt(303)))
1000-element Vector{BigFloat}:
[...]
```
"""
function hma(s::AbstractArray{<:Real}, n::Integer)

    @assert isodd(n) "n must be odd"
    @assert n >= 5 "n must be >= 5"
    @assert length(s) >= n "dataset must be >= than n"

    w = hmaSymmetricWeights(n)
    m = (n-1) ÷ 2
    ls = length(s)

    function hmai(i)
        if i - 1 < m
            u = hmaAsymmetricWeights(m + i, w)[end:-1:1]
            sum(s[1:i + m] .* u)
        elseif i - 1 + m >= ls
            u = hmaAsymmetricWeights(m + ls - i + 1, w)
            sum(s[i-m:ls] .* u)
        else
            sum(s[i-m:i+m] .* w)
        end
    end
    
    map(x -> hmai(x), 1:ls)
end
