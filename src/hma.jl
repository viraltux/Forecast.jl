"""
Implementation of Henderson filter in Julia.
"""
# Ported https://github.com/bpalmer4 Python implementation 
# Reference:
# - https://www.abs.gov.au/websitedbs/d3310114.nsf/4a256353001af3ed4b2562bb00121564/5fc845406def2c3dca256ce100188f8e!OpenDocument


function hmaSymmetricWeights(n::Int)
    m = (n-1)÷2
    m1 = (m+1)^2
    m2 = (m+2)^2
    m3 = (m+3)^2
    d = BigInt(8*(m+2)*(m2-1)*(4*m2-1)*(4*m2-9)*(4*m2-25))

    w = map(x -> round(Float16(BigFloat(315*(m1-(x^2))*(m2-(x^2))*(m3-(x^2))*(3*m2-11*(x^2)-16))/d), digits = 3), range(0, m+1, step = 1))[end-1:-1:1]
    u = vcat(w, w[end-1:-1:1])
    return if mod(n, 2) != 0
        u
    else 
        vcat(u, missing)
    end
end

function hmaAsymmetricWeights(m::Int, w::AbstractArray{<:Number})
    n = length(w)

    @assert m <= n "The m argument must be less than w"
    @assert m >= (n-1)÷2 "The m argument must be greater than (n-1)/2"

    sumResidual = sum(w[range(m + 1, n, step = 1)]) / m
    sumEnd = sum(map(x -> (x - ((m+1.0)/2.0)) * w[x], range(m+1, n, step = 1)))

    ic = if n >= 13 && n < 15
        3.5
    elseif n >= 15
        4.5
    else 1.0
    end
    b2s2 = BigFloat(4/pi/ic^2)
    denominator = 1.0 + ((m*(m-1.0)*(m+1.0) / 12.0 ) * b2s2)

    u = map(r -> round(Float16(w[r] + sumResidual + (((r-1+1.0) - (m+1.0)/2.0) * b2s2) / denominator * sumEnd), digits = 3), 
        range(1, m, step = 1))
    return (u)
end

"""
Package: Forecast

    hma(s, n)

Applies the Henderson moving average filter to dataset `s` with `n`-term.

Information about involved processes and application can be found at the following:

("Time Series Analysis: The Process of Seasonal Adjustment")[https://www.abs.gov.au/websitedbs/d3310114.nsf/4a256353001af3ed4b2562bb00121564/5fc845406def2c3dca256ce100188f8e!OpenDocument]
Australian Bureau of Statistics

# Arguments
- `s`: Observations' support.
- `n`: Observation values.

# Returns
An array of Henderson filter smoothed values provided in `s`.

# Examples
```julia-repl
julia> hma(rand(24), 13)
24-element Array{Float64,1}:
[...]
```
"""
function hma(s::AbstractArray{<:Number}, n::Int)
    @assert isodd(n) "n must be odd"
    @assert n >= 5 "n must be >= 5"
    @assert length(s) >= n "dataset must be >= than n"

    w = hmaSymmetricWeights(n)
    m = (n-1) ÷ 2
    l = length(s)
    map(i -> if i - 1 < m
            u = hmaAsymmetricWeights(m + i, w)[end:-1:1]
            sum(s[1:i + m] .* u)
        elseif i - 1 + m >= l
            u = hmaAsymmetricWeights(m + l - i + 1, w)
            sum(s[i-m:l] .* u)
        else
            sum(s[i-m:i+m] .* w)
        end, range(1, l, step = 1))
end