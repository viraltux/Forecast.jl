"""
Package: Forecast

    sma(x, n)
    sma(x, n, center)

Smooth a vector of data using a simple moving average.

# Arguments
- `x`: Vector of data.
- `n`: Size of the moving average.
- `center`: if true centers the moving averaged values in the response padding with `missing` values, otherwise the padding takes place at the end.

# Returns
Vector of moving average smoothed values containing `missing` values to preserve the size of the original vector.

# Examples
```julia-repl
julia> sma(1:5,3,true)
5-element Array{Any,1}:
  missing
 2.0
 3.0
 4.0
  missing
```
"""
function sma(x, n)
    n == 1 && return x
    x = convert(Vector{Float64},x)
    n = convert(Int64,n)
    sma(x,n)
end

function sma(x::Vector{Float64}, n::Int64)::Vector{Float64}

    n == 1 && return x
    
    N = length(x)
    @assert 1 <= n & n <= N

    res = Vector{Float64}(undef, N-n+1)
    
    # initial moving average value
    res[1] = ma = sum(x[1:n])/n
    for i in 1:N-n
        res[1+i] = ma = ma + (x[n+i] - x[i])/n
    end

    res

end


function sma(x, n, center::Bool)
    n == 1 && return x
    x = convert(Vector{Float64},x)
    n = convert(Int64,n)
    sma(x,n,center)
end


function sma(x::Vector{<:Union{Missing,Real}},
             n::Integer,
             center::Bool)::Vector{Union{Missing,Float64}}

    n == 1 && return x

    N = length(x)
    @assert 1 <= n & n <= N

    res = Vector{Union{Missing,Real}}(missing,N)

    # initial and final value positions
    ivp = findfirst(!ismissing, x)
    fvp = findlast(!ismissing, x)
    
    # using missing values to center ma 
    a = center ? div(n,2) : 0

    # initial moving average value
    ma = sum(x[ivp:n+ivp-1])/n
    res[a+ivp] = ma
    for i in 1:(N-n-ivp-(N-fvp)+1)
        resai = ma + (x[n+ivp+i-1] - x[ivp+i-1])/n
        # missing values are imputed
        res[a+ivp+i] = ismissing(resai) ? ma : resai
        ma = res[a+ivp+i]
    end

    res

end

