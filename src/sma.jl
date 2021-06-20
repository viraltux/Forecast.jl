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
function sma(x::AbstractVector{T}, n::Integer) where T<:Number

    n == 1 && return x
    N = length(x)
    @assert 1 <= n <= N
    
    V = Base.promote_op(/, T, typeof(n))
    res = Vector{V}(undef, N-n+1)
    
    # initial moving average value
    res[1] = ma = sum(x[1:n])/n
    for i in 1:N-n
        res[1+i] = ma = ma + (x[n+i] - x[i]) / n
    end

    res

end

function sma(x::AbstractVector{<:Union{Missing,T}},
             n::Integer,
             center::Bool) where T<:Number

    n == 1 && return x
    N = length(x)
    @assert 1 <= n <= N

    res = sma(collect(skipmissing(x)),n)

    # Missing padding
    ivp = repeat([missing],  findfirst(!ismissing, x)-1 + n÷2)
    fvp = repeat([missing], N-findlast(!ismissing, x)-1 + n-n÷2)

    center ? vcat(ivp,res,fvp) : vcat(res,ivp,fvp)
        
end

