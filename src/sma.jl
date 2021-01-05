"""
Package: Forecast

    sma(x, n; center = true)

Smooth a vector of data using a simple moving average.

    Args:
        `x`: Vector of data.
        `n`: Size of the moving average.
        `center`: centers the moving averaged values in the response.
    Returns:
        Vector of moving average smoothed values containing `missing` values to preserve the size of the original vector.

# Examples
```julia-repl
julia> sma(1:5,3;center=true)
5-element Array{Any,1}:
  missing
 2.0
 3.0
 4.0
  missing
```
"""
function sma(x, n; center = true)

    N = length(x)
    @assert 1 <= n & n <= N

    if n == 1
        return x
    end

    res = Vector{Any}(missing,N)

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
