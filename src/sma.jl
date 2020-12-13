function sma(x, n; center = false)

    if n == 1
        return x
    end

    N = length(x)
    @assert 1 <= n & n <= N
    res = Vector{Any}(missing,N)

    # initial and final value positions
    ivp = findfirst(!ismissing, x)
    fvp = findlast(!ismissing, x)
    
    # using missing values to center ma 
    a = center ? div(n,2) : 0

    # initial moving average value
    ma = sum(x[ivp:n+ivp])/n
    res[a+ivp] = ma
    for i in 1:(N-n-ivp-(N-fvp))
        resai = ma + (x[n+ivp+i] - x[ivp+i])/n
        # missing values are imputed
        res[a+ivp+i] = ismissing(resai) ? ma : resai
        ma = res[a+ivp+i]
    end

    res

end
