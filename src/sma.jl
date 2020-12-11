function sma(x, n; normx = false, normy = false)

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
    a = normx ? 0 : div(n,2)

    # initial moving average value
    ma = sum(x[ivp:n+ivp])/n
    res[a+ivp] = ma
    for i in 1:(N-n-ivp-(N-fvp))
        resai = ma + (x[n+ivp+i] - x[ivp+i])/n
        # missing values are imputed
        res[a+ivp+i] = ismissing(resai) ? ma : resai
        ma = res[a+ivp+i]
    end
    
    if normy
        Nm = N-count(x->ismissing(x),x)
        Nmr = N-count(x->ismissing(x),res)
        d = sum(skipmissing(x))/Nm-sum(skipmissing(res))/Nmr
        res = res .+ d
    end

    res
end


