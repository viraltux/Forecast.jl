function sma(x, n; extended = false)
    #TODO implement extended flag for doing ma as long as x

    if n == 1
        return x
    end

    N = length(x)
    @assert 1 <= n & n <= N
    res = Vector{Any}(missing,N)

    # initial and final value positions
    ivp = findfirst(!ismissing, x)
    fvp = findlast(!ismissing, x)
    
    # using missing values in 1:a to center ma 
    a = div(n-1,2)

    # initial moving average value
    ma = sum(x[ivp:n+ivp])/n
    res[a+ivp+1] = ma
    for i in 2:(N-n+1-ivp-(N-fvp))
        resai = ma + (x[n+i-1+ivp] - x[i-1+ivp])/n
        # missing values are imputed
        res[a+i+ivp] = ismissing(resai) ? ma : resai
        ma = res[a+i+ivp]
    end

    if extended
        b = n-a-1
        aux = hcat(res,x)
        aux = mapslices(x->coalesce(x...),aux,dims=2)[:,1]
        sma2 = collect(skipmissing(sma(aux,n)))
        res[1:a+ivp]=sma2[1:a+ivp]
        N2 = length(sma2)
        res[(N-b+1):N]=sma2[(N2-b+1):N2]
    end

    res
end
