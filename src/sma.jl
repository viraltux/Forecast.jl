function sma(x, n; extend = false, normx = false, normy = false)

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
    
    if extend
        b = n-a
        #TODO optimize by using just the initial and final values
        aux = hcat(res,x)
        aux = mapslices(x->coalesce(x...),aux,dims=2)[:,1]
        sma2 = collect(skipmissing(sma(aux,n)))
        res[1:a+ivp-1]=sma2[1:a+ivp-1]
        N2 = length(sma2)
        res[(N-b+1):N]=sma2[(N2-b+1):N2]
    end

    if normy
        Nm = N-count(x->ismissing(x),x)
        Nmr = N-count(x->ismissing(x),res)
        d = sum(skipmissing(x))/Nm-sum(skipmissing(res))/Nmr
        res = res .+ d
    end

    res
end
