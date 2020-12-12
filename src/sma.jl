function sma(x, n; center = false, extend = false, norm = false)

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

    if extend
        # extending the x axis
        norm01(xv) = (xv .- minimum(xv)) / maximum( (xv .- minimum(xv)))
        normab(xv,a,b) = norm01(xv) * (b-a) .+ a
        ## initial and final value positions
        ivpr = findfirst(!ismissing, res)
        fvpr = findlast(!ismissing, res)

        x1 = collect(normab(ivpr:fvpr,1,N))
        y1 = collect(skipmissing(res))
        yx = zeros(N)
        
        # linear estimation of the y axis
        for i in eachindex(x1[1:(end-1)])
            xi = x1[i];    yi = y1[i]
            xi1 = x1[i+1]; yi1 = y1[i+1]
            for xv in ceil(xi):1:floor(xi1)
                yx[Int(xv)] = (yi1-yi)*(xv-xi)/(xi1-xi) + yi
            end
        end
        res = yx
    end

    if norm
        
    end
    
    res

end
