"""
Package: Forecast

    predict(AR,n; ci = (0.8,.95))

Preedict a multivariate autoregressive series model.
    
The predicted series follows the model:

```math
Xt = \\Phi_0 + \\sum_{i=1}^p \\Phi_i \\cdot X_{t-i} + E
```

# Arguments
`xar`           AR struct coming from an `ar` fitting.
`n`             Number of elements to be predicted
`level`            Confidence Interval for the prediction, its default value is (0.8, 0.95)

# Returns
A multivariate series simulating an AR model each column containing a dimension and ordered by time ascending rows.
"""
function forecast(xar::AR, n::Integer; levels = (0.8,.95))

    @assert n > 0 "n must be greater than 0"
    
    Φ,Φ0,Σ = xar.Φ,xar.Φ0,xar.Σ

    m,p = arsize(Φ)

    x = xar.x isa TimeArray ? values(xar.x) : xar.x
    x0 = compact(x[end:-1:end-p+1,:]')
    E = MvNormal(m,0)
    mu = arsim(Φ,Φ0,x0,E,n)
    se = sqrt.(fvar(Φ,Σ,n))

    # Prediction Intervals
    a1 = levels[1]
    a2 = levels[2]

    z1 = mapslices(x -> Distributions.quantile.(Normal.(0,x),
                                               a1 + (1-a1)/2),
                                               reshape(se,n,m), dims = 2)
    z1 = compact(z1)
    lower1 = mu .- z1
    upper1 = mu .+ z1
    z2 = mapslices(x -> Distributions.quantile.(Normal.(0,x),
                                               a2 + (1-a2)/2),
                                               reshape(se,n,m), dims = 2)
    z2 = compact(z2)
    lower2 = mu .- z2
    upper2 = mu .+ z2
    
   FORECAST(xar, levels,
            mu, hcat(upper1,upper2), hcat(lower1,lower2),
            "Forecasting "*xar.call)

end

"""
Forecast recurrence variance for a linear combination
"""
function fvar(Φ,Σ,n)
    Φ = compact(Φ)

    m,p = arsize(Φ)

    v = zeros(n,m)
    v[1,:] = m <= 1 ? [Σ^2] : diag(Σ^2)
    # reshape and format for matrix operation
    Φ = Φ isa Number ?  reshape([Φ],1,1) : reshape(Φ,m,m*p)
    
    for i in 2:n
        for j in 1:m
            if i <= p
                v[i,j] = (Φ[j,:].^2)' * repeat(v[1:p,j],inner=m) +  v[1,j]
            else
                v[i,j] = (Φ[j,:].^2)' * repeat(v[i-p:i-1,j],inner=m) +  v[1,j]
            end
        end
    end
    
    compact(v)
end
