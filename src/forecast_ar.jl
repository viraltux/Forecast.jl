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
function forecast(xar::AR, n::Integer; xlevels = (0.8,.95))

    Φ,Φ0,Σ = xar.Φ,xar.Φ0,xar.Σ
    nd = max(ndims(Σ),1)
    x0 = compact(xar.x[end-size(Σ,1)[1]+1:end,:])
    E = MvNormal(nd,0)
    mu = arsim(Φ,Φ0,x0,E,n)
    se = sqrt.(fvar(Φ,Σ,n))

    # Prediction Intervals
    a1 = xlevels[1]
    a2 = xlevels[2]

    z1 = mapslices(x -> Distributions.quantile.(Normal.(0,x),
                                               a1 + (1-a1)/2),
                                               reshape(se,n,nd), dims = 2)
    z1 = compact(z1)
    lower1 = mu .- z1
    upper1 = mu .+ z1
    z2 = mapslices(x -> Distributions.quantile.(Normal.(0,x),
                                               a2 + (1-a2)/2),
                                               reshape(se,n,nd), dims = 2)
    z2 = compact(z2)
    lower2 = mu .- z2
    upper2 = mu .+ z2
    
   FORECAST(xar, xlevels,
            mu, hcat(upper1,upper2), hcat(lower1,lower2),
            "")

end

"""
Forecast recurrence variance for a linear combination
"""
function fvar(Φ,Σ,n)
    Φ = compact(Φ)

    d = ndims(Φ)
    sΦ = size(Φ)
    if d == 0
        m,p = 1,1
    elseif d == 1
        m,p = 1,sΦ[1]
    elseif d == 2
        m,p = sΦ[1],1
    elseif d == 3
        m,p = sΦ[1],sΦ[3]
    else
        @error "Φ should have less than 4 dimensions"
    end

    v = zeros(n,m)
    v[1,:] = m <= 1 ? [Σ^2] : diag(Σ^2)
    # reshape and format for matrix operation
    Φ = d == 0 ?  reshape([Φ],1,1) : reshape(Φ,m,m*p)
    
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
