"""
Package: Forecast

    arsim(Φ,n)
    arsim(Φ,Σ,n)
    arsim(Φ,Φ0,x0,n)
    arsim(Φ,Φ0,x0,n,E)
    predict(AR,n)

Preedict a multivariate autoregressive series model.
    
The predicted series follows the model:

```math
Xt = \\Phi_0 + \\sum_{i=1}^p \\Phi_i \\cdot X_{t-i} + E
```

# Arguments
`AR`           AR struct coming from an `ar` fitting.
`n`            Number of elements to be predicted
`ci`           Confidence Interval for the prediction, its default value is (0.8, 0.95)

# Returns
A multivariate series simulating an AR model each column containing a dimension and ordered by time ascending rows.
"""
function predict(arstr::AR, n::Integer, ci = (0.8,.95))
    #TODO implement ci
    
    sΦ = size(arstr.Φ)
    m = sΦ[1]
    p = length(sΦ) > 2 ? sΦ[3] : 1

    e = rand(MvNormal(arstr.Σ),n)
    x = Array{Float64}(undef,n,m)

    # reshape for matrix operation
    Φ = reshape(arstr.Φ,m,m*p) #TODO check size Φ if just two dimensions not reshape
    x0 = reshape(arstr.x0,m*p,1)
    
    for i in 1:n
        x[i,:] =  Φ * x0 + Φ0  + e[:,i]
        x0 = vcat(x[i,:],x0)[1:m*p]
    end
    x
end
