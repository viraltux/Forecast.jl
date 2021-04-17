"""
Package: Forecast

    arsim(Φ,n)
    arsim(Φ,Σ,n)
    arsim(Φ,Φ0,x0,n)
    arsim(Φ,Φ0,x0,n,E)

Simulate a multivariate autoregressive series model.
    
The simulated series follows the model:

```math
Xt = \\Phi_0 + \\sum_{i=1}^p \\Phi_i \\cdot X_{t-i} + E
```

# Arguments
`Φ`
`Φ0`
`Σ`
`x0` 
`n`
`E`            Distribution for the error, if not set a default Gaussian one is used instead.


# Returns
A multivariate series simulating an AR model each column containing a dimension and ordered by time ascending rows.
"""
function arsim(Φ,n)
    sΦ = size(Φ)
    m = sΦ[1]
    p = length(sΦ) > 2 ? sΦ[3] : 1

    Φ0 = zeros(m)
    x0 = rand(m,p)
    arsim(Φ,Φ0,x0,n)
end

function arsim(Φ,Σ,n)
    sΦ = size(Φ)
    m = sΦ[1]
    p = length(sΦ) > 2 ? sΦ[3] : 1

    Φ0 = zeros(m)
    x0 = rand(m,p)
    arsim(Φ,Φ0,x0,n,MvNormal(Σ))
end

function arsim(Φ,Φ0,x0,n)
    arsim(Φ,Φ0,x0,n,MvNormal(size(Φ)[1],1))
end

function arsim(Φ,Φ0,x0,n,E)

    sΦ = size(Φ)
    m = sΦ[1]
    p = length(sΦ) > 2 ? sΦ[3] : 1
    # @assert m == size(Φ0)[1]
    # @assert size(x0) = (p,m)

    e = rand(E,n)
    x = Array{Float64}(undef,n,m)

    # reshape for matrix operation
    Φ = reshape(Φ,d,m*p) #TODO check size Φ if just two dimensions not reshape
    x0 = reshape(x0,m*p,1)
    
    for i in 1:n
        x[i,:] =  Φ * x0 + Φ0  + e[:,i]
        x0 = vcat(x[i,:],x0)[1:m*p]
    end
    x
end
