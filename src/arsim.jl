"""
Package: Forecast

    arsim(Φ,n)
    arsim(Φ,Φ0,n)
    arsim(Φ,Φ0,x0,n)
    arsim(Φ,Φ0,x0,Σ,n)
    arsim(Φ,Φ0,x0,E,n)
    arsim(AR,n::Integer)

Simulate a multivariate autoregressive series model.

The simulated series follows the model:

```math
Xt = \\Phi_0 + \\sum_{i=1}^p \\Phi_i \\cdot X_{t-i} + E
```

# Arguments
`Φ`            Array with dimensions (m,m,p) for the parameters in the AR model.
`Φ0`           Vector size `m` for the constant in the AR model. Default value is 0.
`x0`           Array with dimensions (m,p) for the initial value in the AR model. Default value is a random value from zero to one.
`Σ`            Variance/Covariance matrix for the AR model with a MvNormal distribution for the noise. Default value is an identity Matrix.
`n`            Number of simulations.
x`E`            Distribution for the error.
`AR`           AR struct coming from an `ar` fitting.

# Returns
A multivariate series simulating an AR model each column containing a dimension and ordered by time ascending rows.
"""
function arsim(Φ,n::Integer)
    Φ = Φ isa Number ? [Φ] : Φ
    sΦ = size(Φ)
    m = Φ isa Vector ? 1 : sΦ[1]
    p = length(sΦ) > 2 ? sΦ[3] : 1
    p = Φ isa Vector ? sΦ[1] : p

    Φ0 = zeros(m)
    x0 = rand(m,p)
    arsim(Φ,Φ0,x0,MvNormal(m,1),n::Integer)
end

function arsim(ar_str::AR,n::Integer)
    Φ = ar_str.Φ
    sΦ = size(Φ)
    m = sΦ[1]
    p = length(sΦ) > 2 ? sΦ[3] : 1

    Φ0 = ar_str.Φ0
    x0 = ar_str.x0
    Σ = ar_str.Σ
    arsim(Φ,Φ0,x0,Σ,n::Integer)
end

function arsim(Φ,Φ0,n::Integer)
    sΦ = size(Φ)
    m = sΦ[1]
    p = length(sΦ) > 2 ? sΦ[3] : 1

    x0 = rand(m,p)
    arsim(Φ,Φ0,x0,MvNormal(m,1),n::Integer)
end

function arsim(Φ,Φ0,x0,n::Integer)
    arsim(Φ,Φ0,x0,MvNormal(size(Φ)[1],1),n::Integer)
end

function arsim(Φ,Φ0,x0,Σ::AbstractArray,n::Integer)
    arsim(Φ,Φ0,x0,MvNormal(Σ),n::Integer)
end

function arsim(Φ,Φ0,x0,E::Distribution,n::Integer)

    Φ = Φ isa Number ? [Φ] : Φ
    sΦ = size(Φ)
    m = Φ isa Vector ? 1 : sΦ[1]
    p = length(sΦ) > 2 ? sΦ[3] : 1
    p = Φ isa Vector ? sΦ[1] : p

    # @assert m == size(Φ0)[1]
    # @assert size(x0) = (p,m)

    e = rand(E,n::Integer)
    x = Array{Float64}(undef,n,m)

    # reshape for matrix operation
    Φ = reshape(Φ,m,m*p) #TODO check size Φ if just two dimensions not reshape
    x0 = reshape(x0,m*p,1)
    
    for i in 1:n
        x[i,:] =  Φ * x0 + Φ0  + e[:,i]
        x0 = vcat(x[i,:],x0)[1:m*p]
    end
    x
end
