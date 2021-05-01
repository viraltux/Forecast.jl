"""
Package: Forecast

    arsim(Φ,n)
    arsim(Φ,Φ0,n)
    arsim(Φ,Φ0,x0,n)
    arsim(Φ,Φ0,x0,Σ,n)
    arsim(Φ,Φ0,x0,E,n)
    arsim(AR,n)

Simulate a multivariate autoregressive series model.

The simulated series follows the model:

```math
Xt = \\Phi_0 + \\sum_{i=1}^p \\Phi_i \\cdot X_{t-i} + E
```

# Arguments
- `Φ`:            Array with dimensions (m,m,p) for the parameters in the AR model.
- `Φ0`:           Vector size `m` for the constant in the AR model. Default value is 0.
- `x0`:           Array with dimensions (m,p) for the initial value in the AR model. Default value is a random value from zero to one.
- `Σ`:            Variance/Covariance matrix for the AR model with a MvNormal distribution for the noise. Default value is an identity Matrix.
- `n`:            Number of simulations.
- `E`:            Distribution for the error.
- `AR`:           AR struct coming from an `ar` fitting.

# Returns
A multivariate series simulating an AR model each column containing a dimension and ordered by time ascending rows.

# Examples
```julia-repl
julia> arsim(1,10)
10-element Vector{Float64,1}:
[...]
"""
function arsim(Φ,n::Integer)
    arsim(Φ,nothing,nothing,MvNormal(0,0),n)
end

function arsim(Φ,Φ0,n::Integer)
    arsim(Φ,Φ0,nothing,MvNormal(0,0),n)
end

function arsim(Φ,Φ0,x0,n::Integer)
    arsim(Φ,Φ0,x0,MvNormal(0,0),n)
end

function arsim(Φ,Φ0,x0,Σ::Number,n::Integer)
    arsim(Φ,Φ0,x0,MvNormal(1,Σ),n)
end

function arsim(Φ,Φ0,x0,Σ::Vector,n::Integer)
    arsim(Φ,Φ0,x0,MvNormal(collect(Diagonal(Σ))),n)
end

function arsim(Φ,Φ0,x0,Σ::Matrix,n::Integer)
    arsim(Φ,Φ0,x0,MvNormal(Σ),n)
end

function arsim(ar_str::AR,n::Integer)
    Φ = ar_str.Φ
    Φ0 = ar_str.Φ0
    Σ = ar_str.Σ
    arsim(Φ,Φ0,nothing,Σ,n)
end

function arsim(Φ,Φ0,x0,E::Distribution,n::Integer)
    Φ = compact(Φ)

    m,p = arsize(Φ)
    
    Φ0 = isnothing(Φ0) ? compact(zeros(m)) : compact(Φ0)
    x0 = isnothing(x0) ? compact(rand(m*p)) : compact(x0)
    E = size(E)[1] == 0 ? MvNormal(m,1) : E
    
    e = rand(E,n)
    x = compact(Array{Float64}(undef,n,m))
    
    # reshape and format for matrix operation
    Φ = Φ  isa Number ? reshape([Φ],1,1) : reshape(Φ,m,m*p)
    Φ0 = Φ0 isa Number ? [Φ0] : Φ0
    x0 = x0 isa Number ? [x0] : reshape(x0,m*p,1)

    for i in 1:n
        x[i,:] = Φ * x0 + Φ0  + e[:,i]
        x0 = vcat(x[i,:],x0)[1:m*p]
    end

    compact(x)
end

function arsize(Φ)
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
        @error "Φ should have less the 4 dimensions"
    end
    (m,p)
end
