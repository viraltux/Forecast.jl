"""
Package: Forecast

    arsim(Φ,Φ0,x0,n; Σ,E,fix)
    arsim(AR,n;fix)

Simulate a multivariate autoregressive series model.

The simulated series follows the model:

```math
Xt = \\Phi_0 + \\sum_{i=1}^p \\Phi_i \\cdot X_{t-i} + E
```

# Arguments
- `AR`:           AR struct coming from an `ar` model.
- `Φ`:            Array with dimensions (m,m,p) for the parameters in the AR model.
- `Φ0`:           Vector size `m` for the constant in the AR model. Default value is 0.
- `x0`:           Array with dimensions (m,p) for the initial value in the AR model. Default value is a random value from zero to one.
- `n`:            Number of simulations.

- `fix`:          Matrix{Union{Missing,Float64}} containing values to be fixed in th simulation.
- `Σ2`:           Variance Covariance matrix for the AR model with a MvNormal distribution for the noise. Default value is an identity Matrix.
- `E`:            Distribution for the error.


# Returns
A multivariate series simulating an AR model each column containing a dimension and ordered by time ascending rows.

# Examples
```julia-repl
julia> arsim(1,1,1,10)
10-element Vector{Float64,1}:
[...]
"""
function arsim(xar::AR, n::Integer,
               fix::Union{Nothing,AbstractVector{Union{Missing,T}}} = nothing) where T<:Real
    Φ = compact(xar.Φ)
    Φ0 = compact(xar.Φ0)
    x0 = compact(zeros(size(Φ,2),size(Φ,3)))
    Σ2 = compact(xar.Σ2)
    arsim(Φ,Φ0,x0,n;Σ2,fix)
end

function arsim(Φ::T, Φ0::T, x0::T, n::Integer;
               Σ2::T = T(1), E::Union{Nothing,Distribution} = Normal(T(0),T(1)),
               fix::Union{Nothing,AbstractVector{Union{Missing,T}}} = nothing) where T<:Real

    x = Array{T}(undef,n)

    if Σ2 == T(0) && isnothing(fix)
        for i in 1:n
            x[i] = x0 = Φ * x0 + Φ0
        end
        return x
    end

    if Σ2 == T(0) && !isnothing(fix)
        @assert size(x) == size(fix)
        for i in 1:n
            x[i] = x0 = ismissing(fix[i]) ? Φ * x0 + Φ0 : fix[i]
        end
        return x
    end

    E =  Σ2 != T(1) ? Normal(0,sqrt(Σ2)) : E
    e = rand(E,n)

    for i in 1:n
        x[i] = x0 = Φ * x0 + Φ0 + e[i]
    end

    return x

end

function arsim(Φ::AbstractVector{T}, Φ0::T, x0::AbstractVector{T}, n::Integer;
               Σ2::T = T(1), E::Union{Nothing,Distribution} = Normal(T(0),T(1)),
               fix::Union{Nothing,AbstractVector{Union{Missing,T}}} = nothing) where T<:Real

    np = length(Φ)

    @assert length(Φ) == length(x0)

    x = Array{T}(undef,n)
    
    if Σ2 == T(0) && isnothing(fix)
        for i in 1:n
            x[i] = sum(Φ .* x0) + Φ0
            x0 = vcat(x[i],x0)[1:np]
        end
        return x
    end

    if Σ2 == T(0) && !isnothing(fix)
        @assert size(x) == size(fix)
        for i in 1:n
            x[i] = ismissing(fix[i]) ? sum(Φ .* x0) + Φ0 : fix[i]
            x0 = vcat(x[i],x0)[1:np]
        end
        return x
    end

    E =  Σ2 != T(1) ? Normal(0,sqrt(Σ2)) : E
    e = rand(E,n)

    for i in 1:n
        x[i] = sum(Φ .* x0) + Φ0 + e[i]
        x0 = vcat(x[i],x0)[1:np]
    end
    
    return x
    
end

function arsim(Φ::AbstractArray{T}, Φ0::AbstractVector{T}, x0::AbstractArray{T}, n::Integer;
               Σ2::Matrix{T} = collect(T(1)*I(length(Φ0))),
               E::Union{Nothing,Distribution} = MvNormal(length(Φ0),1),
               fix::Union{Nothing,AbstractMatrix{Union{Missing,T}}} = nothing) where T<:Real

    m = size(Φ,1)
    np = size(Φ,3)

    @assert length(Φ0) == m
    @assert size(x0,1) == m
    @assert size(x0,2) == np

    mnp = m*np
    x0 = reshape(x0,mnp)
    Φ = reshape(Φ,m,mnp)

    x = Array{T}(undef,n,m)

    if Σ2 == zeros(m,m) && isnothing(fix)
        for i in 1:n
            x[i,:] = Φ * x0 + Φ0
            x0 = vcat(x[i,:],x0)[1:mnp]
        end
        
        return x
    end

    if Σ2 == zeros(m,m) && !isnothing(fix)
        @assert size(x) == size(fix)
        for i in 1:n
            x[i,:] =  Φ * x0 + Φ0
            fixi = findall((!ismissing).(fix[i,:]))
            x[i,fixi] = fix[i,fixi]
            x0 = vcat(x[i,:],x0)[1:mnp]
        end
        return x
    end
    
    E =  Σ2 != collect(I(length(Φ0))) ? MvNormal(Σ2) : E
    e = rand(E,n)

    for i in 1:n
        x[i,:] = Φ * x0 + Φ0 + e[:,i]
        x0 = vcat(x[i,:],x0)[1:mnp]
    end
    
    return x
    
end
