"""
Package: Forecast

    forecast(xar, n; levels = (0.8,.95))

Forecast a univariate or multivariate autoregressive model.
    
The forecasting follows the model:

```math
Xt = \\Phi_0 + \\sum_{i=1}^p \\Phi_i \\cdot X_{t-i} + E
```

# Arguments
`xar`           AR struct coming from the `ar` function.
`n`             Number of time periods to be forecasted.
`alpha`         Prediction intervals levels; its default value is (0.8, 0.95)
`fixMean`       Fixes the mean in the forecast with a DataFrame which first 
                column is a timestamp type and missing values indicate values
                to be estimated. Default value is `nothing`.
`fixΣ2`         fixes Σ2 values in the forecast

# Returns
A FORECAST struct
"""
function forecast(xar, n::Integer;
                  alpha::Tuple{Real,Real} = (0.8,.95),
                  fixMean::Union{Nothing,DataFrame} = nothing,
                  fixΣ2::Union{T,AbstractMatrix{T}} = xar.Σ2) where T<:Real

    @assert n > 0 "n must be greater than 0"

    m = xar.ndims
    np = xar.order
    
    Φ,Φ0 = xar.Φ, xar.Φ0
     
    dfts = xar.x = tots(xar.x)
    names_x = names(dfts)[2:end]
    x = Array(dfts[:,2:end])
    ts = dfts[:,1]
    name_ts = names(dfts)[1]

    x0 = compact(x[end:-1:end-np+1,:]')
    PT = promote_type(eltype(Φ),eltype(Φ0))
    x0 = PT.(x0)
    
    Σ20 = compact(zeros(m,m))
    !isnothing(fixMean) && (fixMean = Array(fixMean[:,2:end]))
    mu = arsim(Φ,Φ0,x0,n; Σ2=Σ20, fix=fixMean)
    se = sqrt.(fvar(Φ,compact(fixΣ2),max(n,np)))[1:n,:]
    
    # Prediction Intervals
    a1 = alpha[1]
    a2 = alpha[2]

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

    # mu
    mu_names = ["mean_" * n  for n in names_x]
    mu_df = hcat(DataFrame(nΔt(ts,n),[name_ts]),
                 DataFrame(reshape(mu,:,size(mu,2)),mu_names))

    # upper
    u1_names = ["upper1_" * n  for n in names_x]
    upper1_df = DataFrame(reshape(upper1,:,size(upper1,2)),u1_names)
    u2_names = ["upper2_" * n  for n in names_x]
    upper2_df = DataFrame(reshape(upper2,:,size(upper2,2)),u2_names)
    upper_df = hcat(mu_df[:,1:1], upper1_df, upper2_df)
    
    # lower
    l1_names = ["lower1_" * n  for n in names_x]
    lower1_df = DataFrame(reshape(lower1,:,size(lower1,2)),l1_names)
    l2_names = ["lower2_" * n  for n in names_x]
    lower2_df = DataFrame(reshape(lower2,:,size(lower2,2)),l2_names)
    lower_df = hcat(mu_df[:,1:1],lower1_df, lower2_df)

    # se
    se_names = ["se_" * n  for n in names_x]
    se_df = DataFrame(reshape(se,:,size(se,2)),se_names)
    se_df = hcat(mu_df[:,1:1],se_df)
    
    return FORECAST(xar, alpha,
                    mu_df,
                    upper_df,
                    lower_df,
                    se_df,
                    "Forecasting " * xar.call)
end

"""
Forecast recursive variance/covariance
"""
function fvar(Φ::AbstractArray{T},
              Σ2::AbstractMatrix{T},
              n::Integer) where T<:Real

    m = size(Φ,1)
    np = size(Φ,3)
    
    Φ = reshape(Φ,m,m,np)

    Σ2i = zeros(T,m,m)
    v = Array{T}(undef,n,m)
    
    Φ = reshape(Φ,m,m*np)
    for i in 1:n
        #Σ2i = Φ * Σ2i * Φ' + Σ2
        Σ2i = Φ * repeat(repeat(Σ2i,1,np) * Φ',np,1) + Σ2
        v[i,:] = diag(Σ2i)
    end

    v
end

function fvar(Φ::AbstractVector{T},
              Σ2::T,
              n::Integer) where T<:Real

    np = length(Φ)
    Σ2i = T(0)
    v = Vector{T}(undef,n)
    Φ2 = Φ.^2
    
    for i in 1:n
        v[i] = Σ2i = sum(Φ2 .* Σ2i) + Σ2
    end

    v
end

function fvar(Φ::T,
              Σ2::T,
              n::Integer) where T<:Real

    Σ2i = T(0)
    v = Vector{T}(undef,n)
    Φ2 = Φ^2
    
    for i in 1:n
        v[i] = Σ2i = Φ2 * Σ2i + Σ2
    end

    v
end

