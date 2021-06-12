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

# Returns
A FORECAST struct
"""
function forecast(xar::AR, n::Integer;
                  alpha = (0.8,.95),
                  fix = nothing)

    @assert n > 0 "n must be greater than 0"

    Φ,Φ0,Σ2 = compact(xar.Φ),compact(xar.Φ0),compact(xar.Σ2)
    fix = isnothing(fix) ? fix : Array(fix[:,2:end])
    
    m,np = arsize(Φ)
    
    dfts = xar.x = tots(xar.x)
    names_x = names(dfts)[2:end]
    x = Array(dfts[:,2:end])
    ts = dfts[:,1]
    name_ts = names(dfts)[1]

    x0 = compact(x[end:-1:end-np+1,:]')
    Σ0 = compact(zeros(m,m))
    mu = arsim(Φ,Φ0,x0,n; Σ=Σ0, fix)
    se = sqrt.(fvar(Φ,Σ2,max(n,np)))[1:n,:]
    
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
Forecast recurrence variance for a linear combination
"""
function fvar(Φ,Σ2,n)

    m,np = arsize(Φ)

    #TODO create fvar for univariate Φ
    Φ = Φ isa Number ? [Φ] : reshape(Φ,m,m,np)
    Σ2 = Σ2 isa Number ? [Σ2] : Σ2

    Φ2 = similar(Φ)
    for i in 1:np
        Φ2[:,:,i] =  Φ[:,:,i]^2
    end
    Φ2 = reshape(Φ2,m,m*np)
    
    Σ2i =zeros(m*np,m)
    v = zeros(n,m)

    for i in 1:n
        Σ2i[1:m,:] = Φ2 * Σ2i + Σ2
        v[i,:] = diag(Σ2i)
        Σ2i = circshift(Σ2i,m)
    end

    #TODO create fvar for univariate Φ
    return compact(v)
end
