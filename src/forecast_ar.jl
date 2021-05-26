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
function forecast(xar::AR, n::Integer; alpha = (0.8,.95))

    @assert n > 0 "n must be greater than 0"
    
    Φ,Φ0,Σ = xar.Φ,xar.Φ0,xar.Σ

    m,np = arsize(Φ)

    dfts = xar.x = tots(xar.x)
    names_x = names(dfts)[2:end]
    x = Array(dfts[:,2:end])
    ts = dfts[:,1]
    name_ts = names(dfts)[1]

    x0 = compact(x[end:-1:end-np+1,:]')
    E = MvNormal(m,0)
    mu = arsim(Φ,Φ0,x0,E,n)
    se = sqrt.(fvar(Φ,Σ,n))

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
    mu_names = ["Mean_" * n  for n in names_x]
    mu_df = hcat(DataFrame(nΔt(ts,n),[name_ts]),
                 DataFrame(reshape(mu,:,size(mu,2)),mu_names))

    # upper
    u1_names = ["Upper1_" * n  for n in names_x]
    upper1_df = DataFrame(reshape(upper1,:,size(upper1,2)),u1_names)
    u2_names = ["Upper2_" * n  for n in names_x]
    upper2_df = DataFrame(reshape(upper2,:,size(upper2,2)),u2_names)
    upper_df = hcat(mu_df[:,1:1], upper1_df, upper2_df)
    
    # lower
    l1_names = ["Lower1_" * n  for n in names_x]
    lower1_df = DataFrame(reshape(lower1,:,size(lower1,2)),l1_names)
    l2_names = ["Lower2_" * n  for n in names_x]
    lower2_df = DataFrame(reshape(lower2,:,size(lower2,2)),l2_names)
    lower_df = hcat(mu_df[:,1:1],lower1_df, lower2_df)

    return FORECAST(xar, alpha,
                    mu_df,
                    upper_df,
                    lower_df,
                    "Forecasting " * xar.call)
end

"""
Forecast recurrence variance for a linear combination
"""
function fvar(Φ,Σ,n)
    Φ = compact(Φ)

    m,np = arsize(Φ)

    v = zeros(n,m)
    v[1,:] = m <= 1 ? [Σ^2] : diag(Σ^2)
    # reshape and format for matrix operation
    Φ = Φ isa Number ?  reshape([Φ],1,1) : reshape(Φ,m,m*np)
    
    for i in 2:n
        for j in 1:m
            if i <= np
                v[i,j] = (Φ[j,:].^2)' * repeat(v[1:np,j],inner=m) +  v[1,j]
            else
                v[i,j] = (Φ[j,:].^2)' * repeat(v[i-np:i-1,j],inner=m) +  v[1,j]
            end
        end
    end
    
    compact(v)
end

"""
Transform a FORECAST object value with given function
"""
function transform(fc::FORECAST, f::Function)

    ts_x = tots(fc.model.x)
    names_x = names(ts_x)
    names_x[2:end] = "$(string(f))_" .* names_x[2:end]

    ts_mean = tots(fc.mean)
    ts_lower = tots(fc.lower)
    ts_upper = tots(fc.upper)
    names_mean = names(ts_mean)
    names_mean[2:end] = "$(string(f))_" .* names_mean[2:end]
    names_lower = names(ts_lower)
    names_lower[2:end] = "$(string(f))_" .* names_lower[2:end]
    names_upper = names(ts_upper)
    names_upper[2:end] = "$(string(f))_" .* names_upper[2:end]

    xts = ts_x[:,1]
    x = f.(Array(ts_x[:,2:end]))

    fmean = f.(Array(ts_mean[:,2:end]))
    flower = f.(Array(ts_lower[:,2:end]))
    fupper = f.(Array(ts_upper[:,2:end]))

    n = size(x,1)
    m = size(x,2)

    pfc = deepcopy(fc)
    pfc.model.varnames = names_x[2:end]

    pfc.model.x = hcat(ts_x[:,1:1], DataFrame(x,names_x[2:end]))
    pfc.mean    = hcat(ts_mean[:,1:1], DataFrame(fmean,names_mean[2:end]))
                        

    if size(fmean,2) > 1
        z = fupper .- repeat(fmean,1,2)
        fmean = repeat(pfmean,1,2)
    else
        z = fupper .- fmean
    end
    
    pfc.upper = hcat(ts_mean[:,1:1],
                     DataFrame(fmean .+ z,names_upper[2:end]))
    pfc.lower = hcat(ts_mean[:,1:1],
                     DataFrame(fmean .- z,names_lower[2:end]))
    
    pfc.call = fc.call * "\nData transformed with function: $(f)"

    return(pfc)
end
