"""
Package: Forecast

    ar(x::DataFrame, or, constant = true;             
                     alpha = 1.0, dΦ0 = nothing, dΦ = nothing)

    ar(x::AbstractArray, or, constant = true; 
                         alpha = false, dΦ0 = nothing, dΦ = nothing, varnames = nothing)

Fit a multivariate autoregressive series model.
    
The fitting is done via Ordinary Least Squared and implements the following model:

```math
Xt = \\Phi_0 + \\sum_{i=1}^p \\Phi_i \\cdot X_{t-i} + \\mathcal{N}(\\vec{0},\\Sigma)
```

# Arguments
- `x`: Multivariate series each column containing a dimension and ordered by time ascending rows.
- `or`: Number of parameters Φ to be estimated.
- `constant`: If `true` `ar` estimates Φ0 otherwise it is assume to be zero.
- `alpha`: fixes to zero all coefficients which significance is below its value. It defaults to one.
- `dΦ0`: Tuple containing two Φ0 objects, the first one will act as an original reference to the second one and different values will be fixed in the fitting process to the values in the second Φ0.
- `dΦ`:  Equivalent to dΦ0 but for Φ.
- `varnames`: Names of the dimensions (by default xi where `i` is an integer)

# Returns
An AR object containing the model coefficients, the error sigma matrix, residuals and a collection of information criteria

# Examples
```julia-repl
julia> ar(rand(100,2),2)
AR([...])
"""
function ar(df::DataFrame,
            order::Integer = 1,
            constant::Bool = true;
            alpha::Real = 1.0,
            dΦ0::Tuple = d_ar_dΦ0(size(df,2)-1,constant,Float64),
            dΦ::Tuple = d_ar_dΦ(size(df,2)-1,order,Float64))

    #TODO promote dΦ0 and dΦ types based on df and control for Integer values on x
    ttype = [Date,DateTime,Day,Month,Week]
    ar_df = eltype(df[:,1]) in ttype ? Array{Float64}(df[:,2:end]) : Array{Float64}(df)
    names_df = eltype(df[:,1]) in ttype ? names(df[:,2:end]) : names(df)
    xar = ar(ar_df, order, constant; alpha, dΦ0, dΦ, varnames = names_df)
    xar.x = df

    return xar
    
end


# Default ar values functions
d_ar_dΦ0(m::Integer,constant::Bool,T::Type) = (repeat([T(1)],m),  repeat([T(constant ? 1 : 0)],m))
d_ar_dΦ(m::Integer,order::Integer,T::Type) =
    (r11=reshape(repeat([T(1)],m*m*order),m,m,order); (r11,r11))
d_ar_varnames(m::Integer)::Vector{String} = ["x"*string(i) for i in 1:m]

function ar(x::AbstractArray{T},
            order::Integer = 1,
            constant::Bool = true;
            alpha::Real = 1.0, 
            dΦ0::Tuple{AbstractArray{T},AbstractArray{T}} = d_ar_dΦ0(size(x,2),constant,T),
            dΦ::Tuple{AbstractArray{T},AbstractArray{T}} = d_ar_dΦ(size(x,2),order,T),
            varnames::AbstractVector{String} = d_ar_varnames(size(x,2))) where T<:Real

    @assert 0.0 < alpha <= 1.0
    
    xar = ar_ols(x, order, constant; dΦ0, dΦ, varnames)

    alpha == 1.0 && return xar
    
    m = xar.ndims
    np = xar.order
    
    dΦs = (reshape(xar.Φpv,m,m,np) .<= alpha)
    dΦ0s = (reshape(xar.Φ0pv,m) .<= alpha)

    all(dΦs) && all(dΦ0s) && return xar

    if isnothing(dΦ)
        dΦ = (xar.Φ, xar.Φ .* dΦs)
        dΦ0 = (xar.Φ0, xar.Φ0 .* dΦ0s)
    else
        dΦ = (dΦ[1], dΦ[2] .* dΦs)
        dΦ0 = (dΦ0[1], dΦ0[2] .* dΦ0s)
    end

    return ar_ols(x, order, constant; dΦ0, dΦ, varnames)

end


function ar_ols(x::AbstractArray{T},
                or::Integer,
                constant::Bool;
                dΦ0::Tuple{AbstractArray{T},AbstractArray{T}},
                dΦ::Tuple{AbstractArray{T},AbstractArray{T}},
                varnames::AbstractVector{String}) where T<:Real

    @assert 1 <= ndims(x) <= 2
    @assert 1 <= or < size(x,1)-1

    n = size(x,1)
    m = size(x,2)
    np = m*m*or + m

    x = x[end:-1:1,:]
    nx = x
    for i in 1:or
        nx = hcat(x,nx[vcat(2:end,1),:])
    end
    nx = nx[1:end-or,:]

    Y = nx[:,1:m]
    X = nx[:,m+1:end]
    X = hcat(repeat([1.0],size(X,1)),X)

    # Fixing parameters
    W = Array{T,2}(undef,(or*m+1,m))
    for i in 1:m
        Xi,Yi = fixΦ(X,Y,i,dΦ0,dΦ)
        dW = (Xi'*Xi)\(Xi'*Yi)
        W[:,i] = fixW(dW,i,dΦ0,dΦ)
    end

    # Fitted values
    fitted = X*W
    
    # Prediction error
    e = Y - fitted

    Φ0 = W[1,:]
    Φ = W[2:end,:]
    Φ = Array(reshape(Φ',m,m,or))

    # Maximum Likelihood noise covariance
    k = or*m*m
    @assert n-k > 0 "Insufficient data for the model"
    Σ2 = variance = 1/(n-k)*e'*e

    # ML parameters std. error.
    Φse = sqrt.(abs.(diag(kron(Σ2, (X'*X)^-1))))
    Φse = fixΦse(Φse,dΦ0,dΦ)

    rΦse = reshape(Φse,:,m)'
    Φ0se = rΦse[:,1]
    Φse  = reshape(rΦse[:,2:end],m,m,:)
    p0se = Φ0se
    pse  = Φse

    #fix p
    np = fixnp(dΦ0,dΦ)

    # TODO check weher to use 'or' or 'np' and if 'or' then a fixor function is needed
    # Information Criteria
    lΣ2   = log(norm(Σ2))
    ic = Dict([("AIC",  lΣ2 + 2*np*m^2/n),
               ("AICC", lΣ2 + 2*(np*m^2+1)/(n-(np*m^2+2))),
               ("BIC",  n*log(norm(Σ2)+m)+(m^2*np+m*(m+1)/2)*log(n)),
               ("H&Q",  lΣ2 + 2*log(log(n))*np*m^2/n)])

    # Statistics
    SStot = sum((Y .- mean(Y,dims=1)).^2,dims=1)
    SSres = sum(e .^ 2,dims=1)
    R2 = vec(1 .- (SSres ./ SStot))
    
    pvf(mu,se) = se == 0 ? 1.0 : cdf(Normal(abs(mu),se),0) #1 to make log(1) = 0
    Φpv = pvf.(Φ,Φse) 
    Φ0pv = pvf.(Φ0,Φ0se)
    p0pv = Φ0pv
    ppv = Φpv
    slpv = vec(reshape(sum(log.(vcat(reshape(Φ0pv,:,m),reshape(Φpv,:,m))),dims=1),:,1))

    stats = Dict([(" Variable", varnames),
                  ("R2",    R2),
                  ("R2adj", 1 .- (1 .- R2) * (n-1)/(n-(np-1)/m-1)),
                  ("Fisher's p-test", 1 .- cdf(Chisq(2*np/m),-2*slpv))])    
    
    coefficients = Φ
    ar_constant = Φ0
    stdev = Σ = sqrt.(diag(Σ2))
    
    call = "ar(X, order="*string(or)*
        ", constant="*string(constant)*")"

    AR(varnames,
       or, m,
       Φ,coefficients,
       Φ0,ar_constant,
       Σ2,variance,
       Σ,stdev,
       x[end:-1:1,:],
       fitted,e,
       ic,stats,
       Φse,pse,Φ0se,p0se,
       Φpv,ppv,Φ0pv,p0pv,
       call)

end

"""
Package: Forecast

    fixΦ(X,Y,dΦ0,dΦ)

For a given X and Y OLS matrices returns the X and Y resulting from fixing parameters given dΦ0 and dΦ
"""
function fixΦ(X::AbstractMatrix{T},
              Y::AbstractMatrix{T},
              i::Integer,
              dΦ0::Tuple{AbstractArray{T},AbstractArray{T}},
              dΦ::Tuple{AbstractArray{T},AbstractArray{T}}) where T<:Real

    Φ0, fΦ0 = dΦ0
    Φ, fΦ = dΦ
    (Φ0 == fΦ0) & (Φ == fΦ) && return((X,Y[:,i]))

    
    rΦ0 = reshape(Φ0,:,1)
    rΦ  = reshape(Φ,size(Φ,1),:)
    rΦ  = hcat(rΦ0,rΦ)

    rfΦ0 = reshape(fΦ0,:,1)
    rfΦ = reshape(fΦ,size(fΦ,1),:)
    rfΦ = hcat(rfΦ0,rfΦ)
    
    dc = findall(rΦ[i,:] .!== rfΦ[i,:])

    @assert length(rΦ) > length(dc) "No degrees of freedom"
    
    # Create new Yi
    Yi = Y[:,i]
    for c in dc
        Yi = Yi .- rfΦ[i,c]*X[:,c]
    end
    
    # Create new X
    Xi = drop(X,c=dc)

    (Xi,Yi)
    
end

"""
Package: Forecast

    fixW(W,dΦ0,dΦ)

For a given Weight matrix returns a version with fixed values based on dΦ0 and dΦ
"""
function fixW(dWi::AbstractVector{T},
              i::Integer,
              dΦ0::Tuple{AbstractArray{T},AbstractArray{T}},
              dΦ::Tuple{AbstractArray{T},AbstractArray{T}}) where T<:Real
    
    Φ0, fΦ0 = dΦ0
    Φ, fΦ = dΦ
    (Φ0 == fΦ0) & (Φ == fΦ) && return(dWi)

    rΦ0 = reshape(Φ0,:,1)
    rΦ  = reshape(Φ,size(Φ,1),:)
    rΦ  = hcat(rΦ0,rΦ)

    rfΦ0 = reshape(fΦ0,:,1)
    rfΦ = reshape(fΦ,size(fΦ,1),:)
    rfΦ = hcat(rfΦ0,rfΦ)
    
    dc = findall(rΦ[i,:] .!== rfΦ[i,:])

    for c in dc
        insert!(dWi,c,rfΦ[i,c])
    end

    dWi
    
end


"""
Package: Forecast

    fixnp(dΦ0,dΦ)

return the number of free parameters
"""
function fixnp(dΦ0::Tuple{AbstractArray{T},AbstractArray{T}},
               dΦ::Tuple{AbstractArray{T},AbstractArray{T}}) where T<:Real
    
    Φ0, fΦ0 = dΦ0
    Φ, fΦ = dΦ
    (Φ0 == fΦ0) & (Φ == fΦ) && return(length(Φ) + length(Φ0))

    rΦ0 = reshape(Φ0,:,1)
    rΦ  = reshape(Φ,size(Φ,1),:)
    rΦ  = hcat(rΦ0,rΦ)

    rfΦ0 = reshape(fΦ0,:,1)
    rfΦ = reshape(fΦ,size(fΦ,1),:)
    rfΦ = hcat(rfΦ0,rfΦ)
    
    dc = findall(rΦ .!== rfΦ)

    length(Φ) + length(Φ0) - length(dc)
end


"""
Package: Forecast

    fixΦse(M,dΦ0,dΦ)

For a given `se` matrix returns an version with zeroes based on dΦ0 and dΦ
"""
function fixΦse(Φse::AbstractVector{T},
                dΦ0::Tuple{AbstractArray{T},AbstractArray{T}},
                dΦ::Tuple{AbstractArray{T},AbstractArray{T}}) where T<:Real
    
    Φ0, fΦ0 = dΦ0
    Φ, fΦ = dΦ
    (Φ0 == fΦ0) & (Φ == fΦ) && return(Φse)

    _,m,np = size(Φ)
    
    rΦ  = reshape(hcat(Φ0 ,reshape(Φ ,m,m*np))',m*(m*np+1),1)
    rfΦ = reshape(hcat(fΦ0,reshape(fΦ,m,m*np))',m*(m*np+1),1)

    dc = findall(rΦ .!== rfΦ)
    
    fΦse = Vector{T}(undef, length(Φse))
    fΦse[:] = Φse
    for ci in dc
        c = Tuple(ci)
        fΦse[c[1]] = 0.0
    end

    fΦse

end
