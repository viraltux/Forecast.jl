"""
Package: Forecast

    ar(x::DataFrame, order, constant = true; method = "ols")
    ar(x::AbstractArray, order, constant = true; method = "ols", varnames = nothing)

Fit a multivariate autoregressive series model.
    
The fitting is done via Ordinary Least Squared and implements the following model:

```math
Xt = \\Phi_0 + \\sum_{i=1}^p \\Phi_i \\cdot X_{t-i} + \\mathcal{N}(\\vec{0},\\Sigma)
```

# Arguments
- `x`: Multivariate series each column containing a dimension and ordered by time ascending rows.
- `order`: Number of parameters Φ to be estimated.
- `constant`: If `true` `ar` estimates Φ0 otherwise it is assume to be zero.
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
function ar(df::DataFrame, order::Integer = 1, constant::Bool = true;
            dΦ0 = nothing, dΦ = nothing)

    dfx =df[:,eltype.(eachcol(df)) .<: Real]
    xar = ar_ols(Array(x), order, constant; dΦ0 = dΦ0, dΦ = dΦ, varnames = names(dfx))
    xar.x = df

    return(xar)
    
end

function ar(x::AbstractArray, order::Integer = 1, constant::Bool = true;
            dΦ0 = nothing, dΦ = nothing, varnames = nothing)

    return ar_ols(x, order, constant; dΦ0 = dΦ0, dΦ = dΦ, varnames = varnames)
    
end

function ar_ols(x::AbstractArray, order::Integer, constant::Bool; 
               dΦ0 = nothing, dΦ = nothing, varnames)

    @assert 1 <= ndims(x) <= 2
    @assert 1 <= order < size(x,1)-1

    #x = x[end:-1:1,:]

    n = size(x,1)
    m = size(x,2)
    p = order

    r0 = repeat([0.0],m)
    r1 = repeat([1.0],m)
    dΦ0 = isnothing(dΦ0) ? (r1,r1) : dΦ0
    dΦ0 = constant ? dΦ0 : (r1,r0)
    r11 = reshape(repeat([1.0],p*m*m),m,m,p)
    dΦ  = isnothing(dΦ)  ? (r11,r11) : dΦ

    varnames = isnothing(varnames) ? ["x"*string(i) for i in 1:m] : varnames

    # n = 100; m = 2; p = 1
    # x = hcat(1:100,-1:-1:-100) # 1 100 -100 99 -99
    # x = x + rand(100,2)/1000
    # # #x = hcat(1:100)

    x = x[end:-1:1,:]
    nx = x
    for i in 1:p
        nx = hcat(x,nx[vcat(2:end,1),:])
    end
    nx = nx[1:end-p,:]

    Y = nx[:,1:m]
    X = nx[:,m+1:end]
    X = hcat(repeat([1.0],size(X,1)),X)
    #W=(X'*X)\(X'*Y)
    
    # What is M?
    # M = Array{Float64,3}(undef,(n-p, p+1, m))
    # for i in 1:m
    #     for j in 1:n-p
    #         M[j,:,i] = x[j:j+p,i]
    #     end
    # end
    # Y = M[:,1,:]
    # X = reshape(M[:,2:end,:],(n-p,p*m))
    # X = hcat(repeat([1.0],inner=n-p),X)

    # Fixing parameters
    W = Array{Float64,2}(undef,(p*m+1,m))
    for i in 1:m
        Xi,Yi = fixΦ(X,Y,i,dΦ0,dΦ)
        dW = (Xi'*Xi)\(Xi'*Yi)
        W[:,i] = fixW(dW,i,dΦ0,dΦ)
    end

    # Fitted values
    fitted = X*W
    
    # Prediction error
    residuals = Y - fitted

    # Maximum Likelihood noise covariance
    k = p*m*m + m # +m for Φ0
    Σ2 = 1/(n-k)*residuals'*residuals
 
    # ML parameters std. error.
    Φse = sqrt.(abs.(diag(kron(Σ2, (X'*X)^-1))))
    Φse = fixΦse(Φse,dΦ0,dΦ)

    # Information Criteria
    lΣ2   = log(norm(Σ2))
    ic = Dict([("AIC",  lΣ2 + 2*p*m^2/n),
               ("AICC", lΣ2 + 2*(p*m^2+1)/(n-(p*m^2+2))),
               ("BIC",  n*log(norm(Σ2)+m)+(m^2*p+m*(m+1)/2)*log(n)),
               ("H&Q",  lΣ2 + 2*log(log(n))*p*m^2/n)])

    Φ0 = W[1,:]
    Φ = W[2:end,:]
    Φ = reshape(Φ',m,m,p)

    coefficients = Φ
    ar_constant = Φ0
    stdev = Σ = compact(real(sqrt(Σ2)))

    rΦse =  reshape(Φse,:,m)'
    Φ0se = rΦse[:,1]
    Φse = reshape(rΦse[:,2:end],m,m,:)
    p0se = Φ0se
    pse = Φse
    
    call = "ar(X, order="*string(order)*
        ", constant="*string(constant)*")"
    
    AR(varnames,
       Φ,coefficients,
       Φ0,ar_constant,
       Σ,stdev, 
       x,fitted,residuals,
       ic,Φse,pse,Φ0se,p0se,call)

end

"""
Package: Forecast

    fixΦ(X,Y,dΦ0,dΦ)

For a given X and Y OLS matrices returns the X and Y resulting from fixing parameters given dΦ0 and dΦ
"""
function fixΦ(X,Y,i,dΦ0,dΦ)

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

For a given Weight matrix returns a version with fixed values xbased on dΦ0 and dΦ
"""
function fixW(dWi,i,dΦ0,dΦ)
    
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

    fixΦse(M,dΦ0,dΦ)

For a given SE matrix returns an version with zeroes based on dΦ0 and dΦ
"""
function fixΦse(Φse,dΦ0,dΦ)
    
    Φ0, fΦ0 = dΦ0
    Φ, fΦ = dΦ
    (Φ0 == fΦ0) & (Φ == fΦ) && return(Φse)

    m,p = arsize(Φ)
    
    rΦ  = reshape(hcat(Φ0 ,reshape(Φ ,m,m*p))',m*(m*p+1),1)
    rfΦ = reshape(hcat(fΦ0,reshape(fΦ,m,m*p))',m*(m*p+1),1)

    dc = findall(rΦ .!== rfΦ)
    
    fΦse = Vector{Float64}(undef, length(Φse))
    fΦse[:] = Φse
    for ci in dc
        c = Tuple(ci)
        fΦse[c[1]] = 0.0
    end

    fΦse

end
