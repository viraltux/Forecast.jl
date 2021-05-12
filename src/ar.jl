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


    @assert 1 <= length(size(x)) <= 2
    n = length(x[:,1])
    @assert 1 <= order < n-1

    m = length(size(x)) == 2 ? size(x)[2] : 1 # dimension
    p = order

    
    dΦ0 = constant ? dΦ0 : (repeat([1.0],m),repeat([0.0],m))
    dΦ0 = isnothing(dΦ0) ? (1,1) : dΦ0
    dΦ  = isnothing(dΦ)  ? (1,1) : dΦ

    varnames = isnothing(varnames) ? ["x"*string(i) for i in 1:m] : varnames
    
    M = Array{Float64,3}(undef,(n-p, p+1, m))
    for i in 1:m
        for j in 1:n-p
            M[j,:,i] = x[j+p:-1:j,i]
        end
    end

    Y = M[:,1,:]
    X = reshape(M[:,2:end,:],(n-p,p*m))
    #X = constant ? hcat(repeat([1.0],inner=n-p),X) : X
    X = hcat(repeat([1.0],inner=n-p),X)

    X,Y = fixΦ(X,Y,dΦ0,dΦ)

    W = (X'*X)\(X'*Y)

    # Maximum Likelihood noise covariance
    k = p*m*m
    Σ2 = 1/(n-k)*(Y-X*W)'*(Y-X*W)

    # Fitted values
    fitted = X*W
    
    # ML parameters covariance
    PC = kron(Σ2, (X'*X)^-1)

    # Prediction error
    residuals = Y - fitted

    # Information Criteria
    lΣ2   = log(norm(Σ2))
    IC = Dict([("AIC",  lΣ2 + 2*p*m^2/n),
               ("AICC", lΣ2 + 2*(p*m^2+1)/(n-(p*m^2+2))),
               ("BIC",  n*log(norm(Σ2)+m)+(m^2*p+m*(m+1)/2)*log(n)),
               ("H&Q",  lΣ2 + 2*log(log(n))*p*m^2/n)])

    @bp
    W = fixM(W,dΦ0,dΦ,true,false)

    Φ0 = constant ? W[1,:] : repeat([0.0],inner=m)

    Φ = Array{Float64,3}(undef,(m,m,p))
    for (i,j,k) in zip(repeat(1:m,inner=p),repeat(1:p,m),1:p*m)
        # Φ[:,i,j] = W[constant ? k+1 : k,:]
        Φ[:,i,j] = W[k+1,:]
    end

    PC = fixM(PC,dΦ0,dΦ)

    coefficients = Φ
    ar_constant = Φ0
    stdev = Σ = compact(real(sqrt(Σ2)))

    call = "ar(X, order="*string(order)*
        ", constant="*string(constant)*")"
    
    AR(varnames,
       Φ,coefficients,
       Φ0,ar_constant,
       Σ,stdev, 
       x,fitted,residuals,
       IC,PC,call)

end

"""
Package: Forecast

    fixΦ(X,Y,dΦ0,dΦ)

For a given X and Y OLS matrices returns the X and Y resulting from fixing parameters given dΦ0 and dΦ
"""
function fixΦ(X,Y,dΦ0,dΦ)

    Φ0, fΦ0 = dΦ0
    Φ, fΦ = dΦ
    (Φ0 == fΦ0) & (Φ == fΦ) && return((X,Y))

    rΦ0 = Φ0 isa Number ? Φ0 : reshape(Φ0,:,1)
    rfΦ0 = fΦ0 isa Number ?  fΦ0 : reshape(fΦ0,:,1)
    rΦ = Φ isa Number ? Φ : reshape(Φ,:,1)
    rfΦ = fΦ isa Number ?  fΦ : reshape(fΦ,:,1)

    rΦ = compact([rΦ0; rΦ])
    rfΦ = compact([rfΦ0; rfΦ])

    dc = findall(rΦ .!== rfΦ)

    @assert length(rΦ) > length(dc) "No degrees of freedom"
    
    # Create new Y
    for i in dc
        Y = Y - rfΦ[i]*X[:,i]
    end

    # Create new X
    X = drop(X,c=dc)

    (X,Y)
    
end

# dΦ0 could be part of dΦ, it simply would mean Φ0 is fixed to zero when no cosntant
# is required.

"""
Package: Forecast

    fixM(M,dΦ0,dΦ)

For a given matrix M returns an expanded with zeroes version of M based on dΦ0 and dΦ
"""
function fixM(M,dΦ0,dΦ,row=true,col=true)

    Φ0, fΦ0 = dΦ0
    Φ, fΦ = dΦ
    (Φ0 == fΦ0) & (Φ == fΦ) && return(M)

    rΦ0 = Φ0 isa Number ? Φ0 : reshape(Φ0,:,1)
    rfΦ0 = fΦ0 isa Number ?  fΦ0 : reshape(fΦ0,:,1)
    rΦ = Φ isa Number ? Φ : reshape(Φ,:,1)
    rfΦ = fΦ isa Number ?  fΦ : reshape(fΦ,:,1)

    rΦ = compact([rΦ0; rΦ])
    rfΦ = compact([rfΦ0; rfΦ])
    
    dc = findall(rΦ .!== rfΦ)

    @assert length(rΦ) > length(dc) "No degrees of freedom"

    cM = copy(M)
    cM = reshape(cM,size(cM,1),size(cM,2))
    # Create new Y
    for i in dc
        if !row & col
            cM = insert_col(cM,i,rfΦ[i])
        elseif row & !col
            cM = insert_row(cM,i,rfΦ[i])
        else            
            cM = insert_cross(cM,i,0)
        end
    end

    return(cM)
    
end



