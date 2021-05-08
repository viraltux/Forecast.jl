"""
Package: Forecast

    ar(x::TimeArray, order, constant = true; method = "ols")
    ar(x::DataFrame, order, constant = true; method = "ols")
    ar(x::AbstractArray, order, constant = true; method = "ols", varnames = nothing)

Fit a multivariate autoregressive series model.
    
Currently only the Ordinary Least Squared method is implemented and fits the following model:

```math
Xt = \\Phi_0 + \\sum_{i=1}^p \\Phi_i \\cdot X_{t-i} + \\mathcal{N}(\\vec{0},\\Sigma)
```

# Arguments
- `x`: Multivariate series each column containing a dimension and ordered by time ascending rows.
- `order`: Number of parameters Φ to be estimated.
- `constant`: If `true` `ar` estimates Φ0 otherwise it is assume to be zero.
- `method`: Method to fit the `ar` model. Currently only "ols".
- `varnames`: Names of the dimensions (by default xi where `i` is an integer)

# Returns
An AR object containing the model coefficients, the error sigma matrix, residuals and a collection of information criteria

# Examples
```julia-repl
julia> ar(rand(100,2),2)
AR([...])
"""
function ar(ta::TimeArray, order::Integer = 1, constant::Bool = true; 
            method = "ols")

    vx = values(ta)
    @assert sum((x -> x isa Number).(vx)) == size(vx,1)*size(vx,2)
        "All values in the time series must be numeric"

    ar_ts = ar_ols(vx, order, constant; varnames = string.(colnames(ta)))
    ar_ts.x = ta
    return ar_ts
end

function ar(df::DataFrame, order::Integer = 1, constant::Bool = true;
            method = "ols")

    dfx = df[:,eltype.(eachcol(df)) .<: Real]
    x = Array(dfx)
    @assert sum((x -> x isa Number).(x)) == size(x,1)*size(x,2)
    "All values in the time series must be numeric"

    return ar_ols(x, order, constant; varnames = names(dfx))
    
end

function ar(x::AbstractArray, order::Integer = 1, constant::Bool = true;
            method = "ols", varnames = nothing)

    return ar_ols(x, order, constant; varnames = varnames)
    
end

function ar_ols(x::AbstractArray, order::Integer, constant::Bool; varnames)

    @assert 1 <= length(size(x)) <= 2
    n = length(x[:,1])
    @assert 1 <= order < n-1
    
    m = length(size(x)) == 2 ? size(x)[2] : 1 # dimension
    p = order

    varnames = isnothing(varnames) ? ["x"*string(i) for i in 1:m] : varnames
    
    M = Array{Float64,3}(undef,(n-p, p+1, m))
    for i in 1:m
        for j in 1:n-p
            M[j,:,i] = x[j+p:-1:j,i]
        end
    end

    Y = M[:,1,:]
    X = reshape(M[:,2:end,:],(n-p,p*m))
    X = constant ? hcat(repeat([1.0],inner=n-p),X) : X

    W = (X'*X)\(X'*Y)

    Φ0 = constant ? W[1,:] : repeat([0.0],inner=m)

    Φ = Array{Float64,3}(undef,(m,m,p))
    for (i,j,k) in zip(repeat(1:m,inner=p),repeat(1:p,m),1:p*m)
        Φ[:,i,j] = W[constant ? k+1 : k,:]
    end

    # Maximum Likelihood noise covariance
    k = p*m*m
    Σ2 = 1/(n-k)*(Y-X*W)'*(Y-X*W)

    # ML parameters covariance
    PC = kron(Σ2, (X'*X)^-1)

    # Fitted values
    fitted = X*W
    
    # Prediction error
    residuals = Y - fitted

    # Information Criteria
    lΣ2   = log(norm(Σ2))
    IC = Dict([("AIC",  lΣ2 + 2*p*m^2/n),
               ("AICC", lΣ2 + 2*(p*m^2+1)/(n-(p*m^2+2))),
               ("BIC",  n*log(norm(Σ2)+m)+(m^2*p+m*(m+1)/2)*log(n)),
               ("H&Q",  lΣ2 + 2*log(log(n))*p*m^2/n)])
    
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
