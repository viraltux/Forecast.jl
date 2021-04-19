"""
Package: Forecast

    ar(x::AbstractArray, order::Integer, constant = true; method = "ols")

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

# Returns
An AR object containing the model coefficients, the prediction variance/covariance matrix, residuals a collection of information criteria
"""
function ar(x::AbstractArray, order::Integer, constant = true; method = "ols")

    return ar_ols(x, order, constant)
    
end

function ar_ols(x::AbstractArray, order::Integer, constant = true)

    @assert 1 <= length(size(x)) <= 2
    n = length(x[:,1])
    @assert 1 <= order < n-1
    
    m = length(size(x)) == 2 ? size(x)[2] : 1 # dimension
    p = order

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

    Φ0= constant ? W[1,:] : repeat([0.0],inner=m)

    Φ = Array{Float64,3}(undef,(m,m,p))
    for (i,j,k) in zip(repeat(1:m,inner=p),repeat(1:p,m),1:p*m)
        Φ[:,i,j] = W[constant ? k+1 : k,:]
    end

    # Maximum Likelihood noise covariance
    k = p*m*m
    Σ = 1/(n-k)*(Y-X*W)'*(Y-X*W)

    # ML parameters covariance
    PC = kron(Σ, (X'*X)^-1)

    # Prediction Error
    residuals = Y - X*W

    # `order` most recent values of x time series 
    xt = x[(end-order+1):end,:]

    # Information Criteria
    lΣ   = log(norm(Σ))
    IC = Dict([("AIC",  lΣ + 2*p*m^2/n),
               ("AICC", lΣ + 2*(p*m^2+1)/(n-(p*m^2+2))),
               ("BIC",  n*log(norm(Σ)+m)+(m^2*p+m*(m+1)/2)*log(n)),
               ("H&Q",  lΣ + 2*log(log(n))*p*m^2/n)])
    
    coefficients = Φ
    constant = Φ0
    variance = Σ

    ar_ols(x, order, constant)
    call = "ar(X, order="*string(order)*
        ", constant="*string(constant)*")"
    
    AR(Φ,coefficients,
       Φ0,constant,
       Σ,variance, 
       residuals,
       IC,PC,xt,call)

end
