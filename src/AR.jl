"""
Package: Forecast

Store results from the function `ar`

# Arguments
`Φ::AbstractArray`            Collection of d by d matrices of coefficients
`coeficients::AbstractArray`  Collection of d by d matrices of coefficients
`Φ0::AbstractArray`           Constant
`constant::AbstractArray`     Constant
`Σ::AbstractArray`            Variance/Covariance Matrix
`variance::AbstractArray`     Variance/Covariance Matrix
`residuals::AbstractArray`    Prediction Error
`IC::Dict`                    Collection of Information Criteria
`PC::AbstractArray`           Parameters Variance/Covariance
`xt::AbstractArray`           Initial values x_t, x_{t-1} ... x_{t-order+1}
`call::String`                method called to generate AR
"""
struct AR
    Φ::AbstractArray
    coeficients::AbstractArray
    Φ0::AbstractArray
    constant::AbstractArray
    Σ::AbstractArray
    variance::AbstractArray
    residuals::AbstractArray
    IC::Dict
    PC::AbstractArray
    xt::AbstractArray
    call::String
end
