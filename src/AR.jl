"""
Package: Forecast

Store results from the function `ar`

# Arguments
`Φ::AbstractArray`            Collection of d by d matrices of coefficients
`coefficients::AbstractArray` Collection of d by d matrices of coefficients
`Φ0::AbstractArray`           Constant
`constant::AbstractArray`     Constant
`Σ::AbstractArray`            Noise sigma variance/covariance Matrix
`stdev::AbstractArray`        Noise sigma variance/covariance Matrix
`fitted::AbstractArray`       Fitted values
`residuals::AbstractArray`    Prediction Error
`IC::Dict`                    Collection of Information Criteria
`PC::AbstractArray`           Parameters Variance/Covariance
`xt::AbstractArray`           Initial values x_t, x_{t-1} ... x_{t-order+1}
`call::String`                method called to generate AR
"""
struct AR
    Φ
    coefficients
    Φ0
    constant
    Σ
    stdev
    fitted
    residuals
    IC::Dict
    PC
    call::String
end

function Base.show(io::IO, xar::AR)
    printstyled("Multivariate Autoregressive Model",bold=true,color=:underline); println()
    println(xar.call)
    println(); printstyled("Residuals Summary",bold=true,color=:underline); println()
    display(summarize(xar.residuals))
    println(); printstyled("Φ0 Constant Coefficient",bold=true,color=:underline); println()
    display(xar.Φ0)
    println();printstyled("Φi Coefficients",bold=true,color=:underline); println()
    display(xar.Φ)
    println(); printstyled("Φ Coefficients' Std. Error",bold=true,color=:underline); println()
    sΦi = []
    nd = ndims(xar.Φ)
    sz = size(xar.Φ)
    for i in 1:(nd <= 1 ? 1 : sz[1])
        for j in 0:(nd == 0 ? 1 : ( nd <= 2 ? sz[1] : sz[1]*sz[3]))
            push!(sΦi, "Φ["*string(i)*","*string(j)*"]:")
        end
    end
    display(hcat(sΦi,sqrt.(abs.(diag(xar.PC)))))
    println(); printstyled("Σ Noise Std. Deviation",bold=true,color=:underline); println()
    display(xar.Σ)
    println(); printstyled("Information Criteria",bold=true,color=:underline); println()
    display(xar.IC)
end
