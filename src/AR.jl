"""
Package: Forecast

Store results from the function `ar`

# Arguments
`names`           List of variable names
`Φ`               Collection of d by d matrices of coefficients
`coefficients`    Copy of Φ
`Φ0`              Constant
`constant`        Copy of Φ0
`Σ`               Noise sigma variance/covariance Matrix
`stdev`           Copy of Σ
`x`               Original dataset
`fitted`          Fitted values
`residuals`       Prediction Error
`IC::Dict`        Collection of Information Criteria
`PC`              Parameters Variance/Covariance
`xt`              Initial values x_t, x_{t-1} ... x_{t-order+1}
`call::String`    Method called to generate AR
"""
mutable struct AR
    varnames
    Φ
    coefficients
    Φ0
    constant
    Σ
    stdev
    x
    fitted
    residuals
    IC::Dict
    PC
    call::String
end

function Base.show(io::IO, xar::AR)

    m,p = arsize(xar.Φ)

    printstyled("Multivariate Autoregressive Model\n",bold=true,color=:underline)
    print("\n    ",xar.call,"\n")
    
    printstyled("\nResiduals Summary\n",bold=true,color=:underline)
    xar_sum = summarize(xar.residuals)
    pretty_table(xar_sum.quantiles, nosubheader = true, show_row_number=false)
    pretty_table(xar_sum.moments, nosubheader = true, show_row_number=false)

    printstyled("Φ0 Constant Coefficient\n",bold=true,color=:underline)
    display(xar.Φ0)

    printstyled("\nΦi Coefficients\n",bold=true,color=:underline)
    display(xar.Φ)

    printstyled("\nΦ[dim,pos] Coefficients' Std. Error\n",bold=true,color=:underline)
    sΦi = []
    nd = ndims(xar.Φ)
    sz = nd == 0 ? 1 : size(xar.Φ)
    for i in 1:(nd <= 1 ? 1 : sz[1])
        for j in 0:(nd == 0 ? 1 : ( nd <= 2 ? sz[1] : sz[1]*sz[3]))
            push!(sΦi, "Φ["*string(i)*","*string(j)*"]:")
        end
    end

    constant = length(xar.Φ) != length(xar.PC)
    
    se = vcat(constant ? [] : 0 , sqrt.(abs.(diag(xar.PC))))
    mu = ndims(xar.Φ) == 0 ? xar.Φ : reshape(hcat(xar.Φ0,reshape(xar.Φ,m,m*p))',m*(m*p+1),1)
    mu = abs.(mu)
    
    function sigf(x)
        if x < 0.001 return  "***" end
        if x < 0.01  return  "** " end
        if x < 0.05  return  "*  " end
        if x < 0.1   return  ".  " end
        if x < 1     return  "   " end
    end

    sig = sigf.(cdf.(Normal.(mu,se),0))

    pretty_table(hcat(sΦi,se,sig), vlines = :none, noheader = true, show_row_number=false)
    print("Signif. codes: 0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1\n")
    
    printstyled("\nΣ Noise Std. Deviation\n",bold=true,color=:underline)
    display(xar.Σ)

    printstyled("\nInformation Criteria\n",bold=true,color=:underline)
    pretty_table(DataFrame(xar.IC), nosubheader = true, show_row_number=false)

end


   # Coefficients:
   #             Estimate Std. Error t value Pr(>|t|)    
   # (Intercept)   0.4809     0.0288    16.7   <2e-16 ***
   # ---
   # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  
