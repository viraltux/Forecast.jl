"""
Package: Forecast

Store results from the function `ar`

# Arguments
`varnames`        List of variable names
`Φ`               Collection of d by d matrices of coefficients
`coefficients`    Alias for Φ
`Φ0`              Constant
`constant`        Alias for Φ0
`Σ`               Noise sigma variance/covariance Matrix
`stdev`           Alias for Σ
`x`               Original dataset
`fitted`          Fitted values
`residuals`       Prediction Error
`ic::Dict`        Collection of Information Criteria
`Φse`             Parameters Standard Error
`pse`             Alias for Φse
`Φ0se`            Cosntant Standard Error
`p0se`            Alias for Φ0se
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
    ic::Dict
    Φse
    pse
    Φ0se
    p0se
    call::String
end

function Base.show(io::IO, xar::AR)

    m,p = arsize(xar.Φ)
    printstyled("Multivariate Autoregressive Model\n",bold=true,color=:underline)
    print("\n    ",xar.call,"\n")
    
    printstyled("\nResiduals Summary\n",bold=true,color=:underline)
    xar_sum = summarize(xar.residuals; varnames = xar.varnames)
    pretty_table(xar_sum.quantiles, nosubheader = true, show_row_number=false)
    pretty_table(xar_sum.moments, nosubheader = true, show_row_number=false)

    printstyled("\nCoefficients\n\n",bold=true,color=:underline)
    Φ0 = xar.Φ0
    Φ0se = xar.Φ0se
    printstyled("Φ0\n",bold=true,color=:underline)
    pretty_table(string.(round.(Φ0,digits=3)) .* sigf.(Φ0,Φ0se), tf = tf_matrix, noheader=true)

    for i in 1:p
        printstyled("Φ",i,"\n",bold=true,color=:underline)
        Φi = xar.Φ[:,:,i]
        Φsei = xar.Φse[:,:,i]
        pretty_table(string.(round.(Φi,digits = 3)) .* sigf.(Φi,Φsei), tf = tf_matrix, noheader=true)
    end
    print("Signif. codes: 0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘^’ 0.1 ‘ ’ 1\n")


    printstyled("\nCoefficients' Std. Error\n",bold=true,color=:underline)
    printstyled("\nΦ0se\n",bold=true,color=:underline)
    pretty_table(xar.Φ0se, tf = tf_matrix, noheader=true)
    for i in 1:p
        printstyled("Φ",i,"se\n",bold=true,color=:underline)
        pretty_table(xar.Φse[:,:,i], tf = tf_matrix, noheader=true)    
    end

    printstyled("\nΣ Noise Std. Deviation\n",bold=true,color=:underline)
    Σ = m == 1 ? [xar.Σ] : xar.Σ
    pretty_table(Σ, tf = tf_matrix, noheader=true)
    
    printstyled("\nInformation Criteria\n",bold=true,color=:underline)
    pretty_table(DataFrame(xar.ic), nosubheader = true, show_row_number=false)

end

function arsize(Φ)
    d = ndims(Φ)
    sΦ = size(Φ)
    if d == 0
        m,p = 1,1
    elseif d == 1
        m,p = 1,sΦ[1]
    elseif d == 2
        m,p = sΦ[1],1
    elseif d == 3
        m,p = sΦ[1],sΦ[3]
    else
        @error "Φ should have less the 4 dimensions"
    end
    (m,p)
end

function sigf(mu,se)
    if se == 0.0  return  "fixed" end
    
    pv = cdf(Normal(abs(mu),se),0)
    if pv < 0.001 return  " ***"   end 
    if pv < 0.01  return  " ** "   end 
    if pv < 0.05  return  " *  "   end 
    if pv < 0.1   return  " ^  "   end 
    if pv <= 1    return  "    "   end 
    
end
