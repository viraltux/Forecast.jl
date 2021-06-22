"""
Package: Forecast

Store results from the function `ar`

# Arguments
`varnames`        List of variable names
`Φ`               Collection of d by d matrices of coefficients
`coefficients`    Alias for Φ
`Φ0`              Constant
`constant`        Alias for Φ0
`Σ2`              ML variance/covariance Matrix
`variance`        Alias for Σ2
`Σ`               Variables Standard deviation 
`stdev`           Alias for Σ
`x`               Original dataset
`fitted`          Fitted values
`residuals`       Prediction Error
`ic::Dict`        Collection of Information Criteria
`stats::Dict`     Collection of Statistics
`Φse`             Parameters Standard Error
`pse`             Alias for Φse
`Φ0se`            Constant Standard Error
`p0se`            Alias for Φ0se
`Φpv`             p-value for Parameters
`ppv`             Alias for Φpv
`Φ0pv`            p-value for Constant
`p0pv`            Alias for Φ0pv
`call::String`    Method called to generate AR
"""
mutable struct AR{T<:Real}
    varnames::Vector{String}
    Φ::Array{T,3}
    coefficients::Array{T,3}
    Φ0::Vector{T}
    constant::Vector{T}
    Σ2::Matrix{T}
    variance::Matrix{T}
    Σ::Vector{T}
    stdev::Vector{T}
    x::Union{Array{T},DataFrame}
    fitted::Array{T}
    residuals::Array{T}
    ic::Dict
    stats::Dict
    Φse::Array{T,3}
    pse::Array{T,3}
    Φ0se::Vector{T}
    p0se::Vector{T}
    Φpv::Array{T,3}
    ppv::Array{T,3}
    Φ0pv::Vector{T}
    p0pv::Vector{T}
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
    Φ0pv = xar.Φ0pv
    printstyled("Φ0\n",bold=true,color=:underline)
    pretty_table(string.(round.(Φ0,digits=3)) .* sigf.(Φ0pv), tf = tf_matrix, noheader=true)

    for i in 1:p
        printstyled("Φ",i,"\n",bold=true,color=:underline)
        Φi = xar.Φ[:,:,i]
        Φpvi = xar.Φpv[:,:,i]
        pretty_table(string.(round.(Φi,digits = 3)) .* sigf.(Φpvi), tf = tf_matrix, noheader=true)
    end
    print("Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘^’ 0.1 ‘ ’ 1  and ‘+’ if fixed\n")

    # printstyled("\nCoefficients' Std. Error\n",bold=true,color=:underline)
    # printstyled("\nΦ0se\n",bold=true,color=:underline)
    # pretty_table(xar.Φ0se, tf = tf_matrix, noheader=true)
    # for i in 1:p
    #     printstyled("Φ",i,"se\n",bold=true,color=:underline)
    #     pretty_table(xar.Φse[:,:,i], tf = tf_matrix, noheader=true)    
    # end

    printstyled("\nΣ2 Variance/Covariance Matrix\n",bold=true,color=:underline)
    pretty_table(xar.Σ2, tf = tf_matrix, noheader=true)
    
    printstyled("\nInformation Criteria\n",bold=true,color=:underline)
    pretty_table(DataFrame(xar.ic), nosubheader = true, show_row_number=false)

    printstyled("\nStatistics\n",bold=true,color=:underline)
    pretty_table(DataFrame(xar.stats), nosubheader = true, show_row_number=false)

end

# Return ( num_variables, num_parameters)
function arsize(Φ)
    d = ndims(Φ)
    sΦ = size(Φ)
    if d == 0
        m,np = 1,1
    elseif d == 1
        m,np = 1,sΦ[1]
    elseif d == 2
        m,np = sΦ[1],1
    elseif d == 3
        m,np = sΦ[1],sΦ[3]
    else
        @error "Φ should have less the 4 dimensions"
    end
    (m,np)
end

function sigf(pv::Real)::String

    if pv < 0.001 return  " ***"   end 
    if pv < 0.01  return  " ** "   end 
    if pv < 0.05  return  " *  "   end 
    if pv < 0.1   return  " ^  "   end 
    if pv < 1     return  "    "   end 
    if pv == 1    return  " +  "   end
    
end
