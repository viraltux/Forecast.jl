"""
Package: Forecast

    summarize(x; varnames = nothing)

Return statistical summary for x

The values returned are Minimum, 1st Quantile, Median, Mean, 3rd Quantile and Maxixum.

# Arguments
- `x`: Vector of data.
- `varnames`: Names for the columns to be summarized, it default to automatic numbering
            or the existing names in TimeArray or DataFrame

# Returns
A SUMMARIZE struct

# Examples
```julia-repl
julia> summarize(rand(100,3))
3×7 DataFrame
 Column  Min         1Q        Median    Mean      3Q        Max      
──────────────────────────────────────────────────────────────────────
 x1      0.00339372  0.364233  0.62742   0.567196  0.782794  0.999332
 x2      0.00981478  0.251394  0.471664  0.495316  0.776891  0.989918
 x3      0.0137438   0.213655  0.461344  0.49218   0.774042  0.98995
```
"""
function summarize(x::TimeArray; names = nothing)
    ta = x
    x = values(ta)
    varnames = string.(colnames(ta))
    summarize(x, varnames = varnames)
end

function summarize(x::DataFrame; varnames = nothing)
    summarize(convert(Array,x), varnames = names(x))
end

function summarize(x::Vector; varnames = nothing)
    summarize(reshape(x,length(x),1), varnames = varnames)
end

function summarize(x::Number; varnames = nothing)
    summarize([x,x], varnames = varnames)
end

function summarize(x::AbstractArray; varnames = nothing)

    @assert ndims(x) <= 2 "Data should have two dimensions at most."
    x = compact(x)
    nd = ndims(x)
    n,m = nd == 1 ? (length(x),1) : size(x)

    varnames = isnothing(varnames) ? ["x"*string(i) for i in 1:m] : varnames

    # Quantiles and Mean
    sk = ["Variable", "Min","1Q","Median","Mean","3Q","Max"]
    sv = zeros(m,6)
    s = Array{Any,2}(undef,m,6)

    minimum_Number(x) = minimum(filter(x -> x isa Number, x))
    maximum_Number(x) = maximum(filter(x -> x isa Number, x))
    mean_Number(x) = mean(filter(x -> x isa Number, x))
    quantile_Number(x,q) = quantile(filter(x -> x isa Number, x),q)
    
    sv[nd <= 1 ? 1 : 1:end,1] = nd <= 1 ? minimum(x) : mapslices(minimum_Number,x,dims=[1])
    sv[nd <= 1 ? 1 : 1:end,6] = nd <= 1 ? maximum(x) : mapslices(maximum_Number,x,dims=[1])
    sv[nd <= 1 ? 1 : 1:end,4] = nd <= 1 ? mean(x) : mapslices(mean_Number,x,dims=[1])
    sv[nd <= 1 ? 1 : 1:end,[2,3,5]] = nd <= 1 ? quantile_Number(x,[.25,.5,.75]) :
                             mapslices(x -> quantile_Number(x,[.25,.5,.75]),x,dims=[1])'
    qua = DataFrame(hcat(varnames,sv), :auto)
    DataFrames.rename!(qua, sk)

    # Moments
    sk = ["Variable", "Mean", "Variance", "Skewness", "Kurtosis"]
    sv = zeros(m,4)
    s = Array{Any,2}(undef,m,4)

    var_Number(x) = var(filter(x -> x isa Number, x))
    skewness_Number(x) = var(filter(x -> x isa Number, x))
    kurtosis_Number(x) = var(filter(x -> x isa Number, x))

    sv[nd <= 1 ? 1 : 1:end,1] = nd <= 1 ? mean(x) : mapslices(mean_Number,x,dims=[1])
    sv[nd <= 1 ? 1 : 1:end,2] = nd <= 1 ? var(x) : mapslices(var_Number,x,dims=[1])
    sv[nd <= 1 ? 1 : 1:end,3] = nd <= 1 ? skewness(x) : mapslices(skewness_Number,x,dims=[1])
    sv[nd <= 1 ? 1 : 1:end,4] = nd <= 1 ? kurtosis(x) : mapslices(kurtosis_Number,x,dims=[1])

    mo = DataFrame(hcat(varnames,sv), :auto)
    DataFrames.rename!(mo, sk)

    # Format Check
    tox = string.(typeof.(x))
    tx = sort(unique(tox))
    nt = length(tx)
    mv = length(varnames)
    vt = Array{Any}(undef,mv,nt)
    vt[:] .= 0

    for (i,ty) in enumerate(tx)
        vt[:,i] = count(==(ty), tox, dims=1) 
    end

    tdf = DataFrame(hcat(varnames,vt),:auto)
    fc = DataFrames.rename!(tdf,vcat("Variable",  tx))
    
    SUMMARIZE(qua,mo,fc)

end

