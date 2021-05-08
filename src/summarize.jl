"""
Package: Forecast

    summarize(x; varnames = nothing)

Return statistical summary for x

The values returned are dividen in three sections, the first one shows Minimum, 1st Quantile, Median, Mean, 3rd Quantile, Maxixum and the p-value for the Jarque-Bera Normality Test. The second one show the first four moment; Mean, Variance, Skewness and Kurtosis, an finally a summary with the different types contained in the Array.

# Arguments
- `x`: Array, TimeArray or DataFrame of data.
- `varnames`: Names for the columns to be summarized, it defaults to automatic naming
            or the existing names in TimeArray or DataFrame.

# Returns
A SUMMARIZE struct

# Examples
```julia-repl
julia> summarize(rand(100,3); varnames = ["a","b","c"])
┌──────────┬────────────┬──────────┬──────────┬──────────┬──────────┬──────────┬──────────────┐
│ Variable │        Min │       1Q │   Median │     Mean │       3Q │      Max │ H0 Normality │
├──────────┼────────────┼──────────┼──────────┼──────────┼──────────┼──────────┼──────────────┤
│        a │ 0.00520465 │ 0.205712 │ 0.462199 │ 0.465784 │   0.6913 │  0.97946 │    0.0593599 │
│        b │ 0.00218787 │ 0.247344 │ 0.485465 │ 0.498587 │ 0.723371 │ 0.985226 │    0.0562301 │
│        c │  0.0244256 │ 0.247598 │ 0.530821 │ 0.498689 │ 0.722731 │ 0.967952 │    0.0356495 │
└──────────┴────────────┴──────────┴──────────┴──────────┴──────────┴──────────┴──────────────┘
┌──────────┬──────────┬───────────┬───────────┬───────────┐
│ Variable │     Mean │  Variance │  Skewness │  Kurtosis │
├──────────┼──────────┼───────────┼───────────┼───────────┤
│        a │ 0.465784 │ 0.0823949 │ 0.0823949 │ 0.0823949 │
│        b │ 0.498587 │ 0.0854883 │ 0.0854883 │ 0.0854883 │
│        c │ 0.498689 │ 0.0790597 │ 0.0790597 │ 0.0790597 │
└──────────┴──────────┴───────────┴───────────┴───────────┘
┌──────────┬─────────┐
│ Variable │ Float64 │
├──────────┼─────────┤
│        a │     100 │
│        b │     100 │
│        c │     100 │
└──────────┴─────────┘
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
    sk = ["Variable", "Min","1Q","Median","Mean","3Q","Max","H0 Normality"]
    sv = zeros(m,7)
    s = Array{Any,2}(undef,m,7)

    minimum_Number(x) = minimum(filter(x -> x isa Number, x))
    maximum_Number(x) = maximum(filter(x -> x isa Number, x))
    mean_Number(x) = mean(filter(x -> x isa Number, x))
    quantile_Number(x,q) = quantile(filter(x -> x isa Number, x),q)
    normality(x) = pvalue(JarqueBeraTest(x))

    sv[nd <= 1 ? 1 : 1:end,7] = nd <= 1 ? normality(x) : mapslices(normality,x,dims=[1])
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

