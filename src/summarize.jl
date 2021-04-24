"""
Package: Forecast

    summarize(x; labels = nothing)

Return statistical summary for x

The values returned are Minimum, 1st Quantile, Median, Mean, 3rd Quantile and, Maxixum

# Arguments
- `x`: Vector of data.
- `labels`: 

# Returns
Vector of moving average smoothed values containing `missing` values to preserve the size of the original vector.

# Examples
```julia-repl
julia> sma(1:5,3;center=true)
5-element Array{Any,1}:
  missing
 2.0
 3.0
 4.0
  missing
```
"""

function summarize(x::TimeArray; labels = nothing)
    ta = x
    x = values(ta)
    labels = string.(colnames(ta))
    summarize(x, labels = labels)
end

function summarize(x::Vector; labels = nothing)
    summarize(reshape(x,length(x),1), labels = labels)
end

function summarize(x::Number; labels = nothing)
    summarize([x,x], labels = labels)
end

function summarize(x::AbstractArray; labels = nothing)
    
    x = compact(x)
    nd = ndims(x)
    n,m = nd == 1 ? (length(x),1) : size(x)
    sk = ["Min","1Q","Median","Mean","3Q","Max"]
    sv = zeros(m,6)
    s = Array{Any,2}(undef,m+1,6)
    
    sv[nd <= 1 ? 1 : 1:end,1] = nd <= 1 ? minimum(x) : mapslices(minimum,x,dims=[1])
    sv[nd <= 1 ? 1 : 1:end,6] = nd <= 1 ? maximum(x) : mapslices(maximum,x,dims=[1])
    sv[nd <= 1 ? 1 : 1:end,4] = nd <= 1 ? mean(x) : mapslices(mean,x,dims=[1])
    sv[nd <= 1 ? 1 : 1:end,[2,3,5]] = nd <= 1 ? quantile(x,[.25,.5,.75]) :
                             mapslices(x -> quantile(x,[.25,.5,.75]),x,dims=[1])'
    s[1,:] = sk
    s[2:end,:] = round.(sv,digits=5)

    vlabels = isnothing(labels) ? ["d"*string(i) for i in 1:m] : labels
    @assert length(vlabels) + 1 == size(s,1)[1]
    
    hcat(vcat(["🅂"],vlabels),s)
end
