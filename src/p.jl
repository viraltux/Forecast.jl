"""
Package: Forecast

    function p(dx::{AbstractVector,AbstractArray},
               orderlag::{AbstractVector,AbstractArray} = [[0]])

Return reverse lagged differences of a given order for Vector, Array and TimeSeries.

# Arguments
- `dx`:       Vector or Array of data.
- `orderlag`: Initial values of the inverted series. The format per column for the order and lag values is [[lag_1, lag_2, ..., lag_m], order_2, order_3, ..., order_n]. the default values is an initial constant of zero for an inversion of order one and lag one.

# Returns
Lagged differences Vector or Array of a given order.

# Examples
```julia-repl
julia> x = repeat(1:12,3);
julia> dx = d(x,12,12);
julia> orderlag = vcat([collect(1:12)],repeat([0],11));
pjulia> p(dx,orderlag) ≈ x
true
julia> p(d(x),[[1]]) ≈ x
true

julia> tx = hcat(co2().co2, co2().co2, co2().co2)
julia> for col in eachcol(tx) replace!(col, missing => 0.0) end
julia> p(d(tx), [ [[333.38]] [[333.38]] [[333.38]] ] ) ≈ tx
true
julia> p(d(tx,2), [ [[333.38]] [[333.38]] [[333.38]] ; -0.27 -0.27 -0.27]) ≈ tx
true
```
"""
function p(dx::AbstractVector,
           orderlag::AbstractVector = [[0]])

    format_ol = "orderlag format is [[lag_1, lag_2, ..., lag_m], order_2, order_3, ..., order_n]"
    @assert orderlag[1] isa Array format_ol
    @assert orderlag isa Array format_ol
    
    order = length(orderlag)
    lag = length(orderlag[1])
    N = length(dx)

    if (lag == 0)
        return x
    end
    
    a = findfirst(!ismissing,dx)
    b = findlast(!ismissing,dx)
    dx = dx[a:b]

    if !isnothing(findfirst(ismissing, dx))
        @error "Missing values allowed at the start and end of `dx` but not within."
        if !isnothing(findfirst(isnan, dx))
            @error "NaN values not allowed; use `missing` at the start or end of `dx` but not within."
        end
    end

    N = length(dx)
    @assert 0 <= lag & lag <= (N-1)

    function pxO1(v)
        for i in 2:order
            v = cumsum(vcat(orderlag[i],v))                    
        end
        v
    end

    function px1L(v)
        dxm = Array{Union{Missing,Real},2}(missing, lag, Integer(ceil((N+lag)/lag))+1 )
        for i in 1:lag
            cs = cumsum(vcat(orderlag[1][i],circshift(v,-i+1)[1:lag:end]))
            dxm[i,1:length(cs)] = cs
        end
        return vcat(vec(dxm))[1:N+lag+order-1]
    end
    
    pxOL = (px1L ∘ pxO1)(dx)

    if (a-1 > 0) | (b-N > 0)
        return(vcat(Array{Union{Missing,Real}}(missing,a-1),
                    pxOL,
                    Array{Union{Missing,Real}}(missing,b-N)))
    else
        return(pxOL)
    end
    
end

function p(dx::AbstractArray,
           orderlag::AbstractArray = [[0]])

    nc = size(dx)[2]
    if orderlag == [[0]]
        orderlag = hcat(Array{Any,2}(missing, 1, nc))
        orderlag[1,:] .= [[0]]
    end
    
    format_ol = "orderlag format per column is [[lag_1, lag_2, ..., lag_m], order_2, order_3, ..., order_n] in a type order x column Array{Union{Missing,Real},2}"
    @assert orderlag[1] isa Array format_ol
    @assert orderlag isa Array format_ol

    lag = length(orderlag[1,:][1])

    if (lag == 0)
        return dx
    end

    pdx(idx) = p(dx[:,idx], orderlag[:,idx])

    px = pdx(1)
    for i in 2:nc
        px = hcat(px, pdx(i))
    end
    px
    
end

function p(df::DataFrame,
           orderlag::AbstractArray = [[0]])

    df = tots(df)
    ts = df[:,1]
    dfx = df[:,2:end]
    pnames = "p[" * string(orderlag) * "]_" .* names(dfx)
    
    px = p(Array(dfx), orderlag)

    dfpx = DataFrame(reshape(px,size(px,1),size(px,2)),pnames)
    pdfts = DataFrame(vcat(nΔt(ts,-(size(px,1)-size(df,1))), ts), [names(df)[1]])
    
    return hcat(pdfts,dfpx)

end

function p(fc::FORECAST,
           orderlag::AbstractArray = [[0]])

    ts_x = tots(fc.model.x)
    names_x = names(ts_x)
    
    ts_mean = tots(fc.mean)
    ts_lower = tots(fc.lower)
    ts_upper = tots(fc.upper)
    names_mean = names(ts_mean)
    names_lower = names(ts_lower)
    names_upper = names(ts_upper)

    xts = ts_x[:,1]
    x = Array(ts_x[:,2:end])

    fmean = Array(ts_mean[:,2:end])
    flower = Array(ts_lower[:,2:end])
    fupper = Array(ts_upper[:,2:end])

    n = size(x,1)
    m = size(x,2)

    pfxmean = p(vcat(x,fmean), orderlag)
    pfx = convert(Matrix{Float64}, reshape(pfxmean,:,size(x,2))[1:end-size(fmean,1),:])
    pfmean = convert(Matrix{Float64}, reshape(pfxmean,:,size(x,2))[end-size(fmean,1)+1:end,:])

    fts_x = vcat(nΔt(xts,-(size(pfxmean,1)-size(x,1)-size(fmean,1))),xts)

    fts_mean = nΔt(xts,size(pfmean,1))

    pfc = deepcopy(fc)
    pfc.model.x = hcat(DataFrame(fts_x,[names_x[1]]),
                       DataFrame(pfx,names_x[2:end]))
    pfc.mean    = hcat(DataFrame(fts_mean[:,1:1],[names_mean[1]]),
                       DataFrame(pfmean,names_mean[2:end]))

    if size(fmean,2) > 1
        z = fupper .- repeat(fmean,1,size(fmean,2))
        pfmean = repeat(pfmean,1,size(pfmean,2))
    else
        z = fupper .- fmean
    end
    
    pfc.upper = hcat(DataFrame(fts_mean[:,1:1],[names_upper[1]]),
                     DataFrame(pfmean .+ z,names_upper[2:end]))
    pfc.lower = hcat(DataFrame(fts_mean[:,1:1],[names_lower[1]]),
                     DataFrame(pfmean .- z,names_lower[2:end]))
    
    pfc.call = fc.call * "\nIntegrated with orderlag\n" * string(orderlag)

    return(pfc)
    
end
