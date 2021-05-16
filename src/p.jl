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

function p(dfx::DataFrame,
           orderlag::AbstractArray = [[0]])

    x = dfx[:,eltype.(eachcol(dfx)) .<: Union{Missing,Real}]
    pnames = "p[" * string(order) * "," * string(lag) * "]_" .* names(x)
    x = Array(x)
    
    px = p(x, orderlag)
    ts = eltype(dfx[:,1]) == Date ? dfx[:,1] : nothing
    dts = ts[2] -ts[1]
    extra = size(px,1) - length(ts)
    ts = ts[1]:dts:ts[end] + extra*dts
    
    pdfx = DataFrame(reshape(px,size(px,1),size(px,2)),pnames)
    
    if !isnothing(ts)
        pdfx = hcat(ts[1:nrow(pdfx)], pdfx)
        rename!(pdfx, names(pdfx)[1] .=> [:Timestamp])
    end
       
    return pdfx

end

function p(fc::FORECAST,
           orderlag::AbstractArray = [[0]])

    n = size(fc.model.x,1)
    m = size(fc.model.x,2)
    
    fx  =  [fc.mean fc.lower fc.upper]
    x   =  repeat(fc.model.x,1,size(fx,2))
    xfx =  [x ; fx]
    
    pxfx = p(xfx, orderlag)

    pfc = deepcopy(fc)
    pfc.model.x = pxfx[1:n,1]
    pfc.mean    = pxfx[n+1:end,1]
    pfc.lower   = pxfx[n+1:end,m+1:3*m]
    pfc.upper   = pxfx[n+1:end,3*m+1:end]
    pfc.call = "Integrated with orderlag" * string(orderlag) * " on " * fc.call

    return(pfc)
    
end



