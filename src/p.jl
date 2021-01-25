"""
Package: Forecast

    function p(dx::{AbstractVector,AbstractArray,TimeArray},
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
julia> p(dx,orderlag) ≈ x
true
julia> p(d(x),[[1]]) ≈ x
true

julia> tx = hcat(co2(), co2(), co2());
julia> vtx = values(tx);
julia> tx = TimeArray(timestamp(tx),replace(vtx, missing => 0.0));
julia> values(p(d(tx), [ [[333.38]] [[333.38]] [[333.38]] ] )) ≈ values(tx)
true
julia> values(p(d(tx,2), [ [[333.38]] [[333.38]] [[333.38]] ; -0.27 -0.27 -0.27] )) ≈ values(tx)
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
        dxm = Array{Any,2}(missing, lag, Integer(ceil((N+lag)/lag))+1 )
        for i in 1:lag
            cs = cumsum(vcat(orderlag[1][i],circshift(v,-i+1)[1:lag:end]))
            dxm[i,1:length(cs)] = cs
        end
        return vcat(vec(dxm))[1:N+lag+order-1]
    end
    
    pxOL = (px1L ∘ pxO1)(dx)

    if (a-1 > 0) | (b-N > 0)
        return(vcat(Array{Any}(missing,a-1),
                    pxOL,
                    Array{Any}(missing,b-N)))
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
    
    format_ol = "orderlag format per column is [[lag_1, lag_2, ..., lag_m], order_2, order_3, ..., order_n] in a type order x column Array{Any,2}"
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

function p(dx::TimeArray,
           orderlag::AbstractArray = [[0]])

    vx = values(dx)
    pvx = p(vx, orderlag)

    tsx = timestamp(dx)
    dt = tsx[2]-tsx[1]

    tsx = tsx[1]:dt:tsx[1]+(size(pvx)[1]-1)*dt

    TimeArray(tsx,pvx)
end
