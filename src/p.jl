"""
Package: Forecast

    function p(dx::AbstractVector,
               orderlag::AbstractVector = [[0]];
               center::Bool=false)

Return Inverted lagged differences of a given order for Vector, Array and TimeSeries.

# Arguments
- `dx`:       Vector or Array of data.
- `orderlag`: Initial values of the inverted series. The format for the order and lag values is [[lag_1, lag_2, ..., lag_m], order_2, order_3, ..., order_n]. the default values is an initial constant of zero for an inversion of order one and lag one.
- `center`:   Center the result in the response using Missing values.

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
```
"""
function p(dx::AbstractVector,
           orderlag::AbstractVector = [[0]];
           center::Bool=false)


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

    N = length(dx)
    @assert 0 <= lag & lag <= (N-1)

    function pxO1(v)
        for i in 2:order
            v = cumsum(vcat(orderlag[i],v))                    
        end
        v
    end

    function px1L(v)
        dxm = zeros(lag, Integer(ceil((N+lag)/lag)))
        for i in 1:lag
            cs = cumsum(vcat(orderlag[1][i],circshift(v,-i+1)[1:lag:end]))
            dxm[i,1:length(cs)] = cs
        end
        return vcat(vec(dxm))[1:N+lag+order-1]
    end
    
    pxOL = (px1L ∘ pxO1)(dx)

    if center
        Nl = N-length(pxOL)
        pxOL = vcat(pxOL,Array{Any}(missing,Nl))
        pxOL = circshift(pxOL,div(Nl,2))
    end

    if (a-1 > 0) | (b-N > 0)
        return(vcat(Array{Any}(missing,a-1),
                    pxOL,
                    Array{Any}(missing,b-N)))
    else
        return(pxOL)
    end
    
end
