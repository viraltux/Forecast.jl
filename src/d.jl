"""
Package: Forecast

    function d(x::{AbstractVector,AbstractArray,DataFrame},
               order::Int=1,
               lag::Int=1;
               center::Bool=false)

Return Lagged differences of a given order for Vector, Array and TimeSeries.

# Arguments
- `x`: Vector or Array of data.
- `order`: Order of the differences; number of recursive iterations on the same vector/array.
- `lag`: Lag for the difference.
- `center`: Center the result in the response using Missing values.

# Returns
Lagged differences Vector or Array of a given order.

# Examples
```julia-repl
julia> x = [1,2,3,4,5];
julia> d(x)
4-element Array{Int64,1}:
 1
 1
 1
 1

julia> d(x,2)
3-element Array{Int64,1}:
 0
 0
 0

julia> d(x,1,2)
3-element Array{Int64,1}:
 2
 2
 2

julia> x = reshape(collect(1:20),10,2)
10×2 Array{Int64,2}:
  1  11
  2  12
  3  13
  4  14
  5  15
  6  16
  7  17
  8  18
  9  19
 10  20 

julia> d(x,2,2)
8×2 Array{Any,2}:
 2  2
 2  2
 2  2
 2  2
 2  2
 2  2
 2  2
 2  2

julia> d(co2())
[ Info: Dataset used in Cleveland et al. paper
4608×1 TimeArray{Any,1,Date,Array{Any,1}} 1974-05-17 to 1986-12-30
│            │ A       │
├────────────┼─────────┤
│ 1974-05-17 │ -0.27   │
│ 1974-05-18 │ 0.35    │
│ 1974-05-19 │ 0.18    │
   ⋮
│ 1986-12-28 │ 0.04    │
│ 1986-12-29 │ -0.12   │
│ 1986-12-30 │ missing │
```
"""
function d(x::AbstractVector,
           order::Int=1,
           lag::Int=1;
           center::Bool=false)

    if (lag == 0) | (order == 0)
        return x
    end

    a = findfirst(!ismissing,x)
    ma = a-1
    b = findlast(!ismissing,x)
    mb = length(x)-b
    x = x[a:b]

    N = length(x)
    @assert 0 <= lag & lag <= (N-1)

    function dx1L(v)
        diff(hcat(circshift(v,lag),v),dims=2)[lag+1:end]
    end

    function dxO1(v)
        for i in 2:order
            v = diff(v)
        end
        v
    end

    dxOL = (dxO1 ∘ dx1L)(x)

    if center
        Nl = N-length(dxOL)
        dxOL = vcat(dxOL,Array{Any}(missing,Nl))
        dxOL = circshift(dxOL,div(Nl,2))
    end

    if (ma > 0) | (mb > 0)
        return(vcat(Array{Any}(missing,ma),
                    dxOL,
                    Array{Any}(missing,mb)))
    else
        return(dxOL)
    end
    
end

function d(x::AbstractArray,
           order::Int=1,
           lag::Int=1;
           center::Bool=false)

    a = findfirst(!ismissing,x)[1]
    ma = a-1
    b = findlast(!ismissing,x)[1]
    mb = size(x)[1]-b
    x = x[a:b,:]

    N,M = size(x)
    @assert 0 <= lag & lag <= (N-1)

    if (lag == 0) | (order == 0)
        return x
    end

    function dx1L(v)
        xh = hcat(circshift(v,lag),v)
        xh[(1+lag):end,M+1:2*M] - xh[(1+lag):end,1:M]
    end

    function dxO1(v)
        for i in 2:order
            v = diff(v, dims = 1)
        end
        v
    end

    dxOL = (dxO1 ∘ dx1L)(x)

    if center
        Nl = N - size(dxOL)[1]
        dxOL = vcat(dxOL, Array{Any}(missing,Nl,M))
        dxOL = circshift(dxOL,(div(Nl,2),0))
    end
    
    if (ma > 0) | (mb > 0)
        vcat(Array{Any}(missing,ma,M),
             dxOL,
             Array{Any}(missing,mb,M))
    else
        return(dxOL)
    end

end

function d(dfx::DataFrame,
           order::Int=1,
           lag::Int=1;
           center::Bool=false)

    x = dfx[:,eltype.(eachcol(dfx)) .<: Union{Missing,Real}]
    dnames = "d[" * string(order) * "," * string(lag) * "]_" .* names(x)
    x = Array(x)
    
    dx = d(x, order, lag; center=center)
    timestamp = eltype(dfx[:,1]) == Date ? dfx[:,1] : nothing

    ddfx = DataFrame(reshape(dx,size(dx,1),size(dx,2)),dnames)
    
    if !isnothing(timestamp)
        ddfx = hcat(timestamp[1:nrow(ddfx)], ddfx)
        rename!(ddfx, names(ddfx)[1] .=> [:Timestamp])
    end
       
    return ddfx
end

