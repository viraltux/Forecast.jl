"""
Package: Forecast

    function d(x::{AbstractVector,AbstractArray,TimeArray},
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
4609×1 TimeArray{Any,1,Date,Array{Any,1}} 1974-05-17 to 1986-12-31
│            │ A       │
├────────────┼─────────┤
│ 1974-05-17 │ -0.27   │
│ 1974-05-18 │ 0.35    │
│ 1974-05-19 │ 0.18    │
   ⋮
│ 1986-12-30 │ missing │
│ 1986-12-31 │ missing │

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
    b = findlast(!ismissing,x)
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

    if (a-1 > 0) | (b-N > 0)
        return(vcat(Array{Any}(missing,a-1),
                    dxOL,
                    Array{Any}(missing,b-N)))
    else
        return(dxOL)
    end
    
end

function d(x::AbstractArray,
           order::Int=1,
           lag::Int=1;
           center::Bool=false)

    a = findfirst(!ismissing,x)[1]
    b = findlast(!ismissing,x)[1]
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
            xh = hcat(circshift(v,lag),v)
            xh[(1+1):end,M+1:2*M] - xh[(1+1):end,1:M]
        end
        v
    end

    dxOL = (dxO1 ∘ dx1L)(x)

    if center
        Nl = N - size(dxOL)[1]
        dxOL = vcat(dxOL, Array{Any}(missing,Nl,M))
        dxOL = circshift(dxOL,(div(Nl,2),0))
    end
    
    vcat(Array{Any}(missing,a-1,M),
         dxOL,
         Array{Any}(missing,b-N,M))

end

function d(x::TimeArray,
           order::Int=1,
           lag::Int=1;
           center::Bool=false)

    vx = values(x)
    tsx = timestamp(x)
    replace!(vx, NaN => missing)

    dvx = d(vx, order, lag=lag; center=center)
    
    TimeArray(tsx,dvx)

end
