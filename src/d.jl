"""
Package: Forecast

    function d(x::{AbstractVector,AbstractArray,TimeArray},
               order::Int=1;
               lag::Int=1,
               center::Bool=true,
               pad::Bool=true)

Return Lagged differences of a given order for Vector, Array and TimeSeries.

# Arguments
- `x`: Vector or Array of data.
- `order`: Order of the differences; number of recursive iterations on the same vector/array.
- `lag`: Lag for the difference.
- `center`: Center the result in the response using Missing values.
- `pad`: Includes or removes `missing` pad.

# Returns
Lagged differences Vector or Array of a given order.

# Examples
```julia-repl
julia> x = [1,2,3,4,5];
julia> d(x)
5-element Array{Any,1}:
 1
 1
 1
 1
  missing

julia> d(x,2)
5-element Array{Any,1}:
  missing
 0
 0
 0
  missing

julia> d(x;lag=2,pad=false)
3-element Array{Any,1}:
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

julia> d(x,2;lag=2,pad=false)
6×2 Array{Any,2}:
 0  0
 0  0
 0  0
 0  0
 0  0
 0  0

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
           order::Int=1;
           lag::Int=1,
           center::Bool=true,
           pad::Bool=true)

    N = length(x)
    @assert 0 <= lag & lag <= (N-1)

    if (lag == 0) | (order == 0)
        return x
    end

    dx = Vector{Any}(missing,N)

    # update dx function
    function udx(v)
        if (lag == 1)
            dx[1:(N-lag)] = diff(v)
        else
            dx[1:(N-lag)] = diff(hcat(circshift(v,lag),v),dims=2)[(1+lag):end]
        end
    end

    # update dx with x
    udx(x)

    # update dx with previous dx
    for i in 2:order
        udx(dx)
    end

    
    # center values
    if pad
        return circshift(dx,div(lag,2)+div(order,2))
    else
        ivp = findfirst(!ismissing, dx)
        fvp = findlast(!ismissing, dx)
        return(dx[ivp:fvp])
    end

end


function d(x::AbstractArray,
           order::Int=1;
           lag::Int=1,
           center::Bool=true,
           pad::Bool=true)

    N,M = size(x)
    
    @assert 0 <= lag & lag <= (N-1)

    if (lag == 0) | (order == 0)
        return x
    end

    dx = Array{Any}(missing,N,M)

    # update dx function
    function udx(v)
        if (lag == 1)
            dx[1:(N-lag),:] = diff(v,dims=1)
        else
            xh = hcat(circshift(v,lag),v)
            dx[1:(N-lag),:] = xh[(1+lag):end,M+1:2*M] - xh[(1+lag):end,1:M]
        end
    end

    # update dx with x
    udx(x)

    # update dx with previous dx
    for i in 2:order
        udx(dx)
    end

    # center values
    if pad
        return circshift(dx,div(lag,2)+div(order,2))
    else
        ivp = findfirst(!ismissing, dx)
        fvp = findlast(!ismissing, dx)
        return(dx[ivp:fvp])
    end

end

function d(x::TimeArray,
           order::Int=1;
           lag::Int=1,
           center::Bool=true,
           pad::Bool=true)

    vx = values(x)
    tsx = timestamp(x)
    replace!(vx, NaN => missing)

    dvx = d(vx, order; lag=lag, center=center, pad=pad)
    
    TimeArray(tsx,dvx)

end
