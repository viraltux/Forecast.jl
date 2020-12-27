"""
function d(x::{AbstractVector,AbstractArray},
           order::Int=1;
           lag::Int=1,
           center::Bool=true,
           padding::Bool=true)

    Return lagged (and iterated) differences for a Vector or an Array

    Args:
        `x`: Vector or Array of data.
        `order`: Order of the differences; number of recursive iterations
                 on the same vector/array.
        `lag`: Lag for the difference.
        `center`: Center the result in the response using Missing values.
        `padding`: Includes or removes `missing` padding.
    Returns:
        Lagged differences Vector or Array of a given order

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

julia> d(x;lag=2,padding=false)
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

julia> d(x,2;lag=2,padding=false)
6×2 Array{Any,2}:
 0  0
 0  0
 0  0
 0  0
 0  0
 0  0
```
"""
function d(x::AbstractVector,
           order::Int=1;
           lag::Int=1,
           center::Bool=true,
           padding::Bool=true)

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
    if padding
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
           padding::Bool=true)

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
        # center values
    if padding
        return circshift(dx,div(lag,2)+div(order,2))
    else
        ivp = findfirst(!ismissing, dx)
        fvp = findlast(!ismissing, dx)
        return(dx[ivp:fvp])
    end

end
