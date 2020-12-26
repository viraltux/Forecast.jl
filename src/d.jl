"""
function d(x::AbstractArray,
           order::Int=1;
           lag::Int=1,
           center::Bool=true,
           padding::Bool=true)

    Return lagged (and iterated) differences for an Array

    Args:
        `x`: Array of data.
        `lag`: Lag for the difference.
        `center`: Center the result in the response using Missing values.
        `order`: Order of the differences (number of recursive iterations
                 on the same array).
    Returns:
        Lagged differences Array of a given order

# Examples
"""
function d(x::AbstractVector,
           order::Int=1;
           lag::Int=1,
           center::Bool=true,
           padding::Bool=true)

    N = length(x)
    @assert 0 <= lag & lag <= (N-1)

    if lag == 0
        return x
    end

    dx = Vector{Any}(missing,N)
    
    # using missing to center values
    a = center ? div(lag,2)+1 : 1
    
    if (lag == 1)
        dx[a:(N-lag+a-1)] = diff(x)
    else
        dx[a:(N-lag+a-1)] = diff(hcat(circshift(x,lag),x),dims=2)[(1+lag):end]
    end

    for i in 2:order
        dx = d(dx, order-1; lag=lag, center=center, padding=padding)
    end
    
    dx
end

function d(x::AbstractArray,
           order::Int=1;
           lag::Int=1,
           center::Bool=true,
           padding::Bool=true)

    N,M = size(x)
    @assert 0 <= lag & lag <= (N-1)

    if lag == 0 | order == 0
        return x
    end

    dx = Array{Any}(missing,N,M)

    # using missing to center values
    a = center ? div(lag,2)+1 : 1
    
    if (lag == 1)
        dx[a:(N-lag+a-1),:] = diff(x,dims=1)
    else
        xh = hcat(circshift(x,lag),x)
        dx[a:(N-lag+a-1),:] = xh[(1+lag):end,M+1:2*M] - xh[(1+lag):end,1:M]
    end

    for i in 2:order
        dx = d(dx, order-1; lag=lag, center=center, padding=padding)
    end
    
    dx
end
