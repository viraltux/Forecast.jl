"""
Package: Forecast

    function d(x::{AbstractVector,AbstractArray,DataFrame},
               or::Int=1,
               la::Int=1;
               center::Bool=false)

Return Laged differences of a given or for Vector, Array and TimeSeries.

# Arguments
- `x`: Vector or Array of data.
- `or`: Order of the differences; number of recursive iterations on the same vector/array.
- `la`: Lag for the difference.
- `center`: Center the result in the response using Missing values.

# Returns
Laged differences Vector or Array of a given order.

# Examples
```julia-repl
julia> x = [1,2,3,4,5];
julia> d(x)
 d(x)
4-element Vector{Int64}:
 1
 1
 1
 1

julia> d(x,2)
3-element Vector{Int64}:
 0
 0
 0


julia> d(x,1,2)
3-element Vector{Int64}:
 2
 2
 2

julia> x = reshape(collect(1:20),10,2);

julia> d(x,2,2)
6×2 Matrix{Int64}:
 0  0
 0  0
 0  0
 0  0
 0  0
 0  0

julia> d(d(x,1,2),1,2) == d(x,2,2)
```
"""
function d(x::AbstractArray,
           or::Int=1,
           la::Int=1;
           center::Bool=false)

    (la == 0) | (or == 0) && (return x)

    n = size(x,1)
    nv = size(x,2)

    @assert 0 <= la & la <= (n-1) "Lag must be larger or equal to size(x,1)"
    
    dx = x
    for i in or:-1:1
        dx = (dx .- circshift(dx,la))[1+la:end,:]
    end

    if center
        nl = n-size(dx,1)
        dx = vcat(dx, Array{Any}(missing,nl,size(x,2)))
        dx = circshift(dx,(div(nl,2),0))
    end

    return size(dx,2) == 1 ? dx[:,1] : dx

end

function d(df::DataFrame,
           or::Int=1,
           la::Int=1;
           center::Bool=false)

    df = tots(df)
    dfts = df[:,1:1]
    dfx = df[:,2:end]
    dnames = "d[" * string(or) * "," * string(la) * "]_" .* names(dfx)
    
    dx = d(Array(dfx), or, la; center=center)

    dfdx = DataFrame(reshape(dx,size(dx,1),size(dx,2)),dnames)
    ddfts = df[end-size(dx,1)+1:end,1:1]
    
    return hcat(ddfts,dfdx)

end

