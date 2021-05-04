"""
Package: Forecast

    nextodd(x)

Return the smallest odd integer greater than or equal to `x`.        
"""
function nextodd(x::Real)::Integer
    cx = Integer(ceil(x))
    mod(cx,2)==0 ? cx+1 : cx
end

"""
Package: Forecast

    drop(M;r,c)

Drop rows and columns from a matrix.
"""
function drop(M::AbstractMatrix;
              r=nothing,
              c=nothing)
    s = size(M)
    dr = collect(1:s[1])
    dc = collect(1:s[2])
    isnothing(r) ? nothing : splice!(dr,r)
    isnothing(c) ? nothing : splice!(dc,c)
    M[dr,dc]
end

"""
Package: Forecast

    compact(x)

Standarize input by dropping empty dimensions and returning either a Number or an Array.
"""
function compact(x)
    if x isa Number return x end
    x = Array(x)
    sx = size(x)
    x = reshape(x,sx)
    x = dropdims(x, dims = tuple(findall(size(x) .== 1)...))
    ndims(x) == 0 ? x[1] : x 
end


