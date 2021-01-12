"""
Package: Forecast

    nextodd(x)

Return the smallest odd integer greater than or equal to `x`.        
"""
function nextodd(x::Real)::Integer
    cx = Integer(ceil(x))
    mod(cx,2)==0 ? cx+1 : cx
end

```
Package: Forecast

    drop rows and columns from a matrix
```
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

