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
function drop(M::AbstractMatrix; r=[], c=[])
    s = size(M)
    dr = collect(1:s[1])
    dc = collect(1:s[2])
    splice!(dr,r)
    splice!(dc,c)
    M[dr,dc]
end

"""
Package: Forecast

    compact(x)

Standarize input by dropping empty dimensions and returning either a Number or an Array.
In the case of a DataFrame it removes all non Real columns except if there are columns 
with Date type in which case keeps the first one found and places it as the first column.
"""
function compact(x::Number)
    x
end

function compact(x::AbstractArray)
    x = Array(dropdims(x, dims = tuple(findall(size(x) .== 1)...)))
    ndims(x) == 0 ? x[1] : x 
end

"""
Package: Forecast

    insert_column(M, at, value = 0.0)

Insert a column with specific value at a given position, values are pushed to the right
"""
function insert_column(M::Matrix, at::Integer, value = 0.0)::Matrix
    nr,nc = size(M)
    !(1 <= at <= nc+1) && return(M)

    at == nc+1 && return hcat(M,repeat([value],nr))
    M = hcat(M[:,1:at], M[:,at:nc])
    M[:,at] .= value
    M
end

"""
Package: Forecast

    insert_row(M, at, value = 0.0)

Insert a row with specific value at a given position, values are pushed down
"""
function insert_row(M::Matrix, at::Integer, value = 0.0)::Matrix
    nr,nc = size(M)
    !(1 <= at <= nr+1) && return(M)

    at == nr+1 && return vcat(M,repeat([value],1,nc))
    M = vcat(M[1:at,:], M[at:nr,:])
    M[at,:] .= value
    M
end

"""
Package: Forecast

    insert_row(M, at, value = 0.0)

Insert a row and a column with specific value at a given cross position, 
values are pushed right and down.
"""
function insert_cross(M::Matrix, at::Integer, value = 0.0)::Matrix
    insert_row(insert_column(M,at,value),at,value)
end

