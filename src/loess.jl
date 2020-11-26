## Loess
using Statistics
using Plots
using LinearAlgebra
using DataFrames

# Tricube weight function as described by Cleveland et al.
function Weights(u)::Float64
    @assert u >= 0
    w = (1-u^3)^3
    max(0,w)
end

# qth farthest xi from x
function lambda(x,q,xv)
    isanx = isa(x,Number) 
    isanx ? x = [x] : 0
    x = convert(Array{Float64,1},x)
    q = convert(Int64,q)
    xv = convert(Array{Float64},xv)
    @assert (sort(xv) == xv)
    res = Float64[]
    push!(res,1)
    for xi in x
        xvxi = abs.(xv .- xi)
        n = length(xv)
        k = max(1,q/n)
        q = min(q,n)
        p = sortperm(xvxi)[q] 
        push!(res,abs(xv[p]-xi)*k)
    end
    isanx ? res[1] : res
end

# Neighborhood weight for any xv[i]
# i the index for the focal point
upsilon(x,i,xv,q) = Weights.(abs.(xv[i].-x) / lambda(x,q,xv))

# Least Square Fitting adjusted by Weights and Variance
function lsq(x,y;
             d=2,
             w = repeat([1.0],inner=length(x)),
             v = repeat([1.0],inner=length(x)))

    @assert (d==1) | (d==2)

    A = hcat(x,repeat([1],inner=length(x)))
    b = y
    d == 2 ? A = hcat(x .^ 2, A) : nothing
    # W = Diagonal(w .* v)
    ## (A'*W*A)^-1*A'*W*b #numerical issues
    ## (A'*W*A)\A'*W*b # uses QR but still numerical issues
    ## Considering errors are uncorrelated to simplify the calculations with weigths 
    A = A.*w
    b = b.*w

    ## TODO: Simplification not needed since A will be dense, consider removal.
    A = A[w.!=0,:]
    b = b[w.!=0]
    #(A'*A)\A'*b
    #(A'*A)^-1*A'*b
    #IterativeSolvers.lsmr(A,b)
    pinv(A'*A,rtol=sqrt(eps(real(float(one(eltype(A)))))))*A'*b
end

function loess(xv,yv;
               d=2,
               x = xv,
               v=repeat([1.0],inner=length(xv)),  
               q=round(3/4*length(xv)),
               iter = 3)

    data = dropmissing(DataFrame(hcat(xv,yv)))
    xv = (data.x1.-minimum(data.x1))/(maximum(data.x1)-minimum(data.x1))
    yv = data.x2
    n = length(xv)
    
    res = Float64[]
    qw = min(length(xv),q)
    for loit in 1:iter
        res = Float64[]
        for fp in 1:n
            #j focal point
            # (i,k) is the lsq window and j the focal point within the window
            j = fp
            i = convert(Int64,min(n-(qw-1),max(1,j-floor((qw-1)/2))))
            k = convert(Int64,min(n,i+qw-1))
            wik = upsilon(xv[j],i:k,xv,q)
            reg_param = lsq(xv[i:k],yv[i:k];d=d,w=wik,v=v[i:k])
            A = hcat(xv[i:k],repeat([1],inner=length(i:k)))
            d == 2 ? A = hcat(xv[i:k].^2,A) : nothing
            reg_x = A*reg_param
            push!(res, reg_x[j-i+1])
        end
        yv = res
    end

    DataFrame(x = data.x1, y = res)
end
