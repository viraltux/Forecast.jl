## Loess
using LinearAlgebra
using Statistics
## TODO: consider returning an Array instead a DataFrame since Arrays
## are default types

# Least Square Fitting adjusted by Weights and Variance Factors
function lsq(x,y;
             d=2,
             w = repeat([1.0],inner=length(x)),
             k = repeat([1.0],inner=length(x)))

    @assert (d==1) | (d==2)
    @assert all(w.>=0) & all(k.>0) "weights must be positive or zero and variance factors positive"

    ## Removing zero weighted values from to speed up calculations
    x = x[w.!=0]
    y = y[w.!=0]
    k = k[w.!=0]
    w = w[w.!=0]
    
    # Ax = b
    A = hcat(x,repeat([1],inner=length(x)))
    b = y
    d == 2 ? A = hcat(x .^ 2, A) : nothing
    
    ## Considering errors are uncorrelated to simplify the calculations with weigths 
    A = A .* (w .* (1 ./ k))
    b = b .* (w .* (1 ./ k))

    # Removing pinv since it fails for long series > 10,000 and introduces
    # artifacts (jumps) in the smoothing.
    # Using pinv to avoid near singular errors with recommended tolerance
    # rtol=sqrt(eps(real(float(one(eltype(A))))))
    # pinv(A'*A,rtol=rtol)*A'*b
    (A'*A)\(A'*b)
end

# using DataFrames
# lm(@formula(y~x),DataFrame(x=x,y=y))

function weights(u::Float64;type::String="tricube")::Float64
    @assert type in ["tricube"]
    @assert u >= 0 "values must be positive or zero"
    if type == "tricube"
        # Tricube weight function as described by Cleveland et al.
        w = (1-u^3)^3
        return max(0,w)
    end
end

# distance of qth farthest xi from x
function lambda(x::Float64,q::Int64;xv::Array{Float64})
    @assert (sort(xv) == xv)
    n = length(xv)
    q = min(q,n)
    xvx = abs.(xv .- x)
    qidx = sortperm(xvx)[1:q]
    qdist = abs(xv[last(qidx)]-x)*max(1,q/n)
    (qdist = qdist, qidx = qidx)
end

# Neighborhood weight for any xi 
function upsilon(x::Float64,i::Int64;q::Int64,xv::Array{Float64})
    lambda_qx, qidx = lambda(x,q;xv=xv)
    i in qidx ? weights(abs(xv[i]-x)/lambda_qx) : 0.0
end

function ghat(x::Float64;xv,yv,
              k=repeat([1],inner=length(xv)),
              q,d=2)

    n = length(xv)
    w = upsilon.(x,1:n,q=q,xv=xv)
    #(@isdefined k) ? nothing :
    lsq_x = lsq(xv,yv;d=d,w=w,k=k)
    d == 1 ? [x,1]'*lsq_x : [x^2,x,1]'*lsq_x

end

using NearestNeighbors

function loess(xv,yv;
               d=2,
               k=repeat([1.0],inner=length(xv)),  
               q=Int64(round(3/4*length(xv))),
               iter = 3,
               model = false)
    
    @assert (d==1) | (d==2) "Linear Regression must be of degree 1 or 2"
    @assert length(findall(x -> ismissing(x), xv)) == 0  "xv should not contain missing values"

    myi = findall(x -> !ismissing(x),yv)
    xv = xv[myi]
    yv = yv[myi]

    ghat.(xv;xv=xv,yv=yv,q=q,d=d)
end
