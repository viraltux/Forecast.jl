## Loess
using LinearAlgebra
using Statistics

## TODO: Though much faster than before still much slower than Loess.loess
##       Things to look:
##       1- faster lsq?
##       2- KDTree? - simplify lambda for unidimensional series
##       3- pass A instead xv,xy

## TODO: for d=1 it fails and greatly deviates from Loess.loess

function ghat(x::Float64;xv,yv,
              ik=repeat([1.0],inner=length(xv)),
              q,d=2)

    #lambda q
    n = length(xv)
    q = min(q,n)
    xvx = @. abs(xv-x)
    qidx = sortperm(xvx)[1:q]
    qdist = abs(xv[last(qidx)]-x)*max(1,q/n)

    w = zeros(n)
    #upsilon
    for wi in 1:n
        aq = abs(xv[wi]-x)/qdist
        w[wi] = max((1-aq^3)^3,0)
    end

    # Ax = b
    A = hcat(xv,repeat([1.0],inner=length(xv)))
    b = yv
    d == 2 ? A = hcat(xv .^ 2.0, A) : nothing
    
    ## Considering errors are uncorrelated to simplify the calculations with weigths 
    A = @. A*(w*ik)
    b = @. b*(w*ik)

    lsq_x = (A'*A)\(A'*b)
    
    d == 1 ? [x,1.0]'*lsq_x : [x^2.0,x,1.0]'*lsq_x

end

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

    #ghat.(xv;xv=xv,yv=yv,q=q,d=d)
    res = zeros(length(xv))
    ik = 1.0 ./ k
    for (i,xi) in enumerate(xv)
        res[i] = ghat(xi;xv=xv,yv=yv,ik=ik,q=q,d=d)
    end
    res
end
