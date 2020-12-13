## Loess
using LinearAlgebra
using Statistics


# TODO Check faster implementations in computational methods for local regression
# William S. Cleveland & E. Grosse https://link.springer.com/article/10.1007/BF01890836

function ghat(x::Float64;A,b,d=2,q,ik)
              
    xv = A[:,d]
    yv = b
    ## lambda q
    n = length(xv)
    q = min(q,n)
    xvx = @. abs(xv-x)
    qidx = sortperm(xvx)[1:q]
    qdist = abs(xv[last(qidx)]-x)*max(1,q/n)

    ## upsilon
    w = zeros(n)
    for wi in qidx
        aq = abs(xv[wi]-x)/qdist
        w[wi] = max((1-aq^3)^3,0)
    end
    
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
               predict = xv)
    
    @assert (d==1) | (d==2) "Linear Regression must be of degree 1 or 2"
    @assert length(findall(x -> ismissing(x), xv)) == 0  "xv should not contain missing values"

    myi = findall(x -> !ismissing(x),yv)
    xv = xv[myi]
    yv = yv[myi]

    #ghat.(xv;xv=xv,yv=yv,q=q,d=d)
    res = zeros(length(predict))
    ik = 1.0 ./ k

    ## Ax = b
    A = hcat(xv,repeat([1.0],inner=length(xv)))
    b = yv
    d == 2 ? A = hcat(xv .^ 2.0, A) : nothing

    for (i,xi) in enumerate(predict)
        # TODO find out why expanding this ghat function in this loops causes a
        # drop in performance. For 10,000 values it goes from
        #   6.595928 seconds (210.07 k allocations: 6.534 GiB, 4.96% gc time)
        # to
        #   26.257998 seconds (890.11 M allocations: 20.167 GiB, 7.03% gc time)
        res[i] = ghat(xi;A=A,b=b,d=d,q=q,ik=ik)
    end
    res
end
