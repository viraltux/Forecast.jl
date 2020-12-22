# TODO Check faster implementations in computational methods for local regression
# William S. Cleveland & E. Grosse https://link.springer.com/article/10.1007/BF01890836
function ghat(x::Float64;A,b,d=2,q,rho)
              
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
    
    A = @. A*(w*rho)
    b = @. b*(w*rho)
    
    lsq_x = (A'*A)\(A'*b)

    d == 1 ? [x,1.0]'*lsq_x : [x^2.0,x,1.0]'*lsq_x

end


"""
function loess(xv,yv;
               d=2,
               q=Int64(round(3/4*length(xv))),
               rho=repeat([1.0],inner=length(xv)),  
               predict = xv)

Smooth a vector of observations using locally weighted regressions.

Although loess can be used to smooth observations for any given number of independent variables, this implementation is univariate. The speed of loess can be greatly increased by using fast aproximations for the linear fitting calculations, however this implementation calculates only exact results.

The loess functionality and nomenclature follows the descriptions in:

    "STL: A Seasonal, Trend Decomposition Procedure Based on Loess"
     Robert B. Cleveland, William S. Cleveland, Jean E. McRae, and Irma Terpenning.
     Journal of Official Statistics Vol. 6. No. 1, 1990, pp. 3-73 (c) Statistics Sweden.

    Args:
        `xv`: Observations' support.
        `yv`: Observation values.
        `d`: Degree of the linear fit, it accepts values 1 or 2.
        `q`: As q increases loess becomes smoother, when q tends to infinity loess tends to an ordinary least square poynomial fit of degree `d`. It defaults to the rounding of 3/4 of xv's length.
        `rho`: Weights expressing the reliability of the observations (e.g. if yi had variances sigma^2*ki where ki where known, the rhoi could be 1/ki). It defaults to 1.0.
        `predict`: Vector containing the real values to be predicted, by default predicts xv.
    Returns:
        The loess values for the values contained in `predict`.

# Examples
```julia-repl
julia> loess(rand(5),rand(5))
5-element Array{Float64,1}:
 0.8279856470091445
 0.6498793177642597
 0.017078421335739336
 0.8595786774014016
 0.2564676547690037
```
"""
function loess(xv,yv;
               d=2,
               q=Int64(round(3/4*length(xv))),
               rho=repeat([1.0],inner=length(xv)),  
               predict = xv)
    
    @assert (d==1) | (d==2) "Linear Regression must be of degree 1 or 2"
    @assert length(findall(x -> ismissing(x), xv)) == 0  "xv should not contain missing values"

    myi = findall(x -> !ismissing(x),yv)
    xv = xv[myi]
    yv = yv[myi]
    rho = rho[myi]
    
    res = zeros(length(predict))

    ## Ax = b
    A = hcat(xv,repeat([1.0],inner=length(xv)))
    b = yv
    d == 2 ? A = hcat(xv .^ 2.0, A) : nothing

    for (i,xi) in enumerate(predict)
        res[i] = ghat(xi;A=A,b=b,d=d,q=q,rho=rho)
    end
    res
end
