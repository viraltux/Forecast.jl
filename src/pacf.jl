"""
Package: Forecast

    pacf(x::{AbstractVector,TimeArray},
         lag = Integer(ceil(10*log10(length(x1)))),
         alpha = (0.95,0.99))

Compute the partial auto-correlation for an univariate series.

The results are normalized to preserve homoscedasticity. The distribution used to normalize the data is an approximation of a Fisher Transformation via a Normal Distribution. There is a plot recipe for the returned object, if the type is `cor` the plot will also show confidence intervals for the given alpha values.

# Arguments
- `x`: Vector or uni-dimensional TimeArray of data.
- `lag`: Maximum number of lags.
- `alpha`: A tuple with two thresholds (t1,t2) with t1 <= t2 to plot confidence intervals. The default values are 0.95 and 0.99.

# Returns
A CCF object

# Examples
```julia-repl
julia> x = rand(100);
res = pacf(x1);
plot(res)
```
"""
function pacf(ta::TimeArray;
             lag = Integer(ceil(10*log10(length(x1)))),
             alpha = (0.95,0.99))

    ccf(values(ta), values(ta); type = type, lag = lag, alpha = alpha)
end

function pacf(x::AbstractVector;
             lag = Integer(ceil(10*log10(length(x)))),
             alpha = (0.95,0.99))

    ccf(x,x; type = type, lag = lag, alpha = alpha)
end




function pacf_step(x::AbstractVector;
                   lag = Integer(ceil(10*log10(length(x)))),
                   alpha = (0.95,0.99))

    pacf_res = [1.0]
    lambda = 10^-36
    N = length(x)
    M = zeros(N-lag,lag+1)
    for i in 1:(N-lag)
        M[i,:] = x[i:i+lag]
    end

    for i in 1:lag
        b = M[:,1:i]*pacf_res
        A = M[:,i+1]
        A = hcat(A,repeat([1.0],N-lag))
        lsq = (A'*A+lambda*I(2))\(A'*b)
        push!(pacf_res,-lsq[end-1]) 
    end

    -pacf_res[2:end]

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

function pacf_real(x::AbstractVector;
             lag = Integer(ceil(10*log10(length(x)))),
             alpha = (0.95,0.99))

    
    # prevent singularities
    lambda = 10^-36

    N = length(x)
    A = zeros(N-lag,lag+1)
    for i in 1:(N-lag)
        A[i,:] = x[i:i+lag]
    end
    bt1 = A[:,1]
    A = A[:,2:end]
    A = hcat(A,repeat([1.0],N-lag))

    pacf_res = []
    for i in 2:lag+1
        A1i = drop(A, c=i-1)
        lsq = (A1i'*A1i+lambda*I(lag))\(A1i'*bt1)
        r1 = round.(bt1.-A1i*lsq, digits=10)

        bt0 = A[:,i-1]
        lsq = (A1i'*A1i+lambda*I(lag))\(A1i'*bt0)
        r0 = round.(bt0.-A1i*lsq, digits=10)
        push!(pacf_res,cor(r0.-mean(r0),r1.-mean(r1)))
        
    end
    
    pacf_res
end
