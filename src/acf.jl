"""
Package: Forecast

    acf(x::{AbstractVector,DataFrame},
        type = "cor",
        lag = Integer(ceil(10*log10(length(x)))),
        alpha = (0.95,0.99))

Compute the auto-correlation or auto-covariance for an univariate series.

The results are normalized to preserve homoscedasticity. The distribution used to normalize the data is an approximation of a Fisher Transformation via a Normal Distribution. There is a plot recipe for the returned object, if the type is `cor` the plot will also show confidence intervals for the given alpha values.

# Arguments
- `x`: Vector or uni-dimensional DataFrame of data.
- `type`: Valid values are "cor" for correlation (default) and "cov" for convariance.
- `lag`: Maximum number of lags.
- `alpha`: A tuple with two thresholds (t1,t2) with t1 <= t2 to plot confidence intervals. The default values are 0.95 and 0.99.

# Returns
A CCF object.

# Examples
```julia-repl
julia> x = rand(100);
res = acf(x; type="cor");
plot(res)
```
"""
function acf(df::DataFrame;
             type::String = "cor",
             lag::Integer = Integer(ceil(10*log10(length(x)))),
             alpha::Tuple{Real,Real} = (0.95,0.99))

    x = Array(df[:,eltype.(eachcol(df)) .<: Real])
    ccf(x, x; type, lag, alpha)
end

function acf(x::AbstractVector{<:Real};
             type::String = "cor",
             lag::Integer = Integer(ceil(10*log10(length(x)))),
             alpha::Tuple{Real,Real} = (0.95,0.99))

    ccf(x,x; type, lag, alpha)
end

