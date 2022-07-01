"""
Package: Forecast

    ccf(x1::{AbstractVector,DataFrame},
        x2::{AbstractVector,DataFrame};
        type = "cor",
        lag = Integer(ceil(10*log10(length(x1)))),
        alpha = (0.95,0.99))

Compute the cross-correlation or cros-covariance of two univariate series.

The results are normalized to preserve homoscedasticity. The distribution used to normalize the data is an approximation of a Fisher Transformation via a Normal Distribution. There is a plot recipe for the returned object, if the type is `cor` the plot will also show confidence intervals for the given alpha values.

If, for a given integer `k`, `x2` repeats `x1` values such that x1[t] = x2[t+k] for all `i` then high correlation value will be placed *at the right from the center* in the results. That is, this convention will be represented in the plots as `x1_t = x2_{t+k} -> _____0__k__` meaning x2 behavior can be predicted by x1 in k units.

# Arguments
- `x1`: Vector or uni-dimensional DataFrame of data.
- `x2`: Vector or uni-dimensional DataFrame of data.
- `type`: Valid values are "cor" for correlation (default) and "cov" for convariance.
- `lag`: Maximum number of lags.
- `alpha`: A tuple with two thresholds (t1,t2) with t1 <= t2 to plot confidence intervals. The default values are 0.95 and 0.99.
    
# Returns
A CCF object. 

# Examples
```julia-repl
julia> x1 = rand(100);
x2 = circshift(x1,6);
res = ccf(x1, x2; type="cor");
plot(res)
```
"""
function ccf(df1::DataFrame,
             df2::DataFrame;
             type::String = "cor",
             lag::Integer = Integer(ceil(10*log10(length(x1)))),
             alpha::Tuple{Real,Real} = (0.95,0.99))

    x1 = Array(df1[:,eltype.(eachcol(df1)) .<: Real])
    x2 = Array(df1[:,eltype.(eachcol(df1)) .<: Real])

    ccf(x1, x2; type = type, lag = lag, alpha = alpha)

end

function ccf(x1::AbstractVector{<:Real},
             x2::AbstractVector{<:Real};
             type::String = "cor",
             lag::Integer = Integer(ceil(10*log10(length(x1)))),
             alpha::Tuple{Real,Real} = (0.95,0.99))

    N = length(x1)
    @assert N == length(x2) "Vectors should be of equal size"
    @assert N >= 4 "Vectors should have at least four values"
    @assert isnothing(findfirst(isnan,x1)) "No missing values allowed"
    @assert isnothing(findfirst(isnan,x2)) "No missing values allowed"
    @assert isnothing(findfirst(ismissing,x1)) "No missing values allowed"
    @assert isnothing(findfirst(ismissing,x2)) "No missing values allowed"
    @assert type in  ["cor","cov"] "The options for `type` are `cor` and `cov`"
    @assert 1 <= lag <= N-4
    @assert length(alpha) == 2
    @assert 0.0 < alpha[1] < alpha[2] < 1.0

    # auto-ccf
    auto = x1 == x2
    
    ft = (type == "cor" ? Statistics.cor : Statistics.cov)
    
    x = []

    for i in 1:lag
        push!(x,ft(x1[i+1:N],x2[1:(N-i)]))
    end

    if auto
        kx = collect(N-1:-1:N-lag)
    else 
        x = vcat(ft(x1,x2),x) # adding central value (lag=0)
        x = reverse(x)
        for i in 1:lag
            push!(x,ft(x2[i+1:N],x1[1:(N-i)]))
        end
        kx = vcat(collect(N-lag:N-1), N, collect(N-1:-1:N-lag))
    end

    # ccf normalized with the standard error for the
    # Normal distribution of an approximation for
    # the Fisher transformation
    function fse(v::Real)::AbstractFloat
        1/sqrt(v-3)
    end
    k = fse.(kx)/fse(N)
    ccf_res = x ./ k

    call = "ccf("* (auto ? "x" : "x1, x2") *
        "; type=\""*type*
        "\", lag="*string(lag)*
        ", alpha="*string(alpha)*")"

    # Confidence Intervals
    a1 = alpha[1]
    a2 = alpha[2]
    z1 = Distributions.quantile(Normal(), a1 + (1-a1)/2)
    ci1 = z1*fse(N)
    z2 = Distributions.quantile(Normal(), a2 + (1-a2)/2)
    ci2 = z2*fse(N)
    ci = (ci1,ci2)

    CCF(ccf_res, N, type, lag, alpha, ci, auto, call)

end

