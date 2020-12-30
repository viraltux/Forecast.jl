using Plots
using Statistics
using Distributions

"""
function ccf(x1::AbstractVector,
             x2::AbstractVector;
             type = "cor",
             lag = Integer(10*log10(length(x1))),
             plot::Bool = false,
             alpha = (0.95,0.99))

Compute the cross-correlation or cross-covariance of two univariate series.

    Args:
        `x1`: Vector of data.
        `x2`: Vector of data.
        `type`: Valid values are "cor" for correlation (default) and "cov" for convariance.
        `plot`: Plots the ccf with a confidence interval at 0.95 (alpha)
        `alpha`: Two thresholds to plot confidence intervals
    Returns:
        Vector of cross-correlation or covariance between two vectors 
        plus an optional plot with cofidence intervals

# Examples
```julia-repl
julia> x1 = rand(100);
x2 = circshift(x1,5);
ccf(x1, x2; type="cor", plot=true);

```
"""
function ccf(x1::AbstractVector,
             x2::AbstractVector;
             type = "cor",
             lag = Integer(10*log10(length(x1))),
             plot::Bool = false,
             alpha = (0.95,0.99))

    N = length(x1)
    @assert N == length(x2) "Vectors should be of equal size"
    @assert N >= 4 "Vectors should have at least four values"
    @assert isnothing(findfirst(isnan,x1)) "No missing values allowed"
    @assert isnothing(findfirst(isnan,x2)) "No missing values allowed"
    @assert type in  ["cor","cov"] "The options for `type` are `cor` and `cov`"
    @assert 1 <= lag <= N-4
    @assert length(alpha) == 2
    @assert 0.0 < alpha[1] <= alpha[2] < 1.0

    ft = (type == "cor" ? cor : cov)
    
    x = []
    for i in 0:lag
        push!(x,ft(x1[i+1:N],x2[1:(N-i)]))
    end

    x = reverse(x)
    for i in 1:lag
        push!(x,ft(x2[i+1:N],x1[1:(N-i)]))
    end
    
    kx = vcat(collect(N-lag:N-1), N, collect(N-1:-1:N-lag))
    # ccf normalized with the standard error for the
    # Normal distsribution of an approximation for
    # the Fisher transformation
    function fse(v)::AbstractFloat
        1/sqrt(v-3)
    end
    k = fse.(kx)/fse(N)
    ccf_res = x ./ k

    if plot
        z0 = Distributions.quantile(Normal(), alpha[1] + (1-alpha[1])/2)
        z1 = Distributions.quantile(Normal(), alpha[2] + (1-alpha[2])/2)
        ps = Plots.sticks(ccf_res,
                          linewidth = 100/(lag+1),
                          legend = nothing,
                          xticks = (1:div(lag,5):2*lag+1,
                                    -lag:div(lag,5):lag))
        ci0 = z0*fse(N)
        ci1 = z1*fse(N)
        Plots.hline!([-ci0,ci0], color = :green, linealpha = 0.5, linestyle = :dash)
        Plots.hline!([-ci1,ci1], color = :orange, linealpha = 0.5, linestyle = :dot)
        display(ps)
    end
    ccf_res
end
