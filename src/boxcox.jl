
"""
Package: Forecast

    boxcox(x::Vector)
    boxcox(x::Vector, λ::Float64)  
    boxcox(x::Vector, λ::Vector)

Compute a Box-Cox power transformation given λ or (λ1,λ2) for data containing negative values, or compute an optimal power transformation if no λ or (λ1,λ2) is provided.

```math
x(\\lambda) =
\\begin{cases}
 \\dfrac{x_i^\\lambda - 1}{\\lambda} & \\text{if } \\lambda \\neq 0 \\\\
 \\ln x_i & \\text{if } \\lambda = 0
\\end{cases}
```
for negative values

```math
x(\\boldsymbol{\\lambda}) =
\\begin{cases}
 \\dfrac{(x_i + \\lambda_2)^{\\lambda_1} - 1}{\\lambda_1} & \\text{if } \\lambda_1 \\neq 0 \\\\
 \\ln (x_i + \\lambda_2) & \\text{if } \\lambda_1 = 0
\\end{cases} 
```

# Arguments
- `x`: Vector to be transformed.
- `λ`: Exponent/s for the tranformation

# Returns
A vector with a boxcox tarnsformation for `x` or a Dict with :x boxcox tranformed and the optimal :λ

# Reference
Box, G. E. P. and Cox, D. R. (1964). An analysis of transformations, Journal of the Royal Statistical Society, Series B, 26, 211-252. A

# Examples
```julia-repl
julia> x = rand(100)
julia> bc = boxcox(x)
julia> iboxcox(bc[:x],bc[:λ]) ≈ x

julia> x = rand(100) .- 0.5
julia> bc = boxcox(x)
julia> iboxcox(bc[:x],bc[:λ]) ≈ x
```
"""
function boxcox(x::Vector{Float64})

    f(λ) = -pvalue(JarqueBeraTest(boxcox(x,λ)))
    
    if isnothing(findfirst(x -> x < 0, x))
        
        λ = Optim.minimizer(optimize(f,-3.0,3.0,Brent()))
    else
        mx = minimum(x)
        lower = [-3.0, abs(mx) + 1.]
        upper = [3.0, Inf]
        initial_x = [0.0, abs(mx) + 1.]
        λ = Optim.optimize(f, lower, upper, initial_x, SAMIN(),
                           Optim.Options(iterations=10^6)).minimizer
    end
    
    return Dict(:x => boxcox(x,λ), :λ => λ)

end

boxcox(x::Vector, λ::Float64)  = λ == 0.0 ? log.(x) : (x .^ λ .- 1.0) ./ λ
boxcox(x::Vector, λ::Vector) = λ[1] == 0.0 ? log.(x .+ λ[2]) : ( (x .+ λ[2]) .^ λ[1] .- 1.0) ./ λ[1]

"""
Package: Forecast

    iboxcox(x::Vector, λ::Float64) 
    iboxcox(x::Vector, λ::Vector) 

Compute the inverse transformation of a Box-Cox power transformation given λ.

# Arguments
- `x`: Vector with a boxcox tranformation to be inverted.
- `λ`: Exponent for the inverse tranformation.

# Returns
A vector with witht the inverse transformation of x given λ.

# Reference
Box, G. E. P. and Cox, D. R. (1964). An analysis of transformations, Journal of the Royal Statistical Society, Series B, 26, 211-252. A

# Examples
```julia-repl
julia> x = rand(100)
julia> bc = boxcox(x)
julia> iboxcox(bc[:x],bc[:λ]) ≈ x

julia> x = rand(100) .- 0.5
julia> bc = boxcox(x)
julia> iboxcox(bc[:x],bc[:λ]) ≈ x
```
"""
iboxcox(x::Vector, λ::Float64) = λ == 0.0 ? exp.(x) : (1.0 .+ λ .* x ) .^ (1/λ)
iboxcox(x::Vector, λ::Vector) = λ[1] == 0.0 ? exp.(x) .- λ[2] : (1.0 .+ λ[1] .* x ) .^ (1/λ[1]) .- λ[2]
