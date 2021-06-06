using Optim

boxcox(x::Vector, λ::Float64)  = λ == 0.0 ? log.(x) : (x .^ λ .- 1.0) ./ λ
iboxcox(x::Vector, λ::Float64) = λ == 0.0 ? exp.(x) : (1.0 .+ λ .* x ) .^ (1/λ)

function boxcox(x::Vector)

    f(λ) = -pvalue(JarqueBeraTest(boxcox(x,λ)))
    λ̂ = Optim.minimizer(optimize(f,-3.0,3.0,Brent()))
    Dict(:x => bc(λ̂), :λ => λ̂)

end




