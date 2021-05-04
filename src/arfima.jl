"""
The [Multinomial distribution](http://en.wikipedia.org/wiki/Multinomial_distribution)
generalizes the *binomial distribution*. Consider n independent draws from a Categorical
distribution over a finite set of size k, and let ``X = (X_1, ..., X_k)`` where ``X_i``
represents the number of times the element ``i`` occurs, then the distribution of ``X``
is a multinomial distribution. Each sample of a multinomial distribution is a k-dimensional
integer vector that sums to n.
The probability mass function is given by
```math
f(x; n, p) = \\frac{n!}{x_1! \\cdots x_k!} \\prod_{i=1}^k p_i^{x_i},
\\quad x_1 + \\cdots + x_k = n
```
```julia
Multinomial(n, p)   # Multinomial distribution for n trials with probability vector p
Multinomial(n, k)   # Multinomial distribution for n trials with equal probabilities
                    # over 1:k
```
"""

import Base.*

    function *(B,x::AbstractVector)

    end


function B(power::AbstractFloat=1,
           x::Real,
           xv::AbstractVector)

end

## Generalized Factorial; allows for Real numbers (not Gamma function)
function gfact(x::Real,n::Integer=Int(ceil(x)-1))
    (x == 0) | (n == 0) && return(0.0)
    res = x
    for i in 2:n
        res *= x-(i-1)
    end
    res
end
gfact(6.2,7)
gfact(6) #720
d = 6.2
(d-0)*(d-1)*(d-2)*(d-3)*(d-4)*(d-5)*(d-6)*(d-7)


function O_B(d::Real,n::Integer)
    @assert d > 0
    ombv = [1.0]
    for i in 1:n
        println(ombv[end])
        s = mod(i,2)==0 ? 1 : -1
        push!(ombv, i <= 20 ? s*gfact(d,i)/factorial(i) : 0.0)
    end
    ombv
end

using Polynomials
p = Polynomial(O_B(2.1,22))
q = Polynomial(O_B(2,10))
p ÷ q


O_B(2.1,10)
O_B(2,4)

Polynomial([1,-1])^3

P = Polynomial
∑ = sum
θ = theta
ϕ = phi

# (1-sum([ϕ_1*B^1,ϕ_2*B^2,...,ϕ_p*B^p])*(1-B)^d*X(t) = (1-sum([ϕ_1*B^1,ϕ_2*B^2,...,ϕ_p*B^p])*ϵ(t)

(ϕ1,ϕ2,ϕp) = (1,2,3)

P([1,-ϕ1,-ϕ2,...,-ϕp]))*P(1,-1)^d*X_t = P([1,-θ1,-θ2,...,-θp])*ϵ_t
...

# apply to all suitable windows for X(t) and build a matrix
a0*X_t + a1*X_t-1+ ... + ak*X_t-k = b0*ϵ_t + b1*ϵ_t-1+ ... + bk*ϵ_t-k




