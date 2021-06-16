using Test
using Forecast

import DataFrames: DataFrame

@testset "p" begin

    x = [1,2,3,4,5]
    x0 = [0,1,2,3,4]
    @test p(d(x)) == x0
    @test p(d(x),[x[1]]) == x
    @test p(d(hcat(x,x,x))) == hcat(x0,x0,x0)
    @test p(d(hcat(x,x,x)), [ x[1] x[1] x[1] ] ) == hcat(x,x,x)

    x = rand(100)
    @test p(d(x),[ x[1] ]) ≈ x
    @test p(d(hcat(x,x,x)), [ x[1] x[1] x[1] ] ) ≈ hcat(x,x,x)

    x = [repeat(1:10,10) repeat(1:10,10)]
    dx = d(x,3)
    x0 = [[1,1,0] [1,1,0]]
    @test p(dx, x0) == x

    la = 11
    x = [repeat(1:la,14) repeat(1:la,14)
         1 1]
    dx = d(x,1,la)
    x0 = reshape(repeat(1:la,inner=2),1,2,la)
    @test p(dx, x0) == x
    
    x = repeat(1:2,30)
    dx = d(x,2,2)
    x0 = zeros(2,1,2) # order 2, 1 variable, lag 1
    x0[1,:,:] = collect(1:2)
    @test p(dx,x0) ≈ x

    tx = hcat(co2().co2, co2().co2, co2().co2)
    for col in eachcol(tx) replace!(col, missing => 0.0) end
    @test p(d(tx), [ 333.38 333.38 333.38 ] ) ≈ tx
    @test p(d(tx,2), [ 333.38 333.38 333.38;
                        -0.27  -0.27  -0.27]    ) ≈ tx

    # calculation of π
    x = 0:0.001:1
    y = sqrt.(1 .- x.^2);
    isapprox(4*p(y)[end]/1000 , π, atol = 0.01)

end
