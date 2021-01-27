using Test
using Forecast

import DataFrames: DataFrame
import TimeSeries: TimeArray

@testset "p" begin

    x = [1,2,3,4,5]
    x0 = [0,1,2,3,4]
    @test p(d(x)) == x0
    @test p(d(x); orderlag=[[x[1]]]) == x
    @test p(d(hcat(x,x,x))) == hcat(x0,x0,x0)
    @test p(d(hcat(x,x,x)); orderlag=[ [[x[1]]] [[x[1]]] [[x[1]]] ]  ) == hcat(x,x,x)

    x = rand(100)
    @test p(d(x),orderlag=[[x[1]]]) ≈ x
    @test p(d(hcat(x,x,x)), orderlag=[ [[x[1]]] [[x[1]]] [[x[1]]] ] ) ≈ hcat(x,x,x)

    x = rand(100,3)
    @test p(d(x), orderlag=[ [[x[1,1]]] [[x[1,2]]] [[x[1,3]]] ] ) ≈ x

    x = repeat(1:12,3)
    dx = d(x; order=12, lag=12)
    orderlag = vcat([collect(1:12)],repeat([0],11))
    @test p(dx; orderlag) == x
    @test p(d(x); orderlag = [[1]]) == x

    tx = hcat(co2(), co2(), co2())
    vtx = values(tx)
    tx = TimeArray(timestamp(tx),replace(vtx, missing => 0.0))
    @test values(p(d(tx), orderlag=[ [[333.38]] [[333.38]] [[333.38]] ] )) ≈ values(tx)
    @test values(p(d(tx; order=2), orderlag=[ [[333.38]] [[333.38]] [[333.38]] ; -0.27 -0.27 -0.27] )) ≈ values(tx)

    # estimation of π
    x = 0:0.001:1
    y = sqrt.(1 .- x.^2)
    isapprox(4*p(y)[end]/1000 , π, atol = 0.01)

end
