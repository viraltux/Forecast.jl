using Test
using Forecast

import DataFrames: DataFrame

@testset "p" begin

    x = [1,2,3,4,5]
    x0 = [0,1,2,3,4]
    @test p(d(x)) == x0
    @test p(d(x),[[x[1]]]) == x
    @test p(d(hcat(x,x,x))) == hcat(x0,x0,x0)
    @test p(d(hcat(x,x,x)), [ [[x[1]]] [[x[1]]] [[x[1]]] ] ) == hcat(x,x,x)

    x = rand(100)
    @test p(d(x),[[x[1]]]) ≈ x
    @test p(d(hcat(x,x,x)), [ [[x[1]]] [[x[1]]] [[x[1]]] ] ) ≈ hcat(x,x,x)

    x = repeat(1:12,3)
    dx = d(x,12,12)
    orderlag = vcat([collect(1:12)],repeat([0],11))
    @test p(dx,orderlag) ≈ x
    @test p(d(x),[[1]]) ≈ x

    tx = hcat(co2().co2, co2().co2, co2().co2)
    for col in eachcol(tx) replace!(col, missing => 0.0) end
    @test p(d(tx), [ [[333.38]] [[333.38]] [[333.38]] ] ) ≈ tx
    @test p(d(tx,2), [ [[333.38]] [[333.38]] [[333.38]] ; -0.27 -0.27 -0.27]) ≈ tx

    # estimation of π
    x = 0:0.001:1
    y = sqrt.(1 .- x.^2)
    isapprox(4*p(y)[end]/1000 , π, atol = 0.01)

end
