using Test
using Forecast

@testset "boxcox" begin

    x = rand(100)
    bc = boxcox(x)
    @test iboxcox(bc[:x],bc[:λ]) ≈ x

    x = rand(100) .- 0.5
    bc = boxcox(x)
    @test iboxcox(bc[:x],bc[:λ]) ≈ x

    x = rand(100)
    bc = boxcox(x)
    @test iboxcox(bc[:x],bc[:λ]) ≈ x

    x =  Int64.(round.(100*rand(100)))
    bc = boxcox(x)
    @test iboxcox(bc[:x],bc[:λ]) ≈ x

    x =  Int64.(round.(100*rand(100) .- 50))
    bc = boxcox(x)
    @test iboxcox(bc[:x],bc[:λ]) ≈ x
    
end
