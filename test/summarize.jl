using Test
using Forecast

@testset "summarize" begin

    x = summarize(rand(100))
    @test size(x.quantiles) == (1,8)
    @test size(x.moments) == (1,5)
    @test size(x.format) == (1,2)

    @test x.quantiles isa DataFrame
    @test x.moments isa DataFrame
    @test x.format isa DataFrame
    
    x = summarize(rand(100,2))
    @test size(x.quantiles) == (2,8)
    @test size(x.moments) == (2,5)
    @test size(x.format) == (2,2)
    
end


