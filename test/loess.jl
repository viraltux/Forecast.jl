using Test
using Forecast

@testset "loess" begin
    @test_throws AssertionError loess(rand(5),rand(5); d=3)
    
    x = loess(rand(10),rand(10); predict = rand(20))
    @test length(x) == 20

    x = loess(rand(10), rand(10); d=1, predict= rand(20))
    @test length(x) == 20
end
