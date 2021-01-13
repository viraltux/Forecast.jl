using Test
using Forecast

@testset "sma" begin
    @test_throws AssertionError sma(1:5,10)
    x = sma(1:5,3)
    @test length(x) == 5
    @test length(collect(skipmissing(x))) == 3
    @test findlast(!ismissing,x) == 4
    x = sma(1:5,3;center=false)
    @test length(x) == 5
    @test length(collect(skipmissing(x))) == 3
    @test findlast(!ismissing,x) == 3
end
