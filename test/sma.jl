using Test
using Forecast

@testset "sma" begin

    @test_throws AssertionError sma(1:5,10)
    
    x = sma(1:5,3)
    @test length(x) == 3
    @test x == [2.0, 3.0, 4.0]

    x = sma(1:5,3,true)
    @test length(x) == 5
    @test length(collect(skipmissing(x))) == 3
    @test findlast(!ismissing,x) == 4
    @test collect(skipmissing(x)) == [2.0, 3.0, 4.0]

    x = sma(1:5,3,false)
    @test length(x) == 5
    @test length(collect(skipmissing(x))) == 3
    @test findlast(!ismissing,x) == 3
    @test collect(skipmissing(x)) == [2.0, 3.0, 4.0]

    r = rand(100)
    x = sma(r,10)
    y = collect(skipmissing(sma(r,10,false)))
    z = collect(skipmissing(sma(r,10,true)))
    @test x == y == z
    
end
