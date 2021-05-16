using Test
using Forecast
using TimeSeries

@testset "stl" begin
    @test_throws AssertionError stl(rand(100),10; ns=5)
    
    x = stl(rand(100),10)
    @test x isa STL
    @test x.decomposition isa DataFrame
    @test x.call isa String
end
