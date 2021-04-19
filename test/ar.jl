using Test
using Forecast
using TimeSeries

@testset "ar" begin
    x = ar(rand(100),10)
    @test x isa AR
end
