using Test
using Forecast

@testset "arsim" begin
    x = arsim(1,10)
    @test x isa AbstractArray
end
