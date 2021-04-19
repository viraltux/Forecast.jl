using Test
using Forecast

@testset "predict_ar" begin
    x = predict(ar(rand(100),10), 10)
    @test x isa AbstractArray
end
