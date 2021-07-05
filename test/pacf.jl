using Test
using Forecast

@testset "pacf" begin
    x = rand(100);
    @test_throws AssertionError Forecast.pacf(x; lag = 100)

    res = pacf(x; lag = 20);
    @test res isa CCF
    @test res.call == "pacf(x; type=\"stepwise-real\", lag=20, alpha=(0.95, 0.99))"
    @test res.type == "pacf_stepwise-real"
    @test res.auto
    @test res.lag == 20
    @test res.N == 100
    @test res.alpha == (0.95, 0.99)
    @test res.ci[1] ≈ 0.22335862551831387
    @test res.ci[2] ≈ 0.293542992294061
end
