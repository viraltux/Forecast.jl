using Test
using Forecast

@testset "ccf" begin
    x1 = rand(100);
    x2 = circshift(x1,6);
    @test_throws AssertionError ccf(x1,x2; lag = 100)

    res = ccf(x1, x2; lag = 20, type="cor");
    @test res isa CCF
    @test res.call == "ccf(x1, x2; type=\"cor\", lag=20, alpha=(0.95, 0.99))"
    @test res.type == "cor"
    @test !res.auto
    @test res.lag == 20
    @test res.N == 100
    @test res.alpha == (0.95, 0.99)
    @test res.ci[1] ≈ 0.19900419155027554
    @test res.ci[2] ≈ 0.2615358405397182

    @test res.ccf[20+1+6] ≈ 0.9685785099293228
end
