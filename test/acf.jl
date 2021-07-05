using Test
import Forecast: acf

@testset "ccf" begin
    x = rand(100);
    @test_throws AssertionError acf(x; lag = 100)

    res = acf(x; lag = 20, type="cor");
    @test res isa CCF
    @test res.call == "ccf(x; type=\"cor\", lag=20, alpha=(0.95, 0.99))"
    @test res.type == "cor"
    @test res.auto
    @test res.lag == 20
    @test res.N == 100
    @test res.alpha == (0.95, 0.99)
    @test res.ci[1] ≈ 0.19900419155027554
    @test res.ci[2] ≈ 0.2615358405397182
end
