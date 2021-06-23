using Test
using Forecast
import Forecast: hmaSymmetricWeights

@testset "hma" begin

    # weights (https://www.mathworks.com/help/econ/seasonal-adjustment-using-snxd7m-seasonal-filters.html)
    x = hmaSymmetricWeights(13)
    @test length(x) == 13
    @test x         == [-0.019, -0.028, 0.0, 0.065, 0.147, 0.214, 0.24, 0.214, 0.147, 0.065, 0.0, -0.028, -0.019]  

    # 13-term moving average
    x = hma(sin.(1:100), 13)
    @test length(x) == 100
    @test first(x)  ≈ 0.5634778709798388
    @test last(x)   ≈ -0.664567410322129
    @test x[50]     ≈ -0.044034801763622136

end
