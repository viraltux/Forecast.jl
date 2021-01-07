using Test
using Forecast

@testset "hma" begin
    r_data = rand(24)

    # asserting return of correct weights 
    # as listed here: https://www.mathworks.com/help/econ/seasonal-adjustment-using-snxd7m-seasonal-filters.html
    x = hmaSymmetricWeights(13)
    @test length(x) == 13
    @test x         == Float16.([-0.019, -0.028, 0.0, .066, .147, .214, .24, .214, .147, .066, 0.0, -0.028, -0.019])

    # asserting correct application of the 13-term moving average
    x = hma(r_data, 13)
    @test length(x)             == 24
    @test first(x)              == 0.5508635070863032
    @test last(x)               == 0.5843536698022159
    @test x[length(x)//2|>Int]  == 0.32899917067680573
end
