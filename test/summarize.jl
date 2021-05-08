using Test
using Forecast

@testset "summarize" begin
    
    data = (datetime = [DateTime(2018, 11, 21, 12, 0), DateTime(2018, 11, 21, 13, 0)],
            col1 = [10.2, 11.2],
            col2 = [20.2, 21.2],
            col3 = [30.2, 31.2])
    ta = TimeArray(data; timestamp = :datetime, meta = "Test")
    sta = summarize(ta)
    @test sta isa Forecast.SUMMARIZE
    @test sta.quantiles isa DataFrame
    @test sta.moments isa DataFrame
    @test sta.format isa DataFrame

    x = summarize(rand(100))
    @test size(x.quantiles) == (1,8)
    @test size(x.moments) == (1,5)
    @test size(x.format) == (1,2)

    x = summarize(rand(100,2))
    @test size(x.quantiles) == (2,8)
    @test size(x.moments) == (2,5)
    @test size(x.format) == (2,2)
    
end


