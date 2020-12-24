using Test
using Forecast

@testset "co2" begin
    co2_ta = co2()
    co2_df = co2(true)
    @test co2_ta isa TimeSeries.TimeArray
    @test co2_df isa DataFrames.DataFrame
    @test size(co2_ta) == (4609,)
    @test size(co2_df) == (17166, 17)
end
