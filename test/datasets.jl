using Test
using Forecast
import DataFrames: DataFrame

@testset "co2" begin
    co2_ta = co2()
    co2_df = co2(true)
    @test co2_ta isa TimeArray
    @test co2_df isa DataFrame
    @test size(co2_ta) == (4609,)
    @test size(co2_df) == (17166, 17)
end
