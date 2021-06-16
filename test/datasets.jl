using Test
using Forecast
import DataFrames: DataFrame

@testset "co2" begin
    co2_paper = co2()
    co2_df = co2(true)
    @test co2_paper isa DataFrame
    @test co2_df isa DataFrame
    @test size(co2_paper) == (4609,2)
    @test size(co2_df) == (17166, 17)
end
