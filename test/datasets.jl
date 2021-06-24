using Test
using Forecast
import DataFrames: DataFrame

@testset "datasets: co2" begin

    co2_paper = co2()
    co2_df = co2(true)
    @test co2_paper isa DataFrame
    @test co2_df isa DataFrame
    @test size(co2_paper) == (4609,2)
    @test size(co2_df) == (17166, 17)

end

@testset "datasets: seaborne" begin

    sb = seaborne()
    sb_full = seaborne(true)
    @test sb isa DataFrame
    @test sb_full isa DataFrame
    @test size(sb) == (2199,4)
    @test size(sb_full) == (13209, 10)

end

@testset "datasets: air" begin

    a = air()
    @test a isa DataFrame
    @test size(a) == (144,2)

end

@testset "datasets: quakes" begin

    q = quakes()
    @test q isa DataFrame
    @test size(q) == (71,2)

end
