using Test
using Forecast

@testset "summarize" begin
    data = (datetime = [DateTime(2018, 11, 21, 12, 0), DateTime(2018, 11, 21, 13, 0)],
            col1 = [10.2, 11.2],
            col2 = [20.2, 21.2],
            col3 = [30.2, 31.2])
    ta = TimeArray(data; timestamp = :datetime, meta = "Example")
    sta = summarize(ta)
    @test sta isa Matrix
    @test size(sta) == (4,7)
    @test sta == Any["🅂" "Min" "1Q" "Median" "Mean" "3Q" "Max";
                     "col1" 10.2 10.45 10.7 10.7 10.95 11.2;
                     "col2" 20.2 20.45 20.7 20.7 20.95 21.2;
                     "col3" 30.2 30.45 30.7 30.7 30.95 31.2]
end


