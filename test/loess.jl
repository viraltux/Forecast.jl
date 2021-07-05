using Test
using Forecast

@testset "loess" begin
    @test_throws AssertionError loess(rand(5),rand(5); d=3)
    
    x = loess(rand(10),rand(10); predict = rand(20))
    @test length(x) == 20

    x = loess(rand(10), rand(10); d=1, predict= rand(20))
    @test length(x) == 20

    x = loess(sin.(collect(1:10.)), sin.(collect(1:10.)); d=2)
    @test x == [0.8414709848078965, 0.9092974268256817, 0.14112000805986716, -0.7568024953079282, -0.9589242746631381, -0.2794154981989261, 0.6569865987187891, 0.9893582466233829, 0.4121184852417571, -0.5440211108893702]

    x = loess(sin.(collect(1:10.)), sin.(collect(1:10.)); d=1)
    @test x == [0.8414709848078965, 0.9092974268256815, 0.14112000805986719, -0.7568024953079282, -0.9589242746631382, -0.2794154981989259, 0.6569865987187891, 0.9893582466233817, 0.41211848524175654, -0.5440211108893699]

    x = loess(Int64.(round.(100*rand(100))),
              Int64.(round.(100*rand(100))); predict = Int64.(round.(100*rand(20))))
    @test length(x) == 20

    x = loess(1:2:200,
              vcat(rand(80),repeat([missing],20));
              predict = 1.:2:40)
    @test length(x) == 20
    
end
