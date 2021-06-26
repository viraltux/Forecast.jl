using Test
using Forecast
using Plots, RecipesBase

@testset "plots" begin

    # plot DataFrame
    @test typeof(plot(co2())) == Plots.Plot{Plots.GRBackend}

    # plot CCF
    @test typeof(plot(ccf(rand(100), rand(100)))) == Plots.Plot{Plots.GRBackend}

    # plot STL
    stl_co2 = stl(co2(),365; robust=true, spm=true)
    @test typeof(plot(stl_co2)) == Plots.Plot{Plots.GRBackend}
    
    # plot FORECAST
    @test typeof(plot(forecast(ar(rand(100),1),10))) == Plots.Plot{Plots.GRBackend}

end
