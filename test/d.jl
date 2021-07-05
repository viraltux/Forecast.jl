using Test
using Forecast

import DataFrames: DataFrame

@testset "d" begin

    x = [1,2,3,4,5]
    dx = d(x,center=true)
    mdx = ismissing.(dx)
    @test mdx == [0,0,0,0,1]
    @test dx[.!Bool.(mdx)] == [1,1,1,1]

    dx = d(x,2,center=true)
    mdx = ismissing.(dx)
    @test mdx == [1,0,0,0,1]
    @test dx[.!Bool.(mdx)] == [0,0,0]

    dx = d(x,1,2)
    mdx = ismissing.(dx)
    @test mdx == [0,0,0]
    @test dx[.!Bool.(mdx)] == [2,2,2]

    x = reshape(collect(1:20),10,2)
    dx = d(x,2,2)
    @test dx == repeat([0],6,2)

end
