using Test
using Forecast

import DataFrames: DataFrame
import TimeSeries: TimeArray

@testset "d" begin

    x = [1,2,3,4,5]
    dx = d(x; center=true)
    mdx = ismissing.(dx)
    @test mdx == [0,0,0,0,1]
    @test dx[.!Bool.(mdx)] == [1,1,1,1]

    dx = d(x; order=2, center=true)
    mdx = ismissing.(dx)
    @test mdx == [1,0,0,0,1]
    @test dx[.!Bool.(mdx)] == [0,0,0]

    dx = d(x; order=1, lag=2)
    mdx = ismissing.(dx)
    @test mdx == [0,0,0]
    @test dx[.!Bool.(mdx)] == [2,2,2]

    x = reshape(collect(1:20),10,2)
    dx = d(x; order=2, lag=2)
    mdx = ismissing.(dx)
    @test mdx == [0 0; 0 0; 0 0; 0 0; 0 0; 0 0; 0 0]
    @test dx == [0 0; 0 0; 0 0; 0 0; 0 0; 0 0; 0 0]

    x = TimeArray(Date(1970, 1, 1):Day(1):Date(1970, 1, 7), [1 11; 2 22; 3 33; 4 44; missing missing; missing missing; 5 55])
    dx = values(d(x))
    @test dx[1:3,:] == [1 11; 1 11; 1 11]
    @test ismissing.(dx[4:6,:]) == [1 1; 1 1; 1 1]
    
end
