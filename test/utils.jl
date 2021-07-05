using Test
import Forecast: nextodd, drop, compact, insert_row, insert_column

@testset "utils" begin

    # nextodd
    @test nextodd(5) == 5
    @test nextodd(5.) isa Integer
    @test nextodd(2) == 3

    # drop
    M = [1 2 3
         4 5 6
         7 8 9]
    @test drop(M;r=2,c=2) == [1 3
                              7 9]
    @test drop(M;r=1,c=2) == [4 6
                              7 9]
    M = [1 2 3.
         4 5 6
         7 8 9]
    @test drop(M;r=1,c=3) == [4. 5.
                              7. 8.]

    @test drop(M) == drop(M;r=[],c=[]) == M
    @test drop(M;r=[1,2],c=3) == [7 8]
    @test drop(M;r=[1,2],c=[2,3]) == reshape([7],1,1)

    # compact
    compact(5) ==  5
    compact(5) isa Integer
    compact(5.) ==  5.
    compact(5.) isa AbstractFloat

    x = reshape([1],1,1,1,1,1)
    compact(x) == 1

    x = reshape([1 2],2,1,1,1,1)
    compact(x) == [1,2]
    x = reshape([1 2],1,1,2,1,1)
    compact(x) == [1,2]
    x = reshape([1 2; 3 4],1,1,2,1,2)
    compact(x) == [1 2
                   3 4]

    #insert row / column
    M = zeros(2,2)
    @test insert_row(M, 2, 1) == [0 0
                                  1 1
                                  0 0]
    @test insert_column(M, 2, 1) == [0 1 0
                                     0 1 0]
                            
    @test insert_column(insert_row(M,2,1),2,1) == [0 1 0
                                                   1 1 1
                                                   0 1 0]
end
