using Test
using Forecast

@testset "arsim" begin

    # arsim(Φ,n)
    Φ,n = .1,10
    x = arsim(Φ,n)
    @test x isa Vector
    @test length(x) == 10

    Φ,n = [.1,.2],10
    x = arsim(Φ,n)
    @test x isa Vector
    @test length(x) == 10

    Φ,n = [.1 .2; .3 .4],10
    x = arsim(Φ,n)
    @test x isa Matrix
    @test size(x) == (10,2)

    # arsim(Φ,Φ0,n)
    Φ,Φ0,n = .1,.2,10
    x = arsim(Φ,Φ0,n)
    @test x isa Vector
    @test length(x) == 10

    Φ,Φ0,n = [.1,.2],.3,10
    x = arsim(Φ,Φ0,n)
    @test x isa Vector
    @test length(x) == 10

    Φ,Φ0,n = [.1 .2; .3 .4],[.5, .6],10
    x = arsim(Φ,Φ0,n)
    @test x isa Matrix
    @test size(x) == (10,2)

    # arsim(Φ,Φ0,x0,n)
    Φ,Φ0,x0,n = .1,.2,.3,10
    x = arsim(Φ,Φ0,x0,n)
    @test x isa Vector
    @test length(x) == 10

    Φ,Φ0,x0,n = [.1,.2],.3,[.4,.5],10
    x = arsim(Φ,Φ0,x0,n)
    @test x isa Vector
    @test length(x) == 10

    Φ,Φ0,x0,n = [.1 .2; .3 .4],[.5, .6],[.7,.8],10
    x = arsim(Φ,Φ0,x0,n)
    @test x isa Matrix
    @test size(x) == (10,2)

    # arsim(Φ,Φ0,x0,Σ,n)
    Φ,Φ0,x0,Σ,n = .1,.2,.3,0,3
    x = arsim(Φ,Φ0,x0,Σ,n)
    @test x ≈ [0.23, 0.223, 0.2223]
    
    Φ,Φ0,x0,Σ,n = .1,.2,.3,.4,10
    x = arsim(Φ,Φ0,x0,Σ,n)
    @test x isa Vector
    @test length(x) == 10

    Φ,Φ0,x0,Σ,n = [.1,.2],.3,[.4,.5],.6,10
    x = arsim(Φ,Φ0,x0,Σ,n)
    @test x isa Vector
    @test length(x) == 10

    Φ,Φ0,x0,Σ,n = [.1 .2; .3 .4],[.5, .6],[.7,.8],[.9 .10; .10 .9],10
    x = arsim(Φ,Φ0,x0,Σ,n)
    @test x isa Matrix
    @test size(x) == (10,2)

    
    Φ,Φ0,x0,Σ,n = reshape([[.1 .2; .3 .4] [.5 .6; .7 .8]],2,2,2),
                  [.9, .10],[.11 .12; .13 .14],[.9 .10; .10 .9],10
    x = arsim(Φ,Φ0,x0,Σ,n)
    @test x isa Matrix
    @test size(x) == (10,2)

    
    # arsim(Φ,Φ0,x0,E,n)
    E1 = MvLogNormal(MvNormal(1,1))
    Φ,Φ0,x0,n = .1,.2,.3,10
    x = arsim(Φ,Φ0,x0,E1,n)
    @test x isa Vector
    @test length(x) == 10
    
    Φ,Φ0,x0,n = [.1,.2],.3,[.4,.5],10
    x = arsim(Φ,Φ0,x0,E1,n)
    @test x isa Vector
    @test length(x) == 10

    E2 = MvLogNormal(MvNormal(2,1))
    Φ,Φ0,x0,n = [.1 .2; .3 .4],[.5, .6],[.7,.8],10
    x = arsim(Φ,Φ0,x0,E2,n)
    @test x isa Matrix
    @test size(x) == (10,2)

    # arsim(AR,100)
    Φ,Φ0,Σ,n = .1,.2,.3,1000
    xar = ar(arsim(Φ,Φ0,Σ,n),1)
    x = arsim(xar,100)
    @test x isa Vector
    @test length(x) == 100

end
