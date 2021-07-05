using Test
using Forecast

@testset "arsim" begin

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

    # arsim(Φ,Φ0,x0,Σ2,n)
    Φ,Φ0,x0,n,Σ2 = 1,0,1,3,0
    x = arsim(Φ,Φ0,x0,n;Σ2)
    @test x ≈ [1, 1, 1]

    Φ,Φ0,x0,n,Σ2 = 1,0,1,3,0
    x = arsim(Φ,Φ0,x0,n;Σ2,fix=[missing,2,missing])
    @test x ≈ [1, 2, 2]
    
    Φ,Φ0,x0,Σ2,n = .1,.2,.3,0.,3
    x = arsim(Φ,Φ0,x0,n;Σ2)
    @test x ≈ [0.23, 0.223, 0.2223]
    
    Φ,Φ0,x0,n,Σ2 = .1,.2,.3,10,.4
    x = arsim(Φ,Φ0,x0,n;Σ2)
    @test x isa Vector
    @test length(x) == 10

    Φ,Φ0,x0,n,Σ2 = [.1,.2],.3,[.4,.5],10,.6
    x = arsim(Φ,Φ0,x0,n;Σ2)
    @test x isa Vector
    @test length(x) == 10

    Φ,Φ0,x0,n,Σ2 = [.1 .2; .3 .4],[.5, .6],[.7,.8],10,[.9 .10; .10 .9]
    x = arsim(Φ,Φ0,x0,n;Σ2)
    @test x isa Matrix
    @test size(x) == (10,2)

    Φ,Φ0,x0,n,Σ2 = reshape([[.1 .2; .3 .4] [.5 .6; .7 .8]],2,2,2),
                  [.9, .10],[.11 .12; .13 .14],10,[.9 .10; .10 .9]
    x = arsim(Φ,Φ0,x0,n;Σ2)
    @test x isa Matrix
    @test size(x) == (10,2)
    
    # arsim(Φ,Φ0,x0,E,n)
    E = MvLogNormal(MvNormal(1,1))
    Φ,Φ0,x0,n = .1,.2,.3,10
    x = arsim(Φ,Φ0,x0,n;E)
    @test x isa Vector
    @test length(x) == 10
    
    Φ,Φ0,x0,n = [.1,.2],.3,[.4,.5],10
    x = arsim(Φ,Φ0,x0,n;E)
    @test x isa Vector
    @test length(x) == 10

    E = MvLogNormal(MvNormal(2,1))
    Φ,Φ0,x0,n = [.1 .2; .3 .4],[.5, .6],[.7,.8],10
    x = arsim(Φ,Φ0,x0,n;E)
    @test x isa Matrix
    @test size(x) == (10,2)

    E = MvNormal(2,1)
    Φ,Φ0,x0,n,Σ2 = reshape([[.1 .2; .3 .4] [.5 .6; .7 .8]],2,2,2),
    [.9, .10],[.11 .12; .13 .14],10,[.9 .10; .10 .9]
    x = arsim(Φ,Φ0,x0,n;E)
    @test x isa Matrix
    @test size(x) == (10,2)

    # arsim(AR,100)
    Φ,Φ0,x0,n,Σ2 = .1,.2,.3,1000,.4
    xar = ar(arsim(Φ,Φ0,x0,n;Σ2),1)
    x = arsim(xar,100)
    @test x isa Vector
    @test length(x) == 100

end
