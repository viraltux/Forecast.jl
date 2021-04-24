using Test
using Forecast

@testset "ar" begin
    x = ar(rand(100),10)
    @test x isa AR
    @test x.IC isa Dict
    @test x.call == "ar(X, order=10, constant=true)"
    @test x.constant == x.Φ0
    @test x.constant isa Number
    @test x.stdev ==  x.Σ
    @test x.stdev isa Number
    @test size(x.PC) == (11,11)
    @test x.coefficients == x.Φ
    @test length(x.coefficients) == 10
    @test length(x.fitted) == 100 - 10
    @test length(x.residuals) == 100 - 10

    n = 1000
    atol = .15
    
    # ar(Φ,n)
    Φ,n = .5,n
    x = ar(arsim(Φ,n),1)
    @test isapprox(x.Φ,Φ; atol = atol)
    @test isapprox(x.Φ0,0; atol = atol)
    
    Φ,n = [.1,.2],n
    x = ar(arsim(Φ,n),2)
    @test isapprox(x.Φ,Φ; atol = atol)
    @test isapprox(x.Φ0,0; atol = atol)

    Φ,n = [.1 .2; .3 .4],n
    x = ar(arsim(Φ,n),1)
    @test isapprox(x.Φ,Φ; atol = atol)
    @test isapprox(x.Φ0,zeros(2); atol = atol)
    
    # ar(Φ,Φ0,n)
    Φ,Φ0,n = .1,.2,n
    x = ar(arsim(Φ,Φ0,n),1)
    @test isapprox(x.Φ,Φ; atol = atol)
    @test isapprox(x.Φ0,Φ0; atol = atol)

    Φ,Φ0,n = [.1,.2],.3,n
    x = ar(arsim(Φ,Φ0,n),2)
    @test isapprox(x.Φ,Φ; atol = atol)
    @test isapprox(x.Φ0,Φ0; atol = atol)

    Φ,Φ0,n = [.1 .2; .3 .4],[.5, .6],n
    x = ar(arsim(Φ,Φ0,n),1)
    @test isapprox(x.Φ,Φ; atol = atol)
    @test isapprox(x.Φ0,Φ0; atol = atol)

    # ar(Φ,Φ0,x0,n)
    Φ,Φ0,x0,n = .1,.2,.3,n
    x = ar(arsim(Φ,Φ0,x0,n),1)
    @test isapprox(x.Φ,Φ; atol = atol)
    @test isapprox(x.Φ0,Φ0; atol = atol)

    Φ,Φ0,x0,n = [.1,.2],.3,[.4,.5],n
    x = ar(arsim(Φ,Φ0,x0,n),2)
    @test isapprox(x.Φ,Φ; atol = atol)
    @test isapprox(x.Φ0,Φ0; atol = atol)

    Φ,Φ0,x0,n = [.1 .2; .3 .4],[.5, .6],[.7,.8],n
    x = ar(arsim(Φ,Φ0,x0,n),1)
    @test isapprox(x.Φ,Φ; atol = atol)
    @test isapprox(x.Φ0,Φ0; atol = atol)

    # ar(Φ,Φ0,x0,Σ,100)
    Φ,Φ0,x0,Σ,n = .1,.2,.3,.4,n
    x = ar(arsim(Φ,Φ0,x0,Σ,n),1)
    @test isapprox(x.Φ,Φ; atol = atol)
    @test isapprox(x.Φ0,Φ0; atol = atol)
    @test isapprox(x.Σ,Σ; atol = atol)

    Φ,Φ0,x0,Σ,n = [.1,.2],.3,[.4,.5],.6,n
    x = ar(arsim(Φ,Φ0,x0,Σ,n),2)
    @test isapprox(x.Φ,Φ; atol = atol)
    @test isapprox(x.Φ0,Φ0; atol = atol)
    @test isapprox(x.Σ,Σ; atol = atol)

    Φ,Φ0,x0,Σ,n = [.1 .2; .3 .4],[.5, .6],[.7,.8],[.9 .10; .10 .9],n
    x = ar(arsim(Φ,Φ0,x0,Σ,n),1)
    @test isapprox(x.Φ,Φ; atol = atol)
    @test isapprox(x.Φ0,Φ0; atol = atol)
    @test isapprox(x.Σ,Σ; atol = atol)

    Φ,Φ0,x0,Σ,n = reshape([[.1 -.1; .05 -.1] [.15 -.4; .3 .6]],2,2,2),
    [.9, .10],[.11 .12; .13 .14],[.9 .10; .10 .9],n
    x = ar(arsim(Φ,Φ0,x0,Σ,n),2)
    @test isapprox(x.Φ,Φ; atol = atol)
    @test isapprox(x.Φ0,Φ0; atol = atol)
    @test isapprox(x.Σ,Σ; atol = atol)

    # ar(Φ,Φ0,x0,E,100)
    E1 = MvLogNormal(MvNormal(1,1))

    Φ,Φ0,x0,n = .1,.2,.3,n
    x = ar(log.(arsim(Φ,Φ0,x0,E1,n)),1)
    @test isapprox(x.Φ,Φ; atol = .25)
    @test isapprox(x.Φ0,Φ0; atol = .25)
    @test isapprox(x.Σ,1; atol = .5)
    
    Φ,Φ0,x0,n = [.1,.2],.3,[.4,.5],n
    x = ar(log.(arsim(Φ,Φ0,x0,E1,n)),2)
    @test isapprox(x.Φ,Φ; atol = .25)
    @test isapprox(x.Φ0,Φ0; atol = .25)
    @test isapprox(x.Σ,1; atol = .5)

    E2 = MvLogNormal(MvNormal(2,1))
    Φ,Φ0,x0,n = [.1 .2; .3 .4],[.5, .6],[.7,.8],n
    x = ar(log.(arsim(Φ,Φ0,x0,E2,n)),1)
    @test isapprox(x.Φ,Φ; atol = .25)
    @test isapprox(x.Φ0,Φ0; atol = .25)
    @test isapprox(x.Σ, [1 0; 0 1], atol = 1)
    
end
