using Test
using Forecast
using Distributions: MvLogNormal, MvNormal

@testset "ar" begin
    xar = ar(rand(100),10)
    @test xar isa AR
    @test xar.ic isa Dict
    @test xar.call == "ar(X, order=10, constant=true)"
    @test length(xar.Φ) == 10
    @test length(xar.fitted) == 100 - 10
    @test length(xar.residuals) == 100 - 10

    xar = ar(rand(100), 10; alpha = 0.05)
    @test xar isa AR
    @test xar.ic isa Dict
    @test xar.call == "ar(X, order=10, constant=true)"
    @test length(xar.fitted) == 100 - 10
    @test length(xar.residuals) == 100 - 10
    
    n = 10000
    atol = .2
    
    # ar(Φ,Φ0,x0,n)
    Φ,Φ0,x0,n = .1,.2,.3,n
    xar = ar(arsim(Φ,Φ0,x0,n),1)
    @test isapprox(xar.Φ,Φ; atol = atol)
    @test isapprox(xar.Φ0,Φ0; atol = atol)

    Φ,Φ0,x0,n = [.1,.2],.3,[.4,.5],n
    xar = ar(arsim(Φ,Φ0,x0,n),2)
    @test isapprox(xar.Φ,Φ; atol = atol)
    @test isapprox(xar.Φ0,Φ0; atol = atol)

    Φ,Φ0,x0,n = [.1 .2; .3 .4],[.5, .6],[.7,.8],n
    xar = ar(arsim(Φ,Φ0,x0,n),1)
    @test isapprox(xar.Φ,Φ; atol = 2*atol)
    @test isapprox(xar.Φ0,Φ0; atol = 2*atol)

    # ar(Φ,Φ0,x0,Σ,100)
    Φ,Φ0,x0,Σ2,n = .1,.2,.3,.4,n
    xar = ar(arsim(Φ,Φ0,x0,n;Σ2),1)
    @test isapprox(xar.Φ,Φ; atol = atol)
    @test isapprox(xar.Φ0,Φ0; atol = atol)
    @test isapprox(xar.Σ2,Σ2; atol = atol)

    Φ,Φ0,x0,Σ2,n = [.1,.2],.3,[.4,.5],.6,n
    xar = ar(arsim(Φ,Φ0,x0,n;Σ2),2)
    @test isapprox(xar.Φ,Φ; atol = atol)
    @test isapprox(xar.Φ0,Φ0; atol = atol)
    @test isapprox(xar.Σ2,Σ2; atol = atol)

    Φ,Φ0,x0,Σ2,n = [.1 .2; .3 .4],[.5, .6],[.7,.8],[.9 .10; .10 .9],n
    xar = ar(arsim(Φ,Φ0,x0,n;Σ2),1)
    @test isapprox(xar.Φ,Φ; atol = atol)
    @test isapprox(xar.Φ0,Φ0; atol = atol)
    @test isapprox(xar.Σ2,Σ2; atol = atol)

    Φ,Φ0,x0,Σ2,n = [.1 .2; .3 .4],[.0, .0],[.7,.8],[.9 .10; .10 .9],n
    xar = ar(arsim(Φ,Φ0,x0,n;Σ2),1,false)
    @test isapprox(xar.Φ,Φ; atol = atol)
    @test isapprox(xar.Φ0,Φ0; atol = atol)
    @test isapprox(xar.Σ2,Σ2; atol = atol)
    
    Φ,Φ0,x0,Σ2,n = reshape([[.1 -.1; .05 -.1] [.15 -.4; .3 .6]],2,2,2),
    [.9, .10],[.11 .12; .13 .14],[.9 .10; .10 .9],n
    xar = ar(arsim(Φ,Φ0,x0,n;Σ2),2)
    @test isapprox(xar.Φ,Φ; atol = atol)
    @test isapprox(xar.Φ0,Φ0; atol = atol)
    @test isapprox(xar.Σ2,Σ2; atol = atol)

    Φ,Φ0,x0,Σ2,n = reshape([[.01 -.01 .02; .05 -.01 .03; .01 -.05 .02]
                           [.015 -.04 .06; .03 .03 .06; .01 .04 -.1]],3,3,2),
    [.0, .0, .0],[.11 .12; .2 .13; .14 .6],[1 0.1 0.5; 0.1 1 0.1; 0.5 0.1 1],2*n
    xar = ar(arsim(Φ,Φ0,x0,n;Σ2),2,false)
    @test isapprox(xar.Φ,Φ; atol = atol)
    @test isapprox(xar.Φ0,Φ0; atol = atol)
    @test isapprox(xar.Σ2,Σ2; atol = atol)
    
    # ar(Φ,Φ0,x0,E,100)
    E1 = MvLogNormal(MvNormal(1,1))

    Φ,Φ0,x0,n = .1,.2,.3,n
    xar = ar(log.(arsim(Φ,Φ0,x0,n;E=E1)),1)
    @test xar isa AR
    
    Φ,Φ0,x0,n = [.1,.2],.3,[.4,.5],n
    xar = ar(log.(arsim(Φ,Φ0,x0,n;E=E1)),2)
    @test xar isa AR

    E2 = MvLogNormal(MvNormal(2,1))
    Φ,Φ0,x0,n = [.1 .2; .3 .4],[.5, .6],[.7,.8],n
    xar = ar(log.(arsim(Φ,Φ0,x0,n;E=E2)),1)
    @test xar isa AR

    # ar(Φ,n) fixed
    Φ,Φ0,x0,n = .5,0.0,1.,n
    x = arsim(Φ,Φ0,x0,n)
    xar  = ar(x,10)
    xarfixed = ar(x,10; dΦ=(xar.Φ,xar.Φ), dΦ0=(xar.Φ0, 0.0))
    xarfalse = ar(x,10,false)
    @test xarfixed.Φ == xarfalse.Φ
    @test xarfixed.Φ0 == xarfalse.Φ0
    @test xarfixed.Φse == xarfalse.Φse
    @test xarfixed.Σ2 == xarfalse.Σ2

    Φ,Φ0,x0,n = .5,0.0,1.,n
    x = arsim(Φ,Φ0,x0,n)
    xar  = ar(x,10)
    fΦ = copy(xar.Φ)
    fΦ[3] = 6.0
    fxar = ar(x,10; dΦ=(xar.Φ,fΦ), dΦ0 = (xar.Φ0, 10.0))
    @test fxar.Φ[3] == 6.0
    @test fxar.Φ0 == 10.0
    @test fxar.Φ0se[1] == 0.0
    @test fxar.Φse[3] == 0.0
    
    Φ,Φ0,x0,n = .5,.5,.5,n
    x = arsim(Φ,Φ0,x0,n)
    xar  = ar(x,1)
    fΦ = copy(xar.Φ)
    fΦ = 1.0
    @test_throws AssertionError ar(x,1; dΦ=(xar.Φ,fΦ), dΦ0 = (xar.Φ0, [1.0]))
    
    Φ,Φ0,x0,Σ2,n = reshape([[.01 -.01 .02; .05 -.01 .03; .01 -.05 .02]
                           [.015 -.04 .06; .03 .03 .06; .01 .04 -.1]],3,3,2),
    [.9, .10, .5],[.11 .12; .2 .13; .14 .6],[1 0.1 0.1; 0.1 1 0.1; 0.1 0.1 1],n
    x = arsim(Φ,Φ0,x0,n;Σ2)
    xar  = ar(x,2)
    fΦ = copy(xar.Φ)
    fΦ[1,3,1] = 6.0
    fΦ[3,1,2] = 9.0
    fxar = ar(x,2; dΦ=(xar.Φ,fΦ))
    @test fxar.Φ[1,3,1] == 6.0
    @test fxar.Φ[3,1,2] == 9.0

end
