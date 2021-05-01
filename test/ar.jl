using Test
using Forecast
import Distributions: MvLogNormal, MvNormal

@testset "ar" begin
    xar = ar(rand(100),10)
    @test xar isa AR
    @test xar.IC isa Dict
    @test xar.call == "ar(X, order=10, constant=true)"
    @test xar.constant == xar.Φ0
    @test xar.constant isa Number
    @test xar.stdev ==  xar.Σ
    @test xar.stdev isa Number
    @test size(xar.PC) == (11,11)
    @test xar.coefficients == xar.Φ
    @test length(xar.coefficients) == 10
    @test length(xar.fitted) == 100 - 10
    @test length(xar.residuals) == 100 - 10

    n = 1000
    atol = .15
    
    # ar(Φ,n)
    Φ,n = .5,n
    xar = ar(arsim(Φ,n),1)
    @test isapprox(xar.Φ,Φ; atol = atol)
    @test isapprox(xar.Φ0,0; atol = atol)

    Φ,n = .5,n
    xar = ar(arsim(Φ,n),1,false)
    @test isapprox(xar.Φ,Φ; atol = atol)
    @test isapprox(xar.Φ0,0; atol = atol)
    
    Φ,n = [.1,.2],n
    xar = ar(arsim(Φ,n),2)
    @test isapprox(xar.Φ,Φ; atol = atol)
    @test isapprox(xar.Φ0,0; atol = atol)

    Φ,n = [.1 .2; .3 .4],n
    xar = ar(arsim(Φ,n),1)
    @test isapprox(xar.Φ,Φ; atol = atol)
    @test isapprox(xar.Φ0,zeros(2); atol = atol)
    
    # ar(Φ,Φ0,n)
    Φ,Φ0,n = .1,.2,n
    xar = ar(arsim(Φ,Φ0,n),1)
    @test isapprox(xar.Φ,Φ; atol = atol)
    @test isapprox(xar.Φ0,Φ0; atol = atol)

    Φ,Φ0,n = [.1,.2],.3,n
    xar = ar(arsim(Φ,Φ0,n),2)
    @test isapprox(xar.Φ,Φ; atol = atol)
    @test isapprox(xar.Φ0,Φ0; atol = atol)

    Φ,Φ0,n = [.1 .2; .3 .4],[.5, .6],n
    xar = ar(arsim(Φ,Φ0,n),1)
    @test isapprox(xar.Φ,Φ; atol = atol)
    @test isapprox(xar.Φ0,Φ0; atol = atol)

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
    Φ,Φ0,x0,Σ,n = .1,.2,.3,.4,n
    xar = ar(arsim(Φ,Φ0,x0,Σ,n),1)
    @test isapprox(xar.Φ,Φ; atol = atol)
    @test isapprox(xar.Φ0,Φ0; atol = atol)
    @test isapprox(xar.Σ,Σ; atol = atol)

    Φ,Φ0,x0,Σ,n = [.1,.2],.3,[.4,.5],.6,n
    xar = ar(arsim(Φ,Φ0,x0,Σ,n),2)
    @test isapprox(xar.Φ,Φ; atol = atol)
    @test isapprox(xar.Φ0,Φ0; atol = atol)
    @test isapprox(xar.Σ,Σ; atol = atol)

    Φ,Φ0,x0,Σ,n = [.1 .2; .3 .4],[.5, .6],[.7,.8],[.9 .10; .10 .9],n
    xar = ar(arsim(Φ,Φ0,x0,Σ,n),1)
    @test isapprox(xar.Φ,Φ; atol = atol)
    @test isapprox(xar.Φ0,Φ0; atol = atol)
    @test isapprox(xar.Σ,Σ; atol = atol)

    Φ,Φ0,x0,Σ,n = reshape([[.1 -.1; .05 -.1] [.15 -.4; .3 .6]],2,2,2),
    [.9, .10],[.11 .12; .13 .14],[.9 .10; .10 .9],n
    xar = ar(arsim(Φ,Φ0,x0,Σ,n),2)
    @test isapprox(xar.Φ,Φ; atol = atol)
    @test isapprox(xar.Φ0,Φ0; atol = atol)
    @test isapprox(xar.Σ,Σ; atol = atol)

    Φ,Φ0,x0,Σ,n = reshape([[.01 -.01 .02; .05 -.01 .03; .01 -.05 .02]
                           [.015 -.04 .06; .03 .03 .06; .01 .04 -.1]],3,3,2),
    [.9, .10, .5],[.11 .12; .2 .13; .14 .6],[1 0.1 0.1; 0.1 1 0.1; 0.1 0.1 1],n
    xar = ar(arsim(Φ,Φ0,x0,Σ,n),2)
    @test isapprox(xar.Φ,Φ; atol = 2*atol)
    @test isapprox(xar.Φ0,Φ0; atol = 2*atol)
    @test isapprox(xar.Σ,Σ; atol = 2*atol)
    
    # ar(Φ,Φ0,x0,E,100)
    E1 = MvLogNormal(MvNormal(1,1))

    Φ,Φ0,x0,n = .1,.2,.3,n
    xar = ar(log.(arsim(Φ,Φ0,x0,E1,n)),1)
    
    Φ,Φ0,x0,n = [.1,.2],.3,[.4,.5],n
    xar = ar(log.(arsim(Φ,Φ0,x0,E1,n)),2)

    E2 = MvLogNormal(MvNormal(2,1))
    Φ,Φ0,x0,n = [.1 .2; .3 .4],[.5, .6],[.7,.8],n
    xar = ar(log.(arsim(Φ,Φ0,x0,E2,n)),1)
    
end
