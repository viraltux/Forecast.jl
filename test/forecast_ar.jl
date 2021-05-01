using Test
using Forecast
import Forecast: fvar

@testset "forecast_ar" begin

    Φ,Σ,n = .991,1,10
    @test fvar(Φ,Σ,n) ≈ [1.0, 1.982081, 2.946564090561, 3.8937646086222375, 4.823992240600336, 5.737551123641018, 6.634739945056495, 7.515852039981027, 8.381175487276607, 9.230993203720098]

    Φ,Σ,n = [1],1,10
    @test fvar(Φ,Σ,n) == [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]


    Φ,Σ,n = [.1,.2,.3,.4],1,10
    @test fvar(Φ,Σ,n) ≈ [1.0, 1.01, 1.0504, 1.144936, 1.32812576, 1.3676603616, 1.394658416256, 1.41080916954496, 1.4192363966542336, 1.4235135889899637]

    Φ,Σ,n = [.1 .2;.3 .4],[1. 0.; 0. 1.],10
    @test fvar(Φ,Σ,n) ≈ [1.0 1.0; 1.05 1.25; 1.0525 1.3125; 1.052625 1.328125; 1.05263125 1.33203125; 1.0526315625 1.3330078125; 1.052631578125 1.333251953125; 1.05263157890625 1.33331298828125; 1.0526315789453125 1.3333282470703125; 1.0526315789472656 1.3333320617675781]
    
    Φ,Σ,n = reshape([[.1 .2; .3 .4] [.5 .6; .7 .8]],2,2,2),[.9 .10; .10 .9],10
    @test fvar(Φ,Σ,n) ≈ [0.8200000000000001 0.8200000000000001; 0.8610000000000001 1.0250000000000001; 1.3862100000000002 2.18325; 1.7086381000000002 3.5433225000000004; 1.9315797410000002 5.369766925000001; 2.0836955470100005 7.773667250250002; 2.1876332707261 10.946685724032502; 2.258641072493421 15.13317168071923; 2.307152717757292 20.657155430220854; 2.3402952114566196 27.945878556329376]
    
    # ar(Φ,n)
    Φ,n = .5,n
    xar = ar(arsim(Φ,n),1)
    fxar = forecast(xar,10)
    @test length(fxar.mean) == 10

end
