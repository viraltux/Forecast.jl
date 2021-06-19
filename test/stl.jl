using Test
using Dates
using Forecast

@testset "stl" begin
    @test_throws AssertionError stl(rand(100),10; ns=5)
    
    x = stl(rand(100),10)
    @test x isa STL
    @test x.decomposition isa DataFrame
    @test x.call isa String
    
    x = stl(co2(),365; robust=true, spm=true)
    @test Array(x.decomposition[1:5,:]) == 
    Any[Date("1974-05-17") 3.3336092113804625 330.03520225817715 0.011188530442382216; Date("1974-05-18") 3.2974595453530924 330.0378212346287 -0.22528077998180152; Date("1974-05-19") 3.2594500979279815 330.0404399841902 0.16010991788181173; Date("1974-05-20") 3.219714340570647 330.043058508077 0.3772271513523151; Date("1974-05-21") 3.178405197122059 330.0456768071629 0.1859179957150463]
    
end


