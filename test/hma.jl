# MIT License

# Copyright (c) 2021 Fran Urbano
# Copyright (c) 2021 Val Lyashov

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

using Test
using Forecast
import Forecast: hmaSymmetricWeights

@testset "hma" begin

    # weights
    w = hmaSymmetricWeights(8)
    @test ismissing(w[end])

    w = hmaSymmetricWeights(7)
    @test sum(w) ≈ 1.0

    w = hmaSymmetricWeights(13)
    @test sum(w) ≈ 1.0
    # (https://www.mathworks.com/help/econ/seasonal-adjustment-using-snxd7m-seasonal-filters.html)
    @test length(w) == 13
    @test round.(w,digits=3) ≈ [-0.019, -0.028, 0.0, 0.065, 0.147, 0.214, 0.24,
                                 0.214, 0.147, 0.065, 0.0, -0.028, -0.019]  

    # 13-term moving average
    x = hma(sin.(1:100), 13)
    @test length(x) == 100
    @test first(x)  ≈ 0.5636212576875559
    @test last(x)   ≈ -0.6658500507547161
    @test x[50]     ≈ -0.043655947941143455

end
