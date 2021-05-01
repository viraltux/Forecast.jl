using Test
using Forecast
using Random
using Logging

const tests = [
    "acf",
    "ar",
    "arsim",
    "ccf",
    "d",
    "datasets",
    "forecast_ar",
    "hma",
    "loess",
    "p",
    "pacf",
    "sma",
    "stl",
    "summarize",
]

printstyled("\nTest Summary List:\n", color=:underline)

Random.seed!(36)
Base.disable_logging(Base.CoreLogging.Error) # disable warnings

for t in tests
    @testset "Test $t" begin
        Random.seed!(36)
        include("$t.jl")
        println()
    end
end

println()
