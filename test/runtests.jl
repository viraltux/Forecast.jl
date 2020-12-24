using Test
using Forecast
using Random
using Logging

const tests = [
    "sma",
    "loess",
    "stl",
    "datasets"
]

printstyled("\nTest Summary List:\n", color=:underline)

Random.seed!(36)

# include with hidden stderr
function include_hide(file)
    io = open("/dev/null", "w")
    redirect_stderr(io) do
        logger = ConsoleLogger(stderr)
        with_logger(logger) do
            include(file)
        end
    end
end

for t in tests
    @testset "Test $t" begin
        Random.seed!(36)
        include_hide("$t.jl")
        println()
    end
end

println("\nAmbiguous methods:")
display(detect_ambiguities(Forecast, imported=true))
println()
