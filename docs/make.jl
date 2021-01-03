push!(LOAD_PATH,"../src/")
using Documenter, Forecast

makedocs(sitename="Forecast.jl")

deploydocs(
    repo = "github.com/viraltux/Forecast.jl.git",
    devbranch = "main"
)

