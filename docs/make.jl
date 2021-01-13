push!(LOAD_PATH,"../src/")
ENV["GKSwstype"] = "nul"

using Documenter, Forecast

makedocs(
    sitename="Forecast.jl",
    modules = [Forecast],
    pages = [
        "index.md",
        "man/quickstart.md",
        "man/methods.md",
        "man/datasets.md",
        "man/docstrings.md",		
    ]
)

deploydocs(
    repo = "github.com/viraltux/Forecast.jl.git",
    devbranch = "main",
)

