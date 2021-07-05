push!(LOAD_PATH,"../src/")
ENV["GKSwstype"] = "nul"

using Documenter, Forecast

makedocs(
    sitename="Forecast.jl",
    modules = [Forecast],
    pages = [
        "index.md",
        "man/quickstart.md",
        "Forecast Examples" =>[
            "man/examples/quakes.md",
            "man/examples/airpassengers.md",
            "man/examples/london.md",
            "man/examples/whatif.md",
        ],
        "man/methods.md",
        "man/datasets.md",
        "man/docstrings.md",		
    ]
)

deploydocs(
    repo = "github.com/viraltux/Forecast.jl.git",
    devbranch = "main",
)

