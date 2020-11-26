using Forecast
using Documenter

makedocs(;
    modules=[Forecast],
    authors="Fran Urbano <viraltux@gmail.com> and contributors",
    repo="https://github.com/viraltux/Forecast.jl/blob/{commit}{path}#L{line}",
    sitename="Forecast.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://viraltux.github.io/Forecast.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/viraltux/Forecast.jl",
)
