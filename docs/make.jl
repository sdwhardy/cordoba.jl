using cordoba
using Documenter

DocMeta.setdocmeta!(cordoba, :DocTestSetup, :(using cordoba); recursive=true)

makedocs(;
    modules=[cordoba],
    authors="Su-",
    repo="https://github.com/sdwhardy/cordoba.jl/blob/{commit}{path}#{line}",
    sitename="cordoba.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://sdwhardy.github.io/cordoba.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/sdwhardy/cordoba.jl",
)
