using Documenter
using MackeyFunctors

makedocs(;
    modules=[MackeyFunctors],
    sitename="MackeyFunctors.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://tbrazel.github.io/MackeyFunctors.jl/stable/",
    ),
    pages=[
        "Home" => "index.md",
        "Examples" => "examples.md",
        "API" => "api.md",
    ],
    checkdocs=:exports,
)

deploydocs(;
    repo="github.com/tbrazel/MackeyFunctors.jl.git",
    devbranch="main",
)
