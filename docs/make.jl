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
        "Installation" => "installation.md",
        "Contributing" => "contributing.md",
        "Types and design" => [
            "Mackey Contexts" => "types/MackeyContexts.md",
            "Mackey Functors" => "types/MackeyFunctors.md"
        ],
        "Manual" => [
            "Abstract algebra" => "manual/alg.md",
            "Constructors" => "manual/constructors.md",
            "misc" => "manual/misc.md",
            "Examples" => "examples.md"
        ]
    ],
    checkdocs=:exports,
)

deploydocs(;
    repo="github.com/tbrazel/MackeyFunctors.jl.git",
    devbranch="main",
)
