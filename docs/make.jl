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
            "Mackey Functors" => "types/MackeyFunctors.md",
            "Mackey Functor homomorphisms" => "types/MackeyFunctorHomomorphisms.md"
        ],
        "Manual" => [
            "Constructors" => "manual/constructors.md",
            "misc" => "manual/misc.md",
            "Examples" => "examples.md"
        ]
    ],
    checkdocs=:exports,
    checkdocs_ignored_modules=[MackeyFunctors.AbstractAlgebraLocal],
)

deploydocs(;
    repo="github.com/tbrazel/MackeyFunctors.jl.git",
    devbranch="main",
)
