using Crystallography
using Documenter

makedocs(;
    modules=[Crystallography],
    authors="Qi Zhang <singularitti@outlook.com>",
    repo="https://github.com/singularitti/Crystallography.jl/blob/{commit}{path}#L{line}",
    sitename="Crystallography.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://singularitti.github.io/Crystallography.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/singularitti/Crystallography.jl",
)
