using Documenter, Crystallography

makedocs(;
    modules=[Crystallography],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/singularitti/Crystallography.jl/blob/{commit}{path}#L{line}",
    sitename="Crystallography.jl",
    authors="Qi Zhang",
    assets=String[],
)

deploydocs(;
    repo="github.com/singularitti/Crystallography.jl",
)
