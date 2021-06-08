using CrystallographyBase
using Documenter

DocMeta.setdocmeta!(CrystallographyBase, :DocTestSetup, :(using CrystallographyBase); recursive=true)

makedocs(;
    modules=[CrystallographyBase],
    authors="Qi Zhang <singularitti@outlook.com>",
    repo="https://github.com/MineralsCloud/CrystallographyBase.jl/blob/{commit}{path}#{line}",
    sitename="CrystallographyBase.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://MineralsCloud.github.io/CrystallographyBase.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/MineralsCloud/CrystallographyBase.jl",
)
