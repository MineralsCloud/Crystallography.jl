using Crystallography
using Documenter

DocMeta.setdocmeta!(Crystallography, :DocTestSetup, :(using Crystallography); recursive=true)

makedocs(;
    modules=[Crystallography],
    authors="singularitti <singularitti@outlook.com> and contributors",
    repo="https://github.com/MineralsCloud/Crystallography.jl/blob/{commit}{path}#{line}",
    sitename="Crystallography.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://MineralsCloud.github.io/Crystallography.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/MineralsCloud/Crystallography.jl",
    devbranch="main",
)
