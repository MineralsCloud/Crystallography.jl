using CrystallographyBase
using Documenter

DocMeta.setdocmeta!(CrystallographyBase, :DocTestSetup, :(using CrystallographyBase); recursive=true)

makedocs(;
    modules=[CrystallographyBase],
    authors="singularitti <singularitti@outlook.com> and contributors",
    repo="https://github.com/MineralsCloud/CrystallographyBase.jl/blob/{commit}{path}#{line}",
    sitename="CrystallographyBase.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://MineralsCloud.github.io/CrystallographyBase.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Manual" => [
            "Installation guide" => "installation.md",
        ],
        "API Reference" => "public.md",
        "Developer Docs" => [
            "Contributing" => "developers/contributing.md",
            "Style Guide" => "developers/style.md",
        ],
        "Troubleshooting" => "troubleshooting.md",
    ],
)

deploydocs(;
    repo="github.com/MineralsCloud/CrystallographyBase.jl",
    devbranch="main",
)
