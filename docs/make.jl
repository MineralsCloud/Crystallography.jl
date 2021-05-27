using Crystallography
using Documenter

DocMeta.setdocmeta!(
    Crystallography,
    :DocTestSetup,
    :(using Crystallography);
    recursive = true,
)

makedocs(;
    modules = [Crystallography],
    authors = "Qi Zhang <singularitti@outlook.com>",
    repo = "https://github.com/MineralsCloud/Crystallography.jl/blob/{commit}{path}#{line}",
    sitename = "Crystallography.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://mineralscloud.github.io/Crystallography.jl",
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        "Manual" => Any["Installation"=>"install.md", "Development"=>"develop.md"],
        "API by modules" => Any[
            "`Crystallography` module"=>"Crystallography.md",
            "`Crystallography.Symmetry` module"=>"Symmetry.md",
        ],
    ],
)

deploydocs(; repo = "github.com/MineralsCloud/Crystallography.jl")
