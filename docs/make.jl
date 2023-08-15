using Crystallography
using CrystallographyBase
using CrystallographyCore
using Documenter

# See https://stackoverflow.com/questions/70137119/how-to-include-the-docstring-for-a-function-from-another-package-in-my-julia-doc
DocMeta.setdocmeta!(CrystallographyCore, :DocTestSetup, :(using CrystallographyCore); recursive=true)
DocMeta.setdocmeta!(CrystallographyBase, :DocTestSetup, :(using CrystallographyBase); recursive=true)
DocMeta.setdocmeta!(Crystallography, :DocTestSetup, :(using Crystallography); recursive=true)

makedocs(;
    modules=[CrystallographyCore, CrystallographyBase, Crystallography],
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
        "Manual" => [
            "Installation Guide" => "man/installation.md",
            # "Definitions and conventions" => "man/definitions.md",
            # "Examples" => "man/examples.md",
            "Troubleshooting" => "man/troubleshooting.md",
        ],
        "Reference" => Any[
            "Public API" => "lib/public.md",
            # "Internals" => map(
            #     s -> "lib/internals/$(s)",
            #     sort(readdir(joinpath(@__DIR__, "src/lib/internals")))
            # ),
        ],
        "Developer Docs" => [
            "Contributing" => "developers/contributing.md",
            "Style Guide" => "developers/style-guide.md",
            "Design Principles" => "developers/design-principles.md",
        ],
    ],
)

deploydocs(;
    repo="github.com/MineralsCloud/Crystallography.jl",
    devbranch="main",
)
