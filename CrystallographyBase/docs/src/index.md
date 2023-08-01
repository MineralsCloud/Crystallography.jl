# CrystallographyBase

Documentation for [CrystallographyBase](https://github.com/MineralsCloud/CrystallographyBase.jl).

See the [Index](@ref main-index) for the complete list of documented functions
and types.

The code, which is [hosted on GitHub](https://github.com/MineralsCloud/CrystallographyBase.jl), is tested
using various continuous integration services for its validity.

This repository is created and maintained by
[@singularitti](https://github.com/singularitti), and contributions are highly welcome.

## Package features

- Define crystal lattices and cells
- Compute reciprocal lattices
- Generate supercells and k-point grids
- Enable math utilities for crystal properties

## Installation

The package can be installed with the Julia package manager.
From [the Julia REPL](https://docs.julialang.org/en/v1/stdlib/REPL/), type `]` to enter
the [Pkg mode](https://docs.julialang.org/en/v1/stdlib/REPL/#Pkg-mode) and run:

```julia-repl
pkg> add CrystallographyBase
```

Or, equivalently, via [`Pkg.jl`](https://pkgdocs.julialang.org/v1/):

```@repl
import Pkg; Pkg.add("CrystallographyBase")
```

## Documentation

- [**STABLE**](https://MineralsCloud.github.io/CrystallographyBase.jl/stable) — **documentation of the most recently tagged version.**
- [**DEV**](https://MineralsCloud.github.io/CrystallographyBase.jl/dev) — _documentation of the in-development version._

## Project status

The package is developed for and tested against Julia `v1.6` and above on Linux, macOS, and
Windows.

## Questions and contributions

You can post usage questions on
[our discussion page](https://github.com/MineralsCloud/CrystallographyBase.jl/discussions).

We welcome contributions, feature requests, and suggestions. If you encounter any problems,
please open an [issue](https://github.com/MineralsCloud/CrystallographyBase.jl/issues).
The [Contributing](@ref) page has
a few guidelines that should be followed when opening pull requests and contributing code.

## Manual outline

```@contents
Pages = [
    "man/installation.md",
    "man/troubleshooting.md",
    "man/definitions.md",
    "man/examples.md",
    "developers/contributing.md",
    "developers/style-guide.md",
    "developers/design-principles.md",
]
Depth = 3
```

## Library outline

```@contents
Pages = ["lib/public.md", "lib/internals.md"]
```

### [Index](@id main-index)

```@index
Pages = ["lib/public.md"]
```
