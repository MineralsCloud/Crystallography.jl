# CrystallographyBase

|                                 **Documentation**                                  |                                                                                                 **Build Status**                                                                                                 |                                        **Others**                                         |
| :--------------------------------------------------------------------------------: | :--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------: | :---------------------------------------------------------------------------------------: |
| [![Stable][docs-stable-img]][docs-stable-url] [![Dev][docs-dev-img]][docs-dev-url] | [![Build Status][gha-img]][gha-url] [![Build Status][appveyor-img]][appveyor-url] [![Build Status][cirrus-img]][cirrus-url] [![pipeline status][gitlab-img]][gitlab-url] [![Coverage][codecov-img]][codecov-url] | [![GitHub license][license-img]][license-url] [![Code Style: Blue][style-img]][style-url] |

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://MineralsCloud.github.io/CrystallographyBase.jl/stable
[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://MineralsCloud.github.io/CrystallographyBase.jl/dev
[gha-img]: https://github.com/MineralsCloud/CrystallographyBase.jl/workflows/CI/badge.svg
[gha-url]: https://github.com/MineralsCloud/CrystallographyBase.jl/actions
[appveyor-img]: https://ci.appveyor.com/api/projects/status/github/MineralsCloud/CrystallographyBase.jl?svg=true
[appveyor-url]: https://ci.appveyor.com/project/singularitti/CrystallographyBase-jl
[cirrus-img]: https://api.cirrus-ci.com/github/MineralsCloud/CrystallographyBase.jl.svg
[cirrus-url]: https://cirrus-ci.com/github/MineralsCloud/CrystallographyBase.jl
[gitlab-img]: https://gitlab.com/singularitti/CrystallographyBase.jl/badges/main/pipeline.svg
[gitlab-url]: https://gitlab.com/singularitti/CrystallographyBase.jl/-/pipelines
[codecov-img]: https://codecov.io/gh/MineralsCloud/CrystallographyBase.jl/branch/main/graph/badge.svg
[codecov-url]: https://codecov.io/gh/MineralsCloud/CrystallographyBase.jl
[license-img]: https://img.shields.io/github/license/MineralsCloud/CrystallographyBase.jl
[license-url]: https://github.com/MineralsCloud/CrystallographyBase.jl/blob/main/LICENSE
[style-img]: https://img.shields.io/badge/code%20style-blue-4495d1.svg
[style-url]: https://github.com/invenia/BlueStyle

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

```julia
julia> import Pkg; Pkg.add("CrystallographyBase")
```

## Documentation

- [**STABLE**][docs-stable-url] — **documentation of the most recently tagged version.**
- [**DEV**][docs-dev-url] — _documentation of the in-development version._

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
