# Crystallography

|                                 **Documentation**                                  |                                                                                                 **Build Status**                                                                                                 |                                        **Others**                                         |
| :--------------------------------------------------------------------------------: | :--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------: | :---------------------------------------------------------------------------------------: |
| [![Stable][docs-stable-img]][docs-stable-url] [![Dev][docs-dev-img]][docs-dev-url] | [![Build Status][gha-img]][gha-url] [![Build Status][appveyor-img]][appveyor-url] [![Build Status][cirrus-img]][cirrus-url] [![pipeline status][gitlab-img]][gitlab-url] [![Coverage][codecov-img]][codecov-url] | [![GitHub license][license-img]][license-url] [![Code Style: Blue][style-img]][style-url] |

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://MineralsCloud.github.io/Crystallography.jl/stable
[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://MineralsCloud.github.io/Crystallography.jl/dev
[gha-img]: https://github.com/MineralsCloud/Crystallography.jl/workflows/CI/badge.svg
[gha-url]: https://github.com/MineralsCloud/Crystallography.jl/actions
[appveyor-img]: https://ci.appveyor.com/api/projects/status/github/MineralsCloud/Crystallography.jl?svg=true
[appveyor-url]: https://ci.appveyor.com/project/singularitti/Crystallography-jl
[cirrus-img]: https://api.cirrus-ci.com/github/MineralsCloud/Crystallography.jl.svg
[cirrus-url]: https://cirrus-ci.com/github/MineralsCloud/Crystallography.jl
[gitlab-img]: https://gitlab.com/singularitti/Crystallography.jl/badges/main/pipeline.svg
[gitlab-url]: https://gitlab.com/singularitti/Crystallography.jl/-/pipelines
[codecov-img]: https://codecov.io/gh/MineralsCloud/Crystallography.jl/branch/main/graph/badge.svg
[codecov-url]: https://codecov.io/gh/MineralsCloud/Crystallography.jl
[license-img]: https://img.shields.io/github/license/MineralsCloud/Crystallography.jl
[license-url]: https://github.com/MineralsCloud/Crystallography.jl/blob/main/LICENSE
[style-img]: https://img.shields.io/badge/code%20style-blue-4495d1.svg
[style-url]: https://github.com/invenia/BlueStyle

The code, which is [hosted on GitHub](https://github.com/MineralsCloud/Crystallography.jl), is tested
using various continuous integration services for its validity.

This repository is created and maintained by
[@singularitti](https://github.com/singularitti), and contributions are highly welcome.

## Package features



## Installation

The package can be installed with the Julia package manager.
From [the Julia REPL](https://docs.julialang.org/en/v1/stdlib/REPL/), type `]` to enter
the [Pkg mode](https://docs.julialang.org/en/v1/stdlib/REPL/#Pkg-mode) and run:

```julia-repl
pkg> add Crystallography
```

Or, equivalently, via [`Pkg.jl`](https://pkgdocs.julialang.org/v1/):

```julia
julia> import Pkg; Pkg.add("Crystallography")
```

## Documentation

- [**STABLE**][docs-stable-url] — **documentation of the most recently tagged version.**
- [**DEV**][docs-dev-url] — _documentation of the in-development version._

## Project status

The package is developed for and tested against Julia `v1.6` and above on Linux, macOS, and
Windows.

## Questions and contributions

You can post usage questions on
[our discussion page](https://github.com/MineralsCloud/Crystallography.jl/discussions).

We welcome contributions, feature requests, and suggestions. If you encounter any problems,
please open an [issue](https://github.com/MineralsCloud/Crystallography.jl/issues).
The [Contributing](@ref) page has
a few guidelines that should be followed when opening pull requests and contributing code.
