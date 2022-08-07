# CrystallographyBase

|                                 **Documentation**                                  |                                                                                                 **Build Status**                                                                                                 |                  **LICENSE**                  |
| :--------------------------------------------------------------------------------: | :--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------: | :-------------------------------------------: |
| [![Stable][docs-stable-img]][docs-stable-url] [![Dev][docs-dev-img]][docs-dev-url] | [![Build Status][gha-img]][gha-url] [![Build Status][appveyor-img]][appveyor-url] [![Build Status][cirrus-img]][cirrus-url] [![pipeline status][gitlab-img]][gitlab-url] [![Coverage][codecov-img]][codecov-url] | [![GitHub license][license-img]][license-url] |

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
[gitlab-img]: https://gitlab.com/singularitti/CrystallographyBase.jl/badges/master/pipeline.svg
[gitlab-url]: https://gitlab.com/singularitti/CrystallographyBase.jl/-/pipelines
[codecov-img]: https://codecov.io/gh/MineralsCloud/CrystallographyBase.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/MineralsCloud/CrystallographyBase.jl
[license-img]: https://img.shields.io/github/license/MineralsCloud/CrystallographyBase.jl
[license-url]: https://github.com/MineralsCloud/CrystallographyBase.jl/blob/master/LICENSE

This package provides some basic types and methods for crystallography calculations.
For more features, see [`Crystallography.jl`](https://github.com/MineralsCloud/Crystallography.jl).

The code is [hosted on GitHub](https://github.com/MineralsCloud/CrystallographyBase.jl),
with some continuous integration services to test its validity.

This repository is created and maintained by [singularitti](https://github.com/singularitti).
You are very welcome to contribute.

## Installation

The package can be installed with the Julia package manager.
From the Julia REPL, type `]` to enter the Pkg REPL mode and run:

```
pkg> add CrystallographyBase
```

Or, equivalently, via the `Pkg` API:

```julia
julia> import Pkg; Pkg.add("CrystallographyBase")
```

## Compatibility

- [Julia version: `v1.3` to `v1.7`](https://julialang.org/downloads/)
- Dependencies:
  - [`Combinatorics.jl`](https://github.com/JuliaMath/Combinatorics.jl) `v0.7.0` and above
  - [`CoordinateTransformations.jl`](https://github.com/JuliaGeometry/CoordinateTransformations.jl) `v0.5.1` and above
  - [`Counters.jl`](https://github.com/scheinerman/Counters.jl) `v0.3.0` and above
  - [`EnumX.jl`](https://github.com/fredrikekre/EnumX.jl) `v1.0.0` and above
  - [`Functors.jl`](https://github.com/FluxML/Functors.jl) `v0.1.0` and above
  - [`PeriodicTable.jl`](https://github.com/JuliaPhysics/PeriodicTable.jl) `v0.1.0` and above
  - [`StaticArrays.jl`](https://github.com/JuliaArrays/StaticArrays.jl) `v0.8.3` and above
- OS: macOS, Linux, Windows, and FreeBSD
- Architecture: x86, x64, ARM

## Documentation

- [**STABLE**][docs-stable-url] &mdash; **documentation of the most recently tagged version.**
- [**DEV**][docs-dev-url] &mdash; _documentation of the in-development version._

## Project Status

The package is tested against, and being developed for, Julia `1.6` and above on Linux,
macOS, and Windows.

## Questions and Contributions

Usage questions can be posted on [our discussion page][discussions-url].

Contributions are very welcome, as are feature requests and suggestions. Please open an
[issue][issues-url] if you encounter any problems. The [contributing](@ref) page has
a few guidelines that should be followed when opening pull requests and contributing code.

[discussions-url]: https://github.com/MineralsCloud/CrystallographyBase.jl/discussions
[issues-url]: https://github.com/MineralsCloud/CrystallographyBase.jl/issues
[contrib-url]: https://github.com/MineralsCloud/CrystallographyBase.jl/discussions
