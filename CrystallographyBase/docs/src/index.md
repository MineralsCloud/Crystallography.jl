```@meta
CurrentModule = CrystallographyBase
```

# CrystallographyBase

Documentation for [CrystallographyBase](https://github.com/MineralsCloud/CrystallographyBase.jl).

See the [Index](@ref main-index) for the complete list of documented functions
and types.

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
  - [`StructHelpers.jl`](https://github.com/jw3126/StructHelpers.jl) `v0.1` and above
- OS: macOS, Linux, Windows, and FreeBSD
- Architecture: x86, x64, ARM

## Documentation

- [**STABLE**](https://MineralsCloud.github.io/CrystallographyBase.jl/stable) &mdash; **documentation of the most recently tagged version.**
- [**DEV**](https://MineralsCloud.github.io/CrystallographyBase.jl/dev) &mdash; _documentation of the in-development version._

## Project Status

The package is tested against, and being developed for, Julia `1.6` and above on Linux,
macOS, and Windows.

## Questions and Contributions

Usage questions can be posted on
[our discussion page](https://github.com/MineralsCloud/CrystallographyBase.jl/discussions).

Contributions are very welcome, as are feature requests and suggestions. Please open an
[issue](https://github.com/MineralsCloud/CrystallographyBase.jl/issues)
if you encounter any problems. The [contributing](@ref) page has
a few guidelines that should be followed when opening pull requests and contributing code.

## Manual Outline

```@contents
Pages = [
    "installation.md",
    "contributing.md",
    "troubleshooting.md",
]
Depth = 3
```

## Library Outline

```@contents
Pages = ["public.md"]
```

### [Index](@id main-index)

```@index
Pages = ["public.md"]
```
