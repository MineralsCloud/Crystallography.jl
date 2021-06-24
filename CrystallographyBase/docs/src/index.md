```@meta
CurrentModule = CrystallographyBase
```

# CrystallographyBase

Documentation for [CrystallographyBase](https://github.com/MineralsCloud/CrystallographyBase.jl).

## Package Features

Provides some basic types and methods for crystallography calculations.

See the [Index](@ref main-index) for the complete list of documented functions
and types.

The code is [hosted on GitHub](https://github.com/MineralsCloud/CrystallographyBase.jl),
with some continuous integration services to test its validity.

This repository is created and maintained by [singularitti](https://github.com/singularitti).
You are very welcome to contribute.

## Installation
<p>
`CrystallographyBase` is a &nbsp;
    <a href="https://julialang.org">
        <img src="https://julialang.org/favicon.ico" width="16em">
        Julia Language
    </a>
    &nbsp; package. To install XPS,
    please <a href="https://docs.julialang.org/en/v1/manual/getting-started/">open
    Julia's interactive session (known as REPL)</a> and press <kbd>]</kbd> key in the REPL to use the package mode, then type the following command
</p>

For stable release

```julia
pkg> add CrystallographyBase
```

For current master

```julia
pkg> add CrystallographyBase#master
```

## Compatibility

- [Julia version: `v1.3.0` to `v1.6.1`](https://julialang.org/downloads/)
- Dependencies:
  - [`Combinatorics.jl`](https://github.com/JuliaMath/Combinatorics.jl) `v0.7.0` and above
  - [`Compat.jl`](https://github.com/JuliaLang/Compat.jl) `v2.2.0` and above
  - [`CoordinateTransformations.jl`](https://github.com/JuliaGeometry/CoordinateTransformations.jl) `v0.5.1` and above
  - [`Counters.jl`](https://github.com/scheinerman/Counters.jl) `v0.3.0` and above
  - [`Functors.jl`](https://github.com/FluxML/Functors.jl) `v0.1.0` and above
  - [`Spglib.jl`](https://github.com/singularitti/Spglib.jl) `v0.2.0` and above
  - [`StaticArrays.jl`](https://github.com/JuliaArrays/StaticArrays.jl) `v0.8.3` and above
- OS: macOS, Linux, Windows, and FreeBSD
- Architecture: x86, x64, ARM

## Manual Outline

```@contents
Pages = [
    "installation.md",
    "develop.md",
    "api.md",
]
Depth = 3
```

## [Index](@id main-index)

```@index
```
