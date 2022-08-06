# CrystallographyBase

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://MineralsCloud.github.io/CrystallographyBase.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://MineralsCloud.github.io/CrystallographyBase.jl/dev)
[![Build Status](https://github.com/MineralsCloud/CrystallographyBase.jl/workflows/CI/badge.svg)](https://github.com/MineralsCloud/CrystallographyBase.jl/actions)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/MineralsCloud/CrystallographyBase.jl?svg=true)](https://ci.appveyor.com/project/singularitti/CrystallographyBase-jl)
[![Build Status](https://api.cirrus-ci.com/github/MineralsCloud/CrystallographyBase.jl.svg)](https://cirrus-ci.com/github/MineralsCloud/CrystallographyBase.jl)
[![pipeline status](https://gitlab.com/singularitti/CrystallographyBase.jl/badges/master/pipeline.svg)](https://gitlab.com/singularitti/CrystallographyBase.jl/-/pipelines)
[![Coverage](https://codecov.io/gh/MineralsCloud/CrystallographyBase.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/MineralsCloud/CrystallographyBase.jl)
[![PkgEval](https://JuliaCI.github.io/NanosoldierReports/pkgeval_badges/C/CrystallographyBase.svg)](https://JuliaCI.github.io/NanosoldierReports/pkgeval_badges/report.html)
[![GitHub license](https://img.shields.io/github/license/MineralsCloud/CrystallographyBase.jl)](https://github.com/MineralsCloud/CrystallographyBase.jl/blob/master/LICENSE)

## Package Features

Provides some basic types and methods for crystallography calculations.
For more features, see [`Crystallography.jl`](https://github.com/MineralsCloud/Crystallography.jl).
See [the documentation of the stable version](https://mineralscloud.github.io/CrystallographyBase.jl/stable)
here.

The code is [hosted on GitHub](https://github.com/MineralsCloud/CrystallographyBase.jl),
with some continuous integration services to test its validity.

This repository is created and maintained by [singularitti](https://github.com/singularitti).
You are very welcome to contribute.

## Installation

`CrystallographyBase` is a &nbsp;
<a href="https://julialang.org">
    <img src="https://julialang.org/assets/infra/julia.ico" width="16em">
    Julia Language
</a>
&nbsp; package. To install `CrystallographyBase`,
please <a href="https://docs.julialang.org/en/v1/manual/getting-started/">open
Julia's interactive session (known as REPL)</a> and press <kbd>]</kbd> key in the REPL to use the package mode, then type the following command:

For stable release

```julia
pkg> add CrystallographyBase
```

For current master

```julia
pkg> add CrystallographyBase#master
```

## Compatibility

- [Julia version: `v1.3.0` to `v1.7.2`](https://julialang.org/downloads/)
- Dependencies:
  - [`Combinatorics.jl`](https://github.com/JuliaMath/Combinatorics.jl) `v0.7.0` and above
  - [`CoordinateTransformations.jl`](https://github.com/JuliaGeometry/CoordinateTransformations.jl) `v0.5.1` and above
  - [`Counters.jl`](https://github.com/scheinerman/Counters.jl) `v0.3.0` and above
  - [`Functors.jl`](https://github.com/FluxML/Functors.jl) `v0.1.0` and above
  - [`PeriodicTable.jl`](https://github.com/JuliaPhysics/PeriodicTable.jl) `v0.1.0` and above
  - [`Spglib.jl`](https://github.com/singularitti/Spglib.jl) `v0.2.0` and above
  - [`StaticArrays.jl`](https://github.com/JuliaArrays/StaticArrays.jl) `v0.8.3` and above
- OS: macOS, Linux, Windows, and FreeBSD
- Architecture: x86, x64, ARM
