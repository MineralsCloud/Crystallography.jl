module CrystallographyBase

using Functors: @functor
using LinearAlgebra: dot
using Requires: @require
using StaticArrays: SVector, SMatrix
using StructEquality: @struct_hash_equal_isequal_isapprox

function __init__()
    @require PeriodicTable = "7b2266bf-644c-5ea3-82d8-af4bbd25a884" begin
        @eval using PeriodicTable: Element, elements
        include("crystaldensity.jl")
    end
    @require Mendeleev = "c116f080-063d-490a-9873-2b5b2cce4c34" begin
        @eval using Mendeleev: Element, elements
        include("crystaldensity.jl")
    end
end

# `LinearAlgebra.det` is much slower and more inaccurate than my simple `_det`.
function _det(matrix)
    a, d, g, b, e, h, c, f, i = matrix  # Only works for 3×3 matrices
    return a * e * i + b * f * g + c * d * h - c * e * g - b * d * i - a * f * h
end

include("systems.jl")
include("lattice.jl")
include("cell.jl")
include("reciprocal.jl")
include("geometry.jl")
# include("wignerseitz.jl")
include("metric.jl")
include("volume.jl")
include("transform.jl")

end
