module CrystallographyBase

using LinearAlgebra: dot
using StaticArrays: SVector, MMatrix
using StructEquality: @struct_hash_equal_isequal_isapprox

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
include("metric.jl")
include("volume.jl")
include("transform.jl")
include("show.jl")

end
