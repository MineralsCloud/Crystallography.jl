module CrystallographyBase

using Functors: @functor
using LinearAlgebra: det
using StaticArrays: SVector, SMatrix
using StructEquality: @struct_hash_equal_isequal_isapprox

import LinearAlgebra: dot, norm

include("systems.jl")
include("lattice.jl")
include("cell.jl")
include("reciprocal.jl")
include("geometry.jl")
# include("wignerseitz.jl")
include("miller.jl")
include("metric.jl")
include("volume.jl")
include("transform.jl")

end
