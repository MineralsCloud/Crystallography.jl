module Crystallography

using CoordinateTransformations: IdentityTransformation
using LinearAlgebra: Diagonal, I, cross, det
using Spglib: Cell
using StaticArrays: SVector, SMatrix, SHermitianCompact
using Unitful: AbstractQuantity, ustrip, unit

import LinearAlgebra: dot, norm
import Spglib: basis_vectors

include("lattice.jl")

Base.size(::Lattice) = (3, 3)
Base.length(::Lattice) = 9  # Number of elements
Base.getindex(A::Lattice, i::Integer, j::Integer) = getindex(A.data, i, j)
Base.eltype(::Lattice{T}) where {T} = T

include("arithmetics.jl")
include("Symmetry.jl")
include("transform.jl")

end
