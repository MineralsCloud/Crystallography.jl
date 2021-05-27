module Crystallography

using CoordinateTransformations: IdentityTransformation
using LinearAlgebra: Diagonal, I, cross, det
using Spglib: Cell
using StaticArrays: SVector, SMatrix, SHermitianCompact
using Unitful: AbstractQuantity, ustrip, unit

import LinearAlgebra: dot, norm
import Spglib: basis_vectors

include("lattice.jl")

"""
    cellvolume(a, b, c, α, β, γ)

Calculates the cell volume from 6 cell parameters.
"""
cellvolume(a, b, c, α, β, γ) =
    a * b * c * sqrt(sin(α)^2 - cos(β)^2 - cos(γ)^2 + 2 * cos(α) * cos(β) * cos(γ))
"""
    cellvolume(l::Lattice)
    cellvolume(c::Cell)

Calculates the cell volume from a `Lattice` or a `Cell`.
"""
cellvolume(lattice::Lattice) = abs(det(lattice.data))
cellvolume(cell::Cell) = cellvolume(cell.lattice)

Base.size(::Lattice) = (3, 3)
Base.length(::Lattice) = 9  # Number of elements
Base.getindex(A::Lattice, i::Integer, j::Integer) = getindex(A.data, i, j)
Base.eltype(::Lattice{T}) where {T} = T

include("arithmetics.jl")
include("Symmetry.jl")
include("transform.jl")

end
