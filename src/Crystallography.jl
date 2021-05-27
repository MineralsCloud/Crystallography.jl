module Crystallography

using CoordinateTransformations: IdentityTransformation
using LinearAlgebra: Diagonal, I, cross, det
using Spglib: Cell
using StaticArrays: SVector, SMatrix, SHermitianCompact
using Unitful: AbstractQuantity, ustrip, unit

import LinearAlgebra: dot, norm
import Spglib: basis_vectors

include("lattice.jl")
include("arithmetics.jl")
include("Symmetry.jl")
include("transform.jl")

end
