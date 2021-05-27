module Crystallography

using LinearAlgebra: det
using StaticArrays: SVector, SMatrix

import LinearAlgebra: dot, norm
import Spglib: basis_vectors

include("lattice.jl")
include("arithmetics.jl")
include("Symmetry.jl")
include("transform.jl")

end
