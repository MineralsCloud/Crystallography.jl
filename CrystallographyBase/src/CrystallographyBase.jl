module CrystallographyBase

using Functors: @functor
using LinearAlgebra: det
using StaticArrays: SVector, SMatrix

import LinearAlgebra: dot, norm

include("lattice.jl")
include("cell.jl")
include("reciprocal.jl")
# include("wignerseitz.jl")
include("miller.jl")
include("metric.jl")
include("volume.jl")
include("transform.jl")

end
