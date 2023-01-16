module CrystallographyBase

using Functors: @functor
using LinearAlgebra: det
using Requires: @require
using StaticArrays: SVector, SMatrix
using StructEquality: @struct_hash_equal_isequal_isapprox

import LinearAlgebra: dot, norm

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
