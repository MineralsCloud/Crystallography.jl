export RealSpace, ReciprocalSpace

abstract type AbstractSpace end
struct RealSpace <: AbstractSpace end
struct ReciprocalSpace <: AbstractSpace end

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

include("miller.jl")
include("metric.jl")
include("reciprocal.jl")
include("cartesian.jl")
