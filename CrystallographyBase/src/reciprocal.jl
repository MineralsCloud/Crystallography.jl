using Compat: eachslice
using Counters: counter
using LinearAlgebra: cross
using Spglib: get_ir_reciprocal_mesh

export ReciprocalPoint, ReciprocalLattice, reciprocal_mesh, coordinates, weights

struct ReciprocalLattice{T} <: AbstractLattice{T}
    data::SMatrix{3,3,T,9}
end
function ReciprocalLattice(lattice::Lattice)
    Ω = cellvolume(lattice)
    𝐚, 𝐛, 𝐜 = basis_vectors(lattice)
    return ReciprocalLattice(
        inv(Ω) * transpose(hcat(cross(𝐛, 𝐜), cross(𝐜, 𝐚), cross(𝐚, 𝐛))),
    )
end
@functor ReciprocalLattice

Base.inv(lattice::Lattice) = ReciprocalLattice(lattice)
function Base.inv(lattice::ReciprocalLattice)
    Ω⁻¹ = cellvolume(lattice)
    𝐚⁻¹, 𝐛⁻¹, 𝐜⁻¹ = basis_vectors(lattice)
    return Lattice(inv(Ω⁻¹) * hcat(cross(𝐛⁻¹, 𝐜⁻¹), cross(𝐜⁻¹, 𝐚⁻¹), cross(𝐚⁻¹, 𝐛⁻¹)))
end

basis_vectors(lattice::ReciprocalLattice) = lattice[1, :], lattice[2, :], lattice[3, :]

"""
    ReciprocalPoint(x, y, z, w)

Represent a special point of the 3D Brillouin zone. Each of them has a weight `w`.
"""
struct ReciprocalPoint{T}
    coord::SVector{3,T}
    weight::Float64
end
ReciprocalPoint(coord::AbstractVector{T}, weight) where {T} =
    ReciprocalPoint{T}(SVector{3}(coord), weight)
ReciprocalPoint(x, y, z, w) = ReciprocalPoint(SVector(x, y, z), w)
@functor ReciprocalPoint (coord,)

# See example in https://spglib.github.io/spglib/python-spglib.html#get-ir-reciprocal-mesh
function reciprocal_mesh(
    cell::Cell,
    mesh,
    is_shift = falses(3);
    is_time_reversal = true,
    symprec = 1e-5,
    cartesian = false,
    ir_only = true,
)
    _, mapping, grid = get_ir_reciprocal_mesh(
        cell,
        mesh,
        is_shift;
        is_time_reversal = is_time_reversal,
        symprec = symprec,
    )
    shift = is_shift ./ 2  # true / 2 = 0.5, false / 2 = 0
    mapping = convert(Vector{Int}, mapping)
    weights = counter(mapping)
    total_number = length(mapping)  # Number of all k-points, not only the irreducible ones
    crystal_coord = if ir_only
        map(unique(mapping)) do index
            x, y, z = (grid[:, index] .+ shift) ./ mesh
            weight = weights[index] / total_number
            ReciprocalPoint(x, y, z, weight)
        end
    else
        map(eachslice(grid; dims = 2)) do point
            x, y, z = (point .+ shift) ./ mesh
            weight = 1 / total_number
            ReciprocalPoint(x, y, z, weight)
        end
    end
    if cartesian
        t = CartesianFromFractional(inv(Lattice(cell)))
        return map(t, crystal_coord)
    else
        return crystal_coord
    end
end

coordinates(arr::AbstractArray{<:ReciprocalPoint}) = map(x -> x.coord, arr)

weights(arr::AbstractArray{<:ReciprocalPoint}) = map(x -> x.weight, arr)
