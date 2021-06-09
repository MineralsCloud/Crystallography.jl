using Counters: counter
using LinearAlgebra: cross
using Spglib: get_ir_reciprocal_mesh

export ReciprocalPoint, ReciprocalLattice, reciprocal_mesh, coordinates, weights

struct ReciprocalLattice{T}
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
        map(unique(mapping)) do id
            x, y, z = (grid[:, id+1] .+ shift) ./ mesh  # Add 1 because `mapping` index starts from 0
            weight = weights[id] / total_number  # Should use `id` not `id + 1`!
            ReciprocalPoint(x, y, z, weight)
        end
    else
        mapslices(grid; dims = 1) do point
            x, y, z = (point .+ shift) ./ mesh  # Add 1 because `mapping` index starts from 0
            weight = 1 / total_number
            ReciprocalPoint(x, y, z, weight)
        end |> vec
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

Base.iterate(lattice::ReciprocalLattice) = iterate(lattice.data)
Base.iterate(lattice::ReciprocalLattice, state) = iterate(lattice.data, state)

Base.eltype(::ReciprocalLattice{T}) where {T} = T

Base.length(::ReciprocalLattice) = 9

Base.size(::ReciprocalLattice) = (3, 3)
Base.size(::ReciprocalLattice, dim::Integer) = dim <= 2 ? 3 : 1

Base.IteratorSize(::Type{<:ReciprocalLattice}) = Base.HasShape{2}()

Base.axes(lattice::ReciprocalLattice, dim::Integer) = axes(lattice.data, dim)

Base.getindex(lattice::ReciprocalLattice, i) = getindex(lattice.data, i)
Base.getindex(lattice::ReciprocalLattice, I::Vararg) = getindex(lattice.data, I...)

Base.firstindex(::ReciprocalLattice) = 1

Base.lastindex(::ReciprocalLattice) = 9
