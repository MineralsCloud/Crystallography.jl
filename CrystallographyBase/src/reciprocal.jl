using Brillouin: irrfbz_path
using Compat: eachslice
using Counters: counter
using LinearAlgebra: cross
using Spglib: get_ir_reciprocal_mesh, get_spacegroup_type

export ReciprocalPoint,
    ReciprocalLattice, ReciprocalPath, reciprocal_mesh, coordinates, weights

"""
    ReciprocalLattice(mat::SMatrix)

Construct a `ReciprocalLattice`.

!!! warning
    You should not use this function directly, always use `inv` of a `Lattice`.
"""
struct ReciprocalLattice{T} <: AbstractLattice{T}
    data::SMatrix{3,3,T,9}
end
@functor ReciprocalLattice

"""
    inv(lattice::Lattice)
    inv(lattice::ReciprocalLattice)

Get the reciprocal of a `Lattice` or a `ReciprocalLattice`.
"""
function Base.inv(lattice::Lattice)
    Ω = det(lattice.data)  # Cannot use `cellvolume`, it takes the absolute value!
    𝐚, 𝐛, 𝐜 = basis_vectors(lattice)
    return ReciprocalLattice(
        inv(Ω) * transpose(hcat(cross(𝐛, 𝐜), cross(𝐜, 𝐚), cross(𝐚, 𝐛))),
    )
end
function Base.inv(lattice::ReciprocalLattice)
    Ω⁻¹ = det(lattice.data)  # Cannot use `cellvolume`, it takes the absolute value!
    𝐚⁻¹, 𝐛⁻¹, 𝐜⁻¹ = basis_vectors(lattice)
    return Lattice(inv(Ω⁻¹) * hcat(cross(𝐛⁻¹, 𝐜⁻¹), cross(𝐜⁻¹, 𝐚⁻¹), cross(𝐚⁻¹, 𝐛⁻¹)))
end

"""
    basis_vectors(lattice::ReciprocalLattice)

Get the three basis vectors from a `ReciprocalLattice`.
"""
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
"""
    reciprocal_mesh(cell::Cell, mesh, is_shift; kwargs...)

List the `ReciprocalPoint`s from the mesh of the reciprocal space of a `Cell`.

# Arguments
- `cell::Cell`: the cell.
- `mesh`: a vector of three integers which specify the mesh numbers along reciprocal primitive axis.
- `is_shift=falses(3)`: a vector of three elements specifying whether the mesh is shifted along the axis in half of adjacent mesh points irrespective of the mesh numbers. The allowed values are `0`, `1`, `true`, and `false`.
- `is_time_reversal=true`: whether to impose the time reversal symmetry on the mesh.
- `symprec=1e-5`: distance tolerance in Cartesian coordinates to find crystal symmetry.
- `cartesian=false`: whether to return the reciprocal points in Cartesian coordinates.
- `ir_only=true`: whether to return the reciprocal points only in the irreducible Brillouin zone.
"""
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

"""
    coordinates(arr::AbstractArray{<:ReciprocalPoint})

Get the coordinates of an array of `ReciprocalPoint`s.
"""
coordinates(arr::AbstractArray{<:ReciprocalPoint}) = map(x -> x.coord, arr)

"""
    weights(arr::AbstractArray{<:ReciprocalPoint})

Get the weights of an array of `ReciprocalPoint`s.
"""
weights(arr::AbstractArray{<:ReciprocalPoint}) = map(x -> x.weight, arr)

function Base.show(io::IO, x::ReciprocalPoint)
    if get(io, :compact, false) || get(io, :typeinfo, nothing) == typeof(x)
        Base.show_default(IOContext(io, :limit => true), x)  # From https://github.com/mauro3/Parameters.jl/blob/ecbf8df/src/Parameters.jl#L556
    else
        println(io, string(typeof(x)))
        print(io, " coord = ", x.coord, ", weight = ", x.weight)
    end
end

struct ReciprocalPath{N}
    special_points::Dict{Symbol,SVector{N,Float64}}
    suggested_paths::Vector{Vector{Symbol}}
    lattice::Lattice
end
"""
    ReciprocalPath(lattice::Lattice, spgnum::Integer)
    ReciprocalPath(cell::Cell)

Construct a `ReciprocalPath` from a `Lattice` or a `Cell`.
"""
function ReciprocalPath(lattice::Lattice, spgnum::Integer)
    kpath = irrfbz_path(spgnum, collect(basis_vectors(lattice)))
    return ReciprocalPath(kpath.points, kpath.paths, lattice)
end
function ReciprocalPath(cell::Cell)
    spg = get_spacegroup_type(cell)
    return ReciprocalPath(Lattice(cell), spg.number)
end

"""
    coordinates(path::ReciprocalPath, cartesian = false)

If `cartesian` is `true`, return the coordinates in the Cartesian coordinate system.
"""
coordinates(path::ReciprocalPath, cartesian = false) =
    cartesian ?
    Dict(
        key => CartesianFromFractional(inv(path.lattice))(value) for
        (key, value) in path.special_points
    ) : path.special_points
