using LinearAlgebra: cross

export ReciprocalPoint,
    ReciprocalLattice, MonkhorstPackGrid, reciprocal, coordinates, weights

"""
    ReciprocalLattice(mat::SMatrix)

Construct a `ReciprocalLattice`.

!!! warning
    You should not use this function directly, always use `reciprocal` of a `Lattice`.
"""
@struct_hash_equal_isequal_isapprox struct ReciprocalLattice{T} <: AbstractLattice{T}
    data::SMatrix{3,3,T,9}
end
@functor ReciprocalLattice

"""
    reciprocal(lattice::Lattice)
    reciprocal(lattice::ReciprocalLattice)

Get the reciprocal of a `Lattice` or a `ReciprocalLattice`.
"""
function reciprocal(lattice::Lattice)
    Ω = _det(lattice.data)  # Cannot use `cellvolume`, it takes the absolute value!
    𝐚, 𝐛, 𝐜 = basisvectors(lattice)
    return ReciprocalLattice(
        inv(Ω) * transpose(hcat(cross(𝐛, 𝐜), cross(𝐜, 𝐚), cross(𝐚, 𝐛)))
    )
end
function reciprocal(lattice::ReciprocalLattice)
    Ω⁻¹ = _det(lattice.data)  # Cannot use `cellvolume`, it takes the absolute value!
    𝐚⁻¹, 𝐛⁻¹, 𝐜⁻¹ = basisvectors(lattice)
    return Lattice(inv(Ω⁻¹) * hcat(cross(𝐛⁻¹, 𝐜⁻¹), cross(𝐜⁻¹, 𝐚⁻¹), cross(𝐚⁻¹, 𝐛⁻¹)))
end

"""
    basisvectors(lattice::ReciprocalLattice)

Get the three basis vectors from a `ReciprocalLattice`.
"""
basisvectors(lattice::ReciprocalLattice) = lattice[1, :], lattice[2, :], lattice[3, :]

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

"""
    MonkhorstPackGrid(mesh, is_shift)

Represent the Monkhorst--Pack grid.

# Arguments
- `mesh`: A length-three vector specifying the k-point grid (``nk_1 × nk_2 × nk_3``) as in Monkhorst--Pack grids.
- `is_shift`: A length-three vector specifying whether the grid is displaced by half a grid step in the corresponding directions.
"""
struct MonkhorstPackGrid
    mesh::SVector{3,UInt}
    is_shift::SVector{3,Bool}
    function MonkhorstPackGrid(mesh, is_shift)
        @assert all(mesh .>= 1)
        if eltype(is_shift) != Bool
            is_shift = Bool.(is_shift)
        end
        return new(mesh, is_shift)
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
