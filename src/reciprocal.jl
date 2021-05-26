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
        mat = reciprocal(Lattice(cell.lattice))'
        return map(crystal_coord) do point
            ReciprocalPoint(mat * point.coord, point.weight)
        end
    else
        return crystal_coord
    end
end

coordinates(arr::AbstractArray{<:ReciprocalPoint}) = map(x -> x.coord, arr)

weights(arr::AbstractArray{<:ReciprocalPoint}) = map(x -> x.weight, arr)
