module Symmetry

using LinearAlgebra: I, diagm, det

using CoordinateTransformations: AffineMap, Translation, LinearMap
using LibSymspg: get_symmetry, get_spacegroup, ir_reciprocal_mesh
using StaticArrays: SMatrix, SDiagonal

using Crystallography: CrystalCoordinates, Cell

export SeitzOperator, IdentityOperator
export getsymmetry,
    getspacegroup, irreciprocalmesh, isidentity, istranslation, ispointsymmetry

function getsymmetry(cell::Cell, symprec::AbstractFloat = 1e-5; seitz::Bool = false)
    maps, translations = get_symmetry(
        cell.lattice,
        cell.positions,
        cell.numbers,
        length(cell.numbers),
        symprec,
    )
    return if seitz
        (SeitzOperator(LinearMap(m), Translation(t)) for (m, t) in zip(maps, translations))
    else
        (AffineMap(m, t) for (m, t) in zip(maps, translations))
    end
end

function getspacegroup(cell::Cell, symprec::AbstractFloat = 1e-5)
    return get_spacegroup(
        cell.lattice,
        cell.positions,
        cell.numbers,
        length(cell.numbers),
        symprec,
    )
end # function getspacegroup

function irreciprocalmesh(
    cell::Cell,
    mesh::AbstractVector{Int},
    symprec::AbstractFloat = 1e-5;
    is_shift::AbstractVector{Bool} = falses(3),
    is_time_reversal::Bool = false,
)
    return ir_reciprocal_mesh(
        mesh,
        collect(is_shift),
        is_time_reversal,
        cell.lattice,
        cell.positions,
        cell.numbers,
        length(cell.numbers),
        symprec,
    )
end # function irreciprocalmesh

struct SeitzOperator{T}
    m::SMatrix{4,4,T}
end
SeitzOperator(m::AbstractMatrix) = SeitzOperator(SMatrix{4,4}(m))
function SeitzOperator(l::LinearMap)
    m = l.linear
    @assert size(m) == (3, 3)
    x = diagm(ones(eltype(m), 4))
    x[1:3, 1:3] = m
    return SeitzOperator(x)
end # function PointSymmetryOperator
function SeitzOperator(t::Translation)
    v = t.translation
    @assert length(v) == 3
    x = diagm(ones(eltype(v), 4))
    x[1:3, 4] = v
    return SeitzOperator(x)
end # function TranslationOperator
SeitzOperator(m::LinearMap, t::Translation) =
    SeitzOperator(t) * SeitzOperator(m) * SeitzOperator(inv(t))
SeitzOperator(a::AffineMap) = SeitzOperator(LinearMap(a.linear), Translation(a.translation))
function SeitzOperator(s::SeitzOperator, pos::AbstractVector)
    @assert length(pos) == 3
    t = SeitzOperator(Translation(pos))
    return t * s * inv(t)
end # function SeitzOperator
const IdentityOperator = SeitzOperator(SDiagonal(1, 1, 1, 1))

isidentity(op::SeitzOperator) = op.m == I

function istranslation(op::SeitzOperator)
    m = op.m
    if m[1:3, 1:3] != I || !(iszero(m[4, 1:3]) && isone(m[4, 4]))
        return false
    end
    return true
end # function istranslation

function ispointsymmetry(op::SeitzOperator)
    m = op.m
    if !(iszero(m[4, 1:3]) && iszero(m[1:3, 4]) && isone(m[4, 4]) && abs(det(m[1:3, 1:3])) == 1)
        return false
    end
    return true
end # function ispointsymmetry

Base.getindex(A::SeitzOperator, I::Vararg{Int}) = getindex(A.m, I...)

Base.:*(m::SeitzOperator, c::CrystalCoordinates) = CrystalCoordinates((m.m*[c; 1])[1:3])
Base.:*(a::SeitzOperator, b::SeitzOperator) = SeitzOperator(a.m * b.m)

function Base.convert(::Type{Translation}, op::SeitzOperator)
    @assert(istranslation(op), "operator is not a translation!")
    return Translation(collect(op.m[1:3, 4]))
end # function Base.convert
function Base.convert(::Type{LinearMap}, op::SeitzOperator)
    @assert(ispointsymmetry(op), "operator is not a point symmetry!")
    return LinearMap(collect(op.m[1:3, 1:3]))
end # function Base.convert

Base.inv(op::SeitzOperator) = SeitzOperator(Base.inv(op.m))

end # module Symmetry
