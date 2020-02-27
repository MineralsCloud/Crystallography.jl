module Symmetry

using LinearAlgebra: diagm, I

using CoordinateTransformations: AffineMap, Translation, LinearMap
using LibSymspg: get_symmetry
using StaticArrays: SMatrix

using Crystallography: CrystalCoordinates, Cell

export SeitzOperator
export getsymmetry, isidentity, istranslation, ispointsymmetry

function getsymmetry(cell::Cell, symprec::AbstractFloat = 1e-5)
    maps, translations = get_symmetry(
        cell.lattice,
        cell.positions,
        cell.numbers,
        length(cell.numbers),
        symprec,
    )
    return (AffineMap(m, t) for (m, t) in zip(maps, translations))
end

struct SeitzOperator{T} <: AbstractMatrix{T}
    m::SMatrix{4,4,T}
end
SeitzOperator(m::AbstractMatrix) = SeitzOperator(SMatrix{4,4}(m))
SeitzOperator() = SeitzOperator(ones(Int, 4, 4))
function SeitzOperator(m::LinearMap)
    @assert size(m.linear) == (3, 3)
    x = diagm(ones(eltype(m.linear), 4))
    x[1:3, 1:3] = m.linear
    return SeitzOperator(x)
end # function PointSymmetryOperator
function SeitzOperator(t::Translation)
    @assert length(t.translation) == 3
    x = diagm(ones(eltype(t.translation), 4))
    x[1:3, 4] = t.translation
    return SeitzOperator(x)
end # function TranslationOperator
SeitzOperator(m::LinearMap, t::Translation) = SeitzOperator(t) * SeitzOperator(m) * SeitzOperator(inv(t))
SeitzOperator(a::AffineMap) = SeitzOperator(LinearMap(a.linear), Translation(a.translation))
function SeitzOperator(s::SeitzOperator, pos::AbstractVector)
    @assert length(pos) == 3
    t = SeitzOperator(Translation(pos))
    return t * s * inv(t)
end # function SeitzOperator

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
    if !(iszero(m[4, 1:3]) && iszero(m[1:3, 4]) && isone(m[4, 4]))
        return false
    end
    return true
end # function ispointsymmetry

Base.size(::SeitzOperator) = (4, 4)

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

end # module Symmetry
