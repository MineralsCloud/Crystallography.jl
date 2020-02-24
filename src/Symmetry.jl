module Symmetry

using LinearAlgebra: diagm, I

using CoordinateTransformations: AffineMap, Translation, LinearMap
using LibSymspg: get_symmetry

using Crystallography: CrystalCoordinates, Cell

export SeitzOperator, IdentityOperator, TranslationOperator, PointSymmetryOperator
export getsymmetry

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

struct SeitzOperator{T}
    m::T
    function SeitzOperator{T}(m) where {T}
        @assert(size(m) == (4, 4), "The operator must be of size 4x4!")
        return new(m)
    end
end
SeitzOperator(m::T) where {T} = SeitzOperator{T}(m)
function SeitzOperator(m::AffineMap)
    linear, translation = m.linear, m.translation
    return SeitzOperator(translation ∘ linear ∘ inv(translation))
end # function SeitzOperator
function SeitzOperator(s::SeitzOperator, pos::AbstractVector)
    @assert length(pos) == 3
    t = TranslationOperator(Translation(pos))
    return t * s * inv(t)
end # function SeitzOperator
SeitzOperator(a:AbstractMatrix, pos::AbstractVector) = SeitzOperator(SeitzOperator(a), pos)

IdentityOperator() = SeitzOperator(ones(Int, 4, 4))

function TranslationOperator(v::AbstractVector)
    @assert length(v) == 3
    x = diagm(ones(eltype(v), 4))
    x[1:3, 4] = v
    return SeitzOperator(x)
end # function TranslationOperator
TranslationOperator(t::Translation) = TranslationOperator(t.translation)

function PointSymmetryOperator(m::AbstractMatrix)
    @assert size(m) == (3, 3)
    x = diagm(ones(eltype(m), 4))
    x[1:3, 1:3] = m
    return SeitzOperator(x)
end # function PointSymmetryOperator
PointSymmetryOperator(l::LinearMap) = PointSymmetryOperator(l.linear)

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

# Base.eltype(t::Translation) = eltype(t.translation)
# Base.eltype(m::LinearMap) = eltype(m.linear)

Base.:*(m::SeitzOperator, c::CrystalCoordinates) = CrystalCoordinates((m.m*[c; 1])[1:3])
Base.:*(a::SeitzOperator, b::SeitzOperator) = SeitzOperator(a.m * b.m)

    # Base.convert(::Type{Translation}, op::TranslationOperator) = Translation(collect(op.m[1:3, 4]))

    # Base.inv(op::TranslationOperator) = TranslationOperator(inv(convert(Translation, op)))


end # module Symmetry
