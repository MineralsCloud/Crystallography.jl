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
# function SeitzOperator(m::Union{TranslationOperator,PointSymmetryOperator}, pos::AbstractVector)
#     @assert length(pos) == 3
#     translation = TranslationOperator(Translation(pos))
#     return translation ∘ m ∘ inv(translation)
# end # function SeitzOperator

IdentityOperator() = SeitzOperator(ones(Int, 4, 4))

function TranslationOperator(t::Translation)
    T = eltype(t)
    return SeitzOperator(vcat(hcat(diagm(0 => ones(T, 3)), t.translation), lastrow(T)))
end # function TranslationOperator

function PointSymmetryOperator(m::LinearMap)
    T = eltype(m)
    return SeitzOperator(vcat(hcat(m.linear, zeros(T, 3)), lastrow(T)))
end # function PointSymmetryOperator

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

lastrow(T::Type{<:Real}) = [zeros(T, 3)... one(T)]

# Base.eltype(t::Translation) = eltype(t.translation)
# Base.eltype(m::LinearMap) = eltype(m.linear)

Base.:*(m::SeitzOperator, c::CrystalCoordinates) = CrystalCoordinates((m.m*[c; 1])[1:3])

Base.:∘(a::SeitzOperator, b::SeitzOperator) = SeitzOperator(a.m * b.m)

    # Base.convert(::Type{Translation}, op::TranslationOperator) = Translation(collect(op.m[1:3, 4]))

    # Base.inv(op::TranslationOperator) = TranslationOperator(inv(convert(Translation, op)))


end # module Symmetry
