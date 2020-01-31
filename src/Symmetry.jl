module Symmetry

using LinearAlgebra: diagm, I

using CoordinateTransformations: AffineMap, Translation, LinearMap
using LibSymspg: get_symmetry

using Crystallography: CrystalCoordinates
using Crystallography.Crystals: Cell

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

abstract type SeitzOperator{T<:AbstractMatrix} end

struct IdentityOperator{T} <: SeitzOperator{T}
    m::T
    function IdentityOperator{T}(m) where {T}
        @assert(size(m) == (4, 4), "The operator must be of size 4x4!")
        @assert(m == I, "The matrix is not an identity matrix!")
        return new(m)
    end
end
IdentityOperator(m::T) where {T} = IdentityOperator{T}(m)
IdentityOperator() = IdentityOperator(Matrix{Int}(I, 4, 4))

struct TranslationOperator{T} <: SeitzOperator{T}
    m::T
    function TranslationOperator{T}(m) where {T}
        @assert(size(m) == (4, 4), "The operator must be of size 4x4!")
        @assert(m[1:3, 1:3] == I)
        @assert(iszero(m[4, 1:3]) && isone(m[4, 4]))
        return new(m)
    end
end
TranslationOperator(m::T) where {T} = TranslationOperator{T}(m)
function TranslationOperator(t::Translation)
    T = eltype(t)
    return TranslationOperator(vcat(
        hcat(diagm(0 => ones(T, 3)), t.translation),
        lastrow(T),
    ))
end

struct PointSymmetryOperator{T} <: SeitzOperator{T}
    m::T
    function PointSymmetryOperator{T}(m) where {T}
        @assert(size(m) == (4, 4), "The operator must be of size 4x4!")
        @assert(m[4, :] == [0, 0, 0, 1] && iszero(m[1:3, end]))
        return new(m)
    end
end
PointSymmetryOperator(m::T) where {T} = PointSymmetryOperator{T}(m)
function PointSymmetryOperator(m::LinearMap)
    T = eltype(m)
    return PointSymmetryOperator(vcat(hcat(m.linear, zeros(T, 3)), lastrow(T)))
end

struct GeneralSeitzOperator{T} <: SeitzOperator{T}
    m::T
    function GeneralSeitzOperator{T}(m) where {T}
        @assert(size(m) == (4, 4), "The operator must be of size 4x4!")
        return new(m)
    end
end
GeneralSeitzOperator(m::T) where {T} = GeneralSeitzOperator{T}(m)
function GeneralSeitzOperator(m::AffineMap)
    linear, translation = m.linear, m.translation
    return GeneralSeitzOperator(-translation * linear * translation)
end # function GeneralSeitzOperator

lastrow(T::Type{<:Real}) = [zeros(T, 3)... ones(T, 1)]

Base.eltype(t::Translation) = eltype(t.translation)
Base.eltype(m::LinearMap) = eltype(m.linear)

Base.:*(m::SeitzOperator, c::CrystalCoordinates) = CrystalCoordinates((m.m*[c... 1])[1:3])

end # module Symmetry
