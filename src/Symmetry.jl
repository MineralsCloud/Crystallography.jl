module Symmetry

using LinearAlgebra: diagm, I

using CoordinateTransformations: AffineMap, Translation, LinearMap
using LibSymspg: get_symmetry

export SeitzOperator, IdentityOperator, TranslationOperator, PointSymmetryOperator
export Cell, getsymmetry

struct Cell{L<:AbstractVecOrMat,P<:AbstractVecOrMat,N<:AbstractVector,M<:Union{AbstractVector,Nothing}}
    lattice::L
    positions::P
    numbers::N
    magmoms::M
end
Cell(lattice, positions, numbers) = Cell(lattice, positions, numbers, nothing)

function getsymmetry(cell::Cell, symprec::AbstractFloat = 1e-5)
    maps, translations = get_symmetry(cell.lattice, cell.positions, cell.numbers, length(cell.numbers), symprec)
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

struct TranslationOperator{T} <: SeitzOperator{T}
    m::T
    function TranslationOperator{T}(m) where {T}
        @assert(size(m) == (4, 4), "The operator must be of size 4x4!")
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
        return new(m)
    end
end
PointSymmetryOperator(m::T) where {T} = PointSymmetryOperator{T}(m)
function PointSymmetryOperator(m::LinearMap)
    T = eltype(m)
    return PointSymmetryOperator(vcat(hcat(m.linear, zeros(T, 3)), lastrow(T)))
end

lastrow(T::Type{<:Real}) = [zeros(T, 3)... ones(T, 1)]

Base.eltype(t::Translation) = eltype(t.translation)
Base.eltype(m::LinearMap) = eltype(m.linear)

end # module Symmetry
