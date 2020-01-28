"""
# module Directions



# Examples

```jldoctest
julia>
```
"""
module Directions

using LinearAlgebra: cross, det, dot

using StaticArrays: FieldVector

using Crystallography:
    AbstractSpace,
    RealSpace,
    ReciprocalSpace,
    CrystalCoordinates,
    BravaisLattice,
    Triclinic,
    Monoclinic,
    Orthorhombic,
    Tetragonal,
    Cubic,
    Hexagonal,
    Trigonal,
    RhombohedralCentered

import LinearAlgebra
import Crystallography

export MetricTensor, MillerIndices, MillerBravaisIndices
export directioncosine, directionangle, distance, interplanar_spacing

struct MetricTensor{T}
    m::T
    function MetricTensor{T}(m) where {T}
        @assert(size(m) == (3, 3), "The metric tensor must be of size 3×3!")
        return new(m)
    end
end
MetricTensor(m::T) where {T} = MetricTensor{T}(m)
function MetricTensor(v1::AbstractVector, v2::AbstractVector, v3::AbstractVector)
    vecs = (v1, v2, v3)
    return MetricTensor(map(x -> dot(x...), Iterators.product(vecs, vecs)))
end
function MetricTensor(a::Real, b::Real, c::Real, α::Real, β::Real, γ::Real)
    g12 = a * b * cos(γ)
    g13 = a * c * cos(β)
    g23 = b * c * cos(α)
    return MetricTensor([
        a^2 g12 g13
        g12 b^2 g23
        g13 g23 c^2
    ])
end
MetricTensor(::BravaisLattice, args...) = MetricTensor(args...)  # Triclinic
MetricTensor(::BravaisLattice{C,Monoclinic}, a, b, c, γ) where {C} =
    MetricTensor(a, b, c, π / 2, π / 2, γ)
MetricTensor(::BravaisLattice{C,Orthorhombic}, a, b, c) where {C} =
    MetricTensor(a, b, c, π / 2, π / 2, π / 2)
MetricTensor(::BravaisLattice{C,Tetragonal}, a, c) where {C} =
    MetricTensor(a, a, c, π / 2, π / 2, π / 2)
MetricTensor(::BravaisLattice{C,Cubic}, a) where {C} =
    MetricTensor(a, a, a, π / 2, π / 2, π / 2)
MetricTensor(::BravaisLattice{C,Hexagonal}, a, c) where {C} =
    MetricTensor(a, a, c, π / 2, π / 2, 2π / 3)
MetricTensor(::BravaisLattice{C,Trigonal}, a, c) where {C} =
    MetricTensor(BravaisLattice(Primitive(), Hexagonal()), a, c)
MetricTensor(::BravaisLattice{RhombohedralCentered,Hexagonal}, a, α) =
    MetricTensor(a, a, a, α, α, α)

function directioncosine(a::CrystalCoordinates, g::MetricTensor, b::CrystalCoordinates)
    return dot(a, g, b) / (length(a, g) * length(b, g))
end

directionangle(a::CrystalCoordinates, g::MetricTensor, b::CrystalCoordinates) =
    acos(directioncosine(a, g, b))

distance(a::CrystalCoordinates, g::MetricTensor, b::CrystalCoordinates) where {T} =
    length(b - a, g)

interplanar_spacing(a::CrystalCoordinates, g::MetricTensor) = 1 / length(a, g)

function reciprocalof(mat::AbstractMatrix)
    @assert size(mat) == (3, 3)
    volume = abs(det(mat))
    a1, a2, a3 = mat[1, :], mat[2, :], mat[3, :]
    return 2π / volume * [cross(a2, a3) cross(a3, a1) cross(a1, a2)]
end # function reciprocalof

struct MillerIndices{S<:AbstractSpace,T<:Integer} <: FieldVector{3,T}
    i::T
    j::T
    k::T
    function MillerIndices{S,T}(i, j, k) where {S,T}
        x = [i, j, k]
        i, j, k = iszero(x) ? x : x .÷ gcd(x)
        return new(i, j, k)
    end
end
MillerIndices{S}(i::T, j::T, k::T) where {S,T} = MillerIndices{S,T}(i, j, k)
MillerIndices{S}(x::AbstractVector) where {S} = MillerIndices{S}(x...)
MillerIndices{S}(x::Tuple) where {S} = MillerIndices{S}(collect(x))

struct MillerBravaisIndices{S<:AbstractSpace,T<:Integer} <: FieldVector{4,T}
    i::T
    j::T
    k::T
    l::T
    function MillerBravaisIndices{S,T}(i, j, k, l) where {S,T}
        x = [i, j, k, l]
        i, j, k, l = iszero(x) ? x : x .÷ gcd(x)
        return new(i, j, k, l)
    end
end
MillerBravaisIndices{S}(i::T, j::T, k::T, l::T) where {S,T} =
    MillerBravaisIndices{S,T}(i, j, k, l)
MillerBravaisIndices{S}(x::AbstractVector) where {S} = MillerBravaisIndices{S}(x...)
MillerBravaisIndices{S}(x::Tuple) where {S} = MillerBravaisIndices{S}(collect(x))

Base.convert(::Type{<:MillerIndices}, mb::MillerBravaisIndices) =
    error("Unsupported operation!")
Base.convert(::Type{<:MillerBravaisIndices}, m::MillerIndices) =
    error("Unsupported operation!")
Base.convert(T::Type{<:MillerIndices}, m::MillerIndices) =
    isa(m, T) ? m : error("Unsupported operation!")
Base.convert(T::Type{<:MillerBravaisIndices}, mb::MillerBravaisIndices) =
    isa(mb, T) ? mb : error("Unsupported operation!")
Base.convert(::Type{MillerIndices{T}}, mb::MillerBravaisIndices{T}) where {T<:RealSpace} =
    MillerIndices{T}(2 * mb[1] + mb[2], 2 * mb[2] + mb[1], mb[4])
Base.convert(
    ::Type{MillerIndices{T}},
    mb::MillerBravaisIndices{T},
) where {T<:ReciprocalSpace} = MillerIndices{T}(mb[1], mb[2], mb[4])
Base.convert(::Type{MillerBravaisIndices{T}}, m::MillerIndices{T}) where {T<:RealSpace} =
    MillerBravaisIndices{T}(2 * m[1] - m[2], 2 * m[2] - m[1], -(m[1] + m[2]), 3 * m[3])
Base.convert(
    ::Type{MillerBravaisIndices{T}},
    m::MillerIndices{T},
) where {T<:ReciprocalSpace} = MillerBravaisIndices{T}(m[1], m[2], -(m[1] + m[2]), m[3])

LinearAlgebra.dot(a::CrystalCoordinates, g::MetricTensor, b::CrystalCoordinates) =
    a' * g.m * b

Base.length(a::CrystalCoordinates, g::MetricTensor) = sqrt(dot(a, g, a))

Crystallography.CrystalCoordinates(m::MillerIndices) = CrystalCoordinates(m.i, m.j, m.k)
Crystallography.CrystalCoordinates(mb::MillerBravaisIndices{T}) where {T} =
    CrystalCoordinates(convert(MillerIndices{T}, mb))

end
