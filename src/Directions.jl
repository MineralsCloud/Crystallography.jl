"""
# module Directions



# Examples

```jldoctest
julia>
```
"""
module Directions

using StaticArrays: FieldVector

using Crystallography: AbstractSpace, RealSpace, ReciprocalSpace

import Crystallography

export MillerIndices, MillerBravaisIndices

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

Crystallography.CrystalCoordinates(m::MillerIndices) = CrystalCoordinates(m.i, m.j, m.k)
Crystallography.crystalCoordinates(mb::MillerBravaisIndices{T}) where {T} =
    CrystalCoordinates(convert(MillerIndices{T}, mb))

end
