"""
# module Prelude



# Examples

```jldoctest
julia>
```
"""
module Prelude

using StaticArrays

export AbstractSpace,
    RealSpace,
    ReciprocalSpace,
    CrystalCoordinates,
    CartesianCoordinates

abstract type AbstractSpace end
struct RealSpace <: AbstractSpace end
struct ReciprocalSpace <: AbstractSpace end

abstract type AbstractCoordinates{T} <: FieldVector{3,T} end

struct CrystalCoordinates{S <: AbstractSpace,T} <: AbstractCoordinates{T}
    x::T
    y::T
    z::T
end
CrystalCoordinates{S}(x::T, y::T, z::T) where {S,T} = CrystalCoordinates{S,T}(x, y, z)

struct CartesianCoordinates{T} <: AbstractCoordinates{T}
    x::T
    y::T
    z::T
end

StaticArrays.similar_type(::Type{<: CrystalCoordinates{S}}, ::Type{T}, size::Size{(3,)}) where {S, T} = CrystalCoordinates{S, T}
StaticArrays.similar_type(::Type{<: CartesianCoordinates}, ::Type{T}, size::Size{(3,)}) where {T} = CartesianCoordinates{T}

Base.inv(::Type{RealSpace}) = ReciprocalSpace
Base.inv(::Type{ReciprocalSpace}) = RealSpace

end
