"""
# module Prelude



# Examples

```jldoctest
julia>
```
"""
module Prelude

using StaticArrays: FieldVector

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

for operator in (:+, :-)
    eval(quote
        Base.$operator(a::CrystalCoordinates{T}, b::CrystalCoordinates{T}) where {T} = CrystalCoordinates{T}(mapreduce(collect, $operator, (a, b)))
        Base.$operator(a::CrystalCoordinates{T}) where {T} = CrystalCoordinates{T}($operator(collect(a)))
        Base.$operator(a::CartesianCoordinates, b::CartesianCoordinates) = CartesianCoordinates(mapreduce(collect, $operator, (a, b)))
    end)
end
for operator in (:*, :/)
    eval(quote
        Base.$operator(n::Number, a::CrystalCoordinates{T}) where {T} = CrystalCoordinates{T}($operator(collect(a), n))
        Base.$operator(a::CrystalCoordinates, n::Number) = $operator(n, a)
        Base.$operator(n::Number, a::CartesianCoordinates) = CartesianCoordinates($operator(collect(a), n))
        Base.$operator(a::CartesianCoordinates, n::Number) = $operator(n, a)
    end)
end
Base.map(f, a::CrystalCoordinates{T}) where {T} = CrystalCoordinates{T}(map(f, collect(a)))
Base.map(f, a::CartesianCoordinates) = CartesianCoordinates(map(f, collect(a)))

Base.inv(::Type{RealSpace}) = ReciprocalSpace
Base.inv(::Type{ReciprocalSpace}) = RealSpace

end
