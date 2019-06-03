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
    CrystalCoordinates

abstract type AbstractSpace end
struct RealSpace <: AbstractSpace end
struct ReciprocalSpace <: AbstractSpace end

struct CrystalCoordinates{S <: AbstractSpace,T} <: FieldVector{3,T}
    x::T
    y::T
    z::T
end
CrystalCoordinates{S}(x::T, y::T, z::T) where {S,T} = CrystalCoordinates{S,T}(x, y, z)

for operator in (:+, :-)
    eval(quote
        Base.$operator(a::CrystalCoordinates{T}, b::CrystalCoordinates{T}) where {T} = CrystalCoordinates{T}(collect(a) + collect(b))
    end)
end
for operator in (:*, :/, :÷)
    eval(quote
        Base.$operator(n::Number, a::CrystalCoordinates{T}) where {T} = CrystalCoordinates{T}($operator(collect(a), n))
        Base.$operator(a::T, n::Number) where {T <: CrystalCoordinates} = $operator(n, a)
    end)
end

Base.inv(::Type{RealSpace}) = ReciprocalSpace
Base.inv(::Type{ReciprocalSpace}) = RealSpace

end
