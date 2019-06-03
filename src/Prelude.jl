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

Base.inv(::Type{RealSpace}) = ReciprocalSpace
Base.inv(::Type{ReciprocalSpace}) = RealSpace

end
