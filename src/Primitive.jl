"""
# module Primitive



# Examples

```jldoctest
julia>
```
"""
module Primitive

export SpaceType,
    RealSpace,
    ReciprocalSpace

abstract type SpaceType end
struct RealSpace <: SpaceType end
struct ReciprocalSpace <: SpaceType end

Base.inv(::Type{RealSpace}) = ReciprocalSpace
Base.inv(::Type{ReciprocalSpace}) = RealSpace

end
