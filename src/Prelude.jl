"""
# module Prelude



# Examples

```jldoctest
julia>
```
"""
module Prelude

export AbstractSpace,
    RealSpace,
    ReciprocalSpace

abstract type AbstractSpace end
struct RealSpace <: AbstractSpace end
struct ReciprocalSpace <: AbstractSpace end

Base.inv(::Type{RealSpace}) = ReciprocalSpace
Base.inv(::Type{ReciprocalSpace}) = RealSpace

end
