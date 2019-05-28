"""
# module Directions



# Examples

```jldoctest
julia>
```
"""
module Directions

using LinearAlgebra: dot

export SpaceType,
    RealSpace,
    ReciprocalSpace,
    CrystalDirection,
    MetricTensor

abstract type SpaceType end
struct RealSpace <: SpaceType end
struct ReciprocalSpace <: SpaceType end

struct CrystalDirection{S,T}
    v::T
    function CrystalDirection{S,T}(v) where {S <: SpaceType,T}
        length(v) == 3 || throw(DimensionMismatch("The crystal direction must be of length 3!"))
        eltype(v) <: Integer || error("The crystal direction must be a vector of integers!")
        iszero(collect(v)) && error("The crystal direction must be a non-zero vector!")
        new(v ./ gcd(v...))
    end
end
CrystalDirection{S}(v::T) where {S,T} = CrystalDirection{S,T}(v)

struct MetricTensor{S,T}
    m::T
    function MetricTensor{S,T}(m) where {S <: SpaceType,T}
        size(m) == (3, 3) || throw(DimensionMismatch("The metric tensor must be of size 3x3!"))
        new(m)
    end
end
MetricTensor{S}(m::T) where {S,T} = MetricTensor{S,T}(m)
function MetricTensor{S}(v1::AbstractVector, v2::AbstractVector, v3::AbstractVector) where {S <: SpaceType}
    vecs = (v1, v2, v3)
    MetricTensor{S}(map(x->dot(x...), Iterators.product(vecs, vecs)))
end

end
