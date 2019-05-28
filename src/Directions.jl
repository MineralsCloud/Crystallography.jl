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
    MetricTensor,
    inv

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
function MetricTensor{RealSpace}(a, b, c, α, β, γ)
    g12 = a * b * cos(γ)
    g13 = a * c * cos(β)
    g23 = b * c * cos(α)
    MetricTensor{RealSpace}([
        a^2 g12 g13
        g12 b^2 g23
        g13 g23 c^2
    ])
end
MetricTensor{ReciprocalSpace}(a, b, c, α, β, γ) = inv(MetricTensor{RealSpace}(a, b, c, α, β, γ))

Base.inv(::Type{RealSpace}) = ReciprocalSpace
Base.inv(::Type{ReciprocalSpace}) = RealSpace
Base.inv(T::Type{<: MetricTensor}) = MetricTensor{inv(first(T.parameters))}
Base.inv(g::MetricTensor) = inv(typeof(g))(inv(g.m))

end
