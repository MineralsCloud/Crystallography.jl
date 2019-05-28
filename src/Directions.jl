"""
# module Directions



# Examples

```jldoctest
julia>
```
"""
module Directions

using LinearAlgebra

using CoordinateTransformations: Translation

using Crystallography.BravaisLattices

export SpaceType,
    RealSpace,
    ReciprocalSpace,
    MillerIndices,
    MetricTensor,
    directioncosine,
    directionangle

abstract type SpaceType end
struct RealSpace <: SpaceType end
struct ReciprocalSpace <: SpaceType end

struct MillerIndices{S,T}
    v::T
    function MillerIndices{S,T}(v) where {S <: SpaceType,T}
        length(v) == 3 || throw(DimensionMismatch("The Miller indices must be of length 3!"))
        eltype(v) <: Integer || error("The Miller indices must be a vector of integers!")
        iszero(collect(v)) && error("The Miller indices must be a non-zero vector!")
        new(v ./ gcd(v...))
    end
end
MillerIndices{S}(v::T) where {S,T} = MillerIndices{S,T}(v)

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
function MetricTensor{RealSpace}(a::Real, b::Real, c::Real, α::Real, β::Real, γ::Real)
    g12 = a * b * cos(γ)
    g13 = a * c * cos(β)
    g23 = b * c * cos(α)
    MetricTensor{RealSpace}([
        a^2 g12 g13
        g12 b^2 g23
        g13 g23 c^2
    ])
end
MetricTensor{RealSpace}(::Type{Triclinic}, args...) = MetricTensor{RealSpace}(args...)
MetricTensor{RealSpace}(::Type{Monoclinic}, a, b, c, γ) = MetricTensor{RealSpace}(a, b, c, π / 2, π / 2, γ)
MetricTensor{RealSpace}(::Type{Orthorhombic}, a, b, c) = MetricTensor{RealSpace}(a, b, c, π / 2, π / 2, π / 2)
MetricTensor{RealSpace}(::Type{Tetragonal}, a, c) = MetricTensor{RealSpace}(a, a, c, π / 2, π / 2, π / 2)
MetricTensor{RealSpace}(::Type{Tetragonal}, a) = MetricTensor{RealSpace}(a, a, a, π / 2, π / 2, π / 2)
MetricTensor{RealSpace}(::Type{Hexagonal}, a, c) = MetricTensor{RealSpace}(a, a, c, π / 2, π / 2, 2π / 3)
MetricTensor{RealSpace}(::Type{Trigonal}, a, c) = MetricTensor{RealSpace}(Hexagonal, a, c)
MetricTensor{RealSpace}(::Type{BravaisLattice{RhombohedralCentered, Hexagonal}}, a, α) = MetricTensor{RealSpace}(a, a, a, α, α, α)
MetricTensor{ReciprocalSpace}(args...) = inv(MetricTensor{RealSpace}(args...))

function directioncosine(a::Translation, g::MetricTensor, b::Translation)
    dot(a, g, b) / (length(a, g) * length(b, g))
end

directionangle(a::Translation, g::MetricTensor, b::Translation) = acos(directioncosine(a, g, b))

LinearAlgebra.dot(a::Translation, g::MetricTensor, b::Translation) = a.translation' * g.m * b.translation

Base.length(a::Translation, g::MetricTensor) = sqrt(dot(a, g, a))

Base.inv(::Type{RealSpace}) = ReciprocalSpace
Base.inv(::Type{ReciprocalSpace}) = RealSpace
Base.inv(T::Type{<: MetricTensor}) = MetricTensor{inv(first(T.parameters))}
Base.inv(g::MetricTensor) = inv(typeof(g))(inv(g.m))

end
