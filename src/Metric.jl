"""
# module Metric



# Examples

```jldoctest
julia>
```
"""
module Metric

using LinearAlgebra

using Crystallography

export MetricTensor,
    directioncosine,
    directionangle,
    distance

struct MetricTensor{S,T}
    m::T
    function MetricTensor{S,T}(m) where {S <: AbstractSpace,T}
        size(m) == (3, 3) || throw(DimensionMismatch("The metric tensor must be of size 3x3!"))
        new(m)
    end
end
MetricTensor{S}(m::T) where {S,T} = MetricTensor{S,T}(m)
function MetricTensor{S}(v1::AbstractVector, v2::AbstractVector, v3::AbstractVector) where {S <: AbstractSpace}
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
MetricTensor{RealSpace}(::Type{BravaisLattice{RhombohedralCentered,Hexagonal}}, a, α) = MetricTensor{RealSpace}(a, a, a, α, α, α)
MetricTensor{ReciprocalSpace}(args...) = inv(MetricTensor{RealSpace}(args...))

function directioncosine(a::CrystalCoordinates{T}, g::MetricTensor{T}, b::CrystalCoordinates{T}) where {T}
    dot(a, g, b) / (length(a, g) * length(b, g))
end

directionangle(a::CrystalCoordinates{T}, g::MetricTensor{T}, b::CrystalCoordinates{T}) where {T} = acos(directioncosine(a, g, b))

distance(a::CrystalCoordinates{T}, g::MetricTensor{T}, b::CrystalCoordinates{T}) where {T} = length(b - a, g)

Base.inv(::Type{MetricTensor{T}}) where {T} = MetricTensor{inv(T)}
Base.inv(g::MetricTensor) = inv(typeof(g))(inv(g.m))

LinearAlgebra.dot(a::CrystalCoordinates{T}, g::MetricTensor{T}, b::CrystalCoordinates{T}) where {T} = a' * g.m * b

Base.length(a::CrystalCoordinates{T}, g::MetricTensor{T}) where {T} = sqrt(dot(a, g, a))

end
