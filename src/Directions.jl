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
using StaticArrays: FieldVector

using Crystallography.BravaisLattices

export SpaceType,
    RealSpace,
    ReciprocalSpace,
    MillerIndices,
    MillerBravaisIndices,
    MetricTensor,
    directioncosine,
    directionangle

abstract type SpaceType end
struct RealSpace <: SpaceType end
struct ReciprocalSpace <: SpaceType end

struct MillerIndices{S,T <: Integer} <: FieldVector{3,T}
    i::T
    j::T
    k::T
    function MillerIndices{S,T}(i, j, k) where {S <: SpaceType,T <: Integer}
        x = [i, j, k]
        i, j, k = iszero(x) ? x : x .÷ gcd(x)
        new(i, j, k)
    end
end
MillerIndices{S}(i::T, j::T, k::T) where {S <: SpaceType,T <: Integer} = MillerIndices{S,T}(i, j, k)
MillerIndices{S}(x::AbstractVector{T}) where {S <: SpaceType,T <: Integer} = MillerIndices{S}(x...)
MillerIndices{S}(x::Tuple) where {S} = MillerIndices{S}(collect(x))

struct MillerBravaisIndices{S,T <: Integer} <: FieldVector{4,T}
    i::T
    j::T
    k::T
    l::T
    function MillerBravaisIndices{S,T}(i, j, k, l) where {S <: SpaceType,T <: Integer}
        x = [i, j, k, l]
        i, j, k, l = iszero(x) ? x : x .÷ gcd(x)
        new(i, j, k, l)
    end
end
MillerBravaisIndices{S}(i::T, j::T, k::T, l::T) where {S <: SpaceType,T <: Integer} = MillerBravaisIndices{S,T}(i, j, k, l)
MillerBravaisIndices{S}(x::AbstractVector{T}) where {S <: SpaceType,T <: Integer} = MillerBravaisIndices{S}(x...)
MillerBravaisIndices{S}(x::Tuple) where {S} = MillerBravaisIndices{S}(collect(x))

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
MetricTensor{RealSpace}(::Type{BravaisLattice{RhombohedralCentered,Hexagonal}}, a, α) = MetricTensor{RealSpace}(a, a, a, α, α, α)
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

function Base.getproperty(x::MillerIndices{RealSpace}, name::Symbol)
    getfield(x, Dict(:u => :i, :v => :j, :w => :k)[name])
end
function Base.getproperty(x::MillerIndices{ReciprocalSpace}, name::Symbol)
    getfield(x, Dict(:h => :i, :k => :j, :l => :k)[name])
end
function Base.getproperty(x::MillerBravaisIndices{RealSpace}, name::Symbol)
    getfield(x, Dict(:u => :i, :v => :j, :t => :k, :w => :l)[name])
end
function Base.getproperty(x::MillerBravaisIndices{ReciprocalSpace}, name::Symbol)
    getfield(x, Dict(:h => :i, :k => :j, :i => :k, :l => :l)[name])
end


end
