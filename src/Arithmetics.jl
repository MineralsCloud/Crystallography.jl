module Arithmetics

using CoordinateTransformations: IdentityTransformation
using EponymTuples: @eponymargs
using LinearAlgebra: I, cross, det, dot, norm
using StaticArrays: SVector, SMatrix, SHermitianCompact

using Crystallography: CellParameters, Lattice, Cell, destruct, cellvolume

import LinearAlgebra
import Crystallography

export RealSpace,
    ReciprocalSpace,
    MetricTensor,
    Miller,
    MillerBravais,
    CrystalFromCartesian,
    CartesianFromCrystal
export directioncosine,
    directionangle,
    distance,
    interplanar_spacing,
    reciprocal,
    @m_str,
    @mb_str

abstract type AbstractSpace end
struct RealSpace <: AbstractSpace end
struct ReciprocalSpace <: AbstractSpace end

struct MetricTensor{T} <: AbstractMatrix{T}
    data::SHermitianCompact{3,T}
end
MetricTensor(m::AbstractMatrix) = MetricTensor(SHermitianCompact{3}(m))
function MetricTensor(v1::AbstractVector, v2::AbstractVector, v3::AbstractVector)
    vecs = (v1, v2, v3)
    return MetricTensor([dot(vecs[i], vecs[j]) for i in 1:3, j in 1:3])
end
function MetricTensor(@eponymargs(a, b, c, α, β, γ))
    g12 = a * b * cos(γ)
    g13 = a * c * cos(β)
    g23 = b * c * cos(α)
    return MetricTensor(SHermitianCompact(SVector(a^2, g12, g13, b^2, g23, c^2)))
end

struct Miller{S<:AbstractSpace} <: AbstractVector{Int}
    data::SVector{3,Int}
    Miller{S}(v) where {S} = new(iszero(v) ? v : v .÷ gcd(v))
end
Miller{S}(i, j, k) where {S} = Miller{S}([i, j, k])

struct MillerBravais{S<:AbstractSpace} <: AbstractVector{Int}
    data::SVector{4,Int}
    function MillerBravais{S}(v) where {S}
        @assert(
            v[3] == -v[1] - v[2],
            "the 3rd index of `MillerBravais` should equal to the negation of the first two!"
        )
        return new(iszero(v) ? v : v .÷ gcd(v))
    end
end
MillerBravais{S}(i, j, k, l) where {S} = MillerBravais{S}([i, j, k, l])

# This is a helper type and should not be exported!
const INDICES = Union{Miller,MillerBravais}

# This is a helper function and should not be exported!
function _indices_str(r::Regex, s::AbstractString, ::Type{T}) where {T<:INDICES}
    m = match(r, strip(s))
    isnothing(m) && error("not a valid expression!")
    brackets = first(m.captures) * last(m.captures)
    x = (parse(Int, x) for x in m.captures[2:(end-1)])
    if brackets ∈ ("()", "{}")
        return T{ReciprocalSpace}(x...)
    elseif brackets ∈ ("[]", "<>")
        return T{RealSpace}(x...)
    else
        error("not a valid expression!")
    end
end # function _indices_str

macro m_str(s)
    r = r"([({[<])\s*([-+]?[0-9]+)[\s,]+([-+]?[0-9]+)[\s,]+([-+]?[0-9]+)[\s,]*([>\]})])"
    _indices_str(r, s, Miller)
end

macro mb_str(s)
    r = r"([({[<])\s*([-+]?[0-9]+)[\s,]+([-+]?[0-9]+)[\s,]+([-+]?[0-9]+)[\s,]+([-+]?[0-9]+)[\s,]*([>\]})])"
    _indices_str(r, s, MillerBravais)
end

struct CartesianFromCrystal
    m::SMatrix{3,3}
end
struct CrystalFromCartesian
    m::SMatrix{3,3}
end
for T in (:CartesianFromCrystal, :CrystalFromCartesian)
    eval(quote
        $T(lattice::Lattice) = $T(convert(Matrix{eltype(lattice)}, lattice))
    end)
end
function CartesianFromCrystal(@eponymargs(a, b, c, α, β, γ))
    v = cellvolume(CellParameters(a, b, c, α, β, γ))
    x, y = sin(γ), cos(γ)
    return CartesianFromCrystal([
        a b * y -c / x^2 * (_F(γ, α, β) + _F(β, α, γ) * y)
        0 b * x -c * _F(β, α, γ) / x
        0 0 v / (a * b * x)
    ])
end # function CartesianFromCrystal
function CrystalFromCartesian(@eponymargs(a, b, c, α, β, γ))  # This is wrong
    v = cellvolume(CellParameters(a, b, c, α, β, γ))
    x = sin(γ)
    return CrystalFromCartesian([
        1 / a -1 / (a * tan(γ)) b * c * _F(γ, α, β) / (v * x)
        0 1 / (b * x) a * c * _F(β, γ, α) / (v * x)
        0 0 a * b * x / v
    ])
end # function CrystalFromCartesian

# This is a helper function and should not be exported!
_F(α, β, γ) = cos(α) * cos(β) - cos(γ)

(x::CartesianFromCrystal)(v::AbstractVector) = x.m * v
(x::CrystalFromCartesian)(v::AbstractVector) = inv(x.m) * v

"""
    cellvolume(g::MetricTensor)

Calculates the cell volume from a `MetricTensor`.
"""
Crystallography.cellvolume(g::MetricTensor) = sqrt(det(g.data))  # `sqrt` is always positive!

function Crystallography.Lattice(@eponymargs(a, b, c, α, β, γ))
    # From https://github.com/LaurentRDC/crystals/blob/dbb3a92/crystals/lattice.py#L321-L354
    v = cellvolume(CellParameters(1, 1, 1, α, β, γ))
    # reciprocal lattice
    a_recip = sin(α) / (a * v)
    csg = (cos(α) * cos(β) - cos(γ)) / (sin(α) * sin(β))
    sg = sqrt(1 - csg^2)
    a1 = [1 / a_recip, -csg / sg / a_recip, cos(β) * a]
    a2 = [0, b * sin(α), b * cos(α)]
    a3 = [0, 0, c]
    return Lattice(a1, a2, a3)
end # function Lattice

function reciprocal(lattice::Lattice)
    volume = cellvolume(lattice)
    a1, a2, a3 = destruct(lattice)
    return 1 / volume * [cross(a2, a3) cross(a3, a1) cross(a1, a2)]
end # function reciprocal

directioncosine(a::AbstractVector, g::MetricTensor, b::AbstractVector) =
    dot(a, g, b) / (norm(a, g) * norm(b, g))

directionangle(a::AbstractVector, g::MetricTensor, b::AbstractVector) =
    acos(directioncosine(a, g, b))

distance(a::AbstractVector, g::MetricTensor, b::AbstractVector) = norm(b - a, g)

interplanar_spacing(a::AbstractVector, g::MetricTensor) = 1 / norm(a, g)

Base.size(::MetricTensor) = (3, 3)
Base.size(::Miller) = (3,)
Base.size(::MillerBravais) = (4,)

Base.getindex(A::MetricTensor, I::Vararg{Int}) = getindex(A.data, I...)
Base.getindex(A::Union{Miller,MillerBravais}, i::Int) = getindex(A.data, i)

Base.inv(g::MetricTensor) = MetricTensor(inv(g.data))
Base.inv(x::CrystalFromCartesian) = CartesianFromCrystal(inv(x.m))
Base.inv(x::CartesianFromCrystal) = CrystalFromCartesian(inv(x.m))
Base.:∘(x::CrystalFromCartesian, y::CartesianFromCrystal) =
    x.m * y.m ≈ I ? IdentityTransformation() : error("undefined!")

Base.convert(::Type{T}, x::T) where {T<:INDICES} = x
Base.convert(::Type{Miller{T}}, mb::MillerBravais{T}) where {T<:RealSpace} =
    Miller{T}(2 * mb[1] + mb[2], 2 * mb[2] + mb[1], mb[4])
Base.convert(::Type{Miller{T}}, mb::MillerBravais{T}) where {T<:ReciprocalSpace} =
    Miller{T}(mb[1], mb[2], mb[4])
Base.convert(::Type{MillerBravais{T}}, m::Miller{T}) where {T<:RealSpace} =
    MillerBravais{T}(2 * m[1] - m[2], 2 * m[2] - m[1], -(m[1] + m[2]), 3 * m[3])
Base.convert(::Type{MillerBravais{T}}, m::Miller{T}) where {T<:ReciprocalSpace} =
    MillerBravais{T}(m[1], m[2], -(m[1] + m[2]), m[3])
function Base.convert(::Type{CellParameters}, g::MetricTensor)
    data = g.data
    a2, b2, c2, ab, ac, bc =
        data[1, 1], data[2, 2], data[3, 3], data[1, 2], data[1, 3], data[2, 3]
    a, b, c = map(sqrt, (a2, b2, c2))
    γ, β, α = acos(ab / (a * b)), acos(ac / (a * c)), acos(bc / (b * c))
    return CellParameters(a, b, c, α, β, γ)
end # function Base.convert
Base.convert(::Type{Lattice}, g::MetricTensor) = Lattice(convert(CellParameters, g))

LinearAlgebra.dot(a::AbstractVector, g::MetricTensor, b::AbstractVector) = a' * g.data * b
LinearAlgebra.norm(a::AbstractVector, g::MetricTensor) = sqrt(dot(a, g, a))

end # module Arithmetics
