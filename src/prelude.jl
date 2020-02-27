using LinearAlgebra: cross, det, dot, norm

using StaticArrays: FieldVector, SVector, SMatrix, SHermitianCompact, Size
using SymPy

import LinearAlgebra
import StaticArrays

export AbstractSpace,
    RealSpace,
    ReciprocalSpace,
    CrystalCoordinates,
    CartesianCoordinates,
    CrystalSystem,
    Oblique,
    Rectangular,
    Square,
    Triclinic,
    Monoclinic,
    Orthorhombic,
    Tetragonal,
    Cubic,
    Trigonal,
    Hexagonal,
    Centering,
    BaseCentered,
    Primitive,
    BodyCentered,
    FaceCentered,
    RhombohedralCentered,
    BaseCentered,
    BravaisLattice,
    MetricTensor,
    MillerIndices,
    MillerBravaisIndices,
    Cell,
    LatticeConstants,
    AxisAngles,
    CellParameters,
    Lattice

export pearsonsymbol,
    arithmeticclass,
    centeringof,
    crystalsystem,
    dimensionof,
    directioncosine,
    directionangle,
    distance,
    interplanar_spacing,
    cellvolume,
    reciprocalof

abstract type AbstractSpace end
struct RealSpace <: AbstractSpace end
struct ReciprocalSpace <: AbstractSpace end

abstract type CrystalSystem{N} end
struct Oblique <: CrystalSystem{2} end
struct Rectangular <: CrystalSystem{2} end
struct Square <: CrystalSystem{2} end
struct Triclinic <: CrystalSystem{3} end
struct Monoclinic <: CrystalSystem{3} end
struct Orthorhombic <: CrystalSystem{3} end
struct Tetragonal <: CrystalSystem{3} end
struct Cubic <: CrystalSystem{3} end
struct Trigonal <: CrystalSystem{3} end
struct Hexagonal{N} <: CrystalSystem{N} end  # Could both be 2D or 3D
Hexagonal(N::Int = 3) =
    N ∈ (2, 3) ? Hexagonal{N}() : throw(ArgumentError("hexagonal must be 2D or 3D!"))

abstract type Centering end
struct Primitive <: Centering end
struct BodyCentered <: Centering end
struct FaceCentered <: Centering end
struct RhombohedralCentered <: Centering end
struct BaseCentered{T} <: Centering end
BaseCentered(T::Symbol) = T ∈ (:A, :B, :C) ? BaseCentered{T}() :
    throw(ArgumentError("centering must be either :A, :B, or :C!"))

struct BravaisLattice{A<:CrystalSystem,B<:Centering} end
BravaisLattice(::A, ::B) where {A,B} = BravaisLattice{A,B}()

pearsonsymbol(::Triclinic) = "a"
pearsonsymbol(::Monoclinic) = "m"
pearsonsymbol(::Orthorhombic) = "o"
pearsonsymbol(::Tetragonal) = "t"
pearsonsymbol(::Cubic) = "c"
pearsonsymbol(::Hexagonal) = "h"
pearsonsymbol(::Trigonal) = "h"
pearsonsymbol(::Primitive) = "P"
pearsonsymbol(::BaseCentered{T}) where {T} = string(T)
pearsonsymbol(::BodyCentered) = "I"
pearsonsymbol(::FaceCentered) = "F"
pearsonsymbol(::RhombohedralCentered) = "R"
pearsonsymbol(::BravaisLattice{A,B}) where {A,B} = pearsonsymbol(A()) * pearsonsymbol(B())

arithmeticclass(::Oblique) = "2"
arithmeticclass(::Rectangular) = "2mm"
arithmeticclass(::Square) = "4mm"
arithmeticclass(::Hexagonal{2}) = "6mm"
arithmeticclass(::Triclinic) = "-1"
arithmeticclass(::Monoclinic) = "2/m"
arithmeticclass(::Orthorhombic) = "mmm"
arithmeticclass(::Tetragonal) = "4/mmm"
arithmeticclass(::Hexagonal{3}) = "6/mmm"
arithmeticclass(::Cubic) = "m-3m"
arithmeticclass(::Primitive) = "P"
arithmeticclass(::BaseCentered) = "S"
arithmeticclass(::BodyCentered) = "I"
arithmeticclass(::FaceCentered) = "F"
arithmeticclass(::RhombohedralCentered) = "R"
arithmeticclass(::BravaisLattice{A,B}) where {A,B} =
    arithmeticclass(A()) * arithmeticclass(B())
arithmeticclass(::BravaisLattice{Hexagonal{3},RhombohedralCentered}) = "-3mR"
arithmeticclass(::BravaisLattice{Oblique}) = "2p"
arithmeticclass(::BravaisLattice{Rectangular,Primitive}) = "2mmp"
arithmeticclass(::BravaisLattice{Rectangular,<:BaseCentered}) = "2mmc"
arithmeticclass(::BravaisLattice{Square}) = "4mmp"
arithmeticclass(::BravaisLattice{Hexagonal{2}}) = "6mmh"

centeringof(::BravaisLattice{C,T}) where {C,T} = T()

dimensionof(::CrystalSystem{N}) where {N} = N
dimensionof(::BravaisLattice{C}) where {C} = dimensionof(C())

abstract type AbstractCoordinates{T} <: FieldVector{3,T} end
struct CrystalCoordinates{T} <: AbstractCoordinates{T}
    x::T
    y::T
    z::T
end
struct CartesianCoordinates{T} <: AbstractCoordinates{T}
    x::T
    y::T
    z::T
end

struct LatticeConstants{T} <: FieldVector{3,T}
    a::T
    b::T
    c::T
    function LatticeConstants{T}(a, b, c) where {T}
        @assert all((a, b, c) .> 0)
        return new(a, b, c)
    end
end
LatticeConstants(a::T, b::T, c::T) where {T} = LatticeConstants{T}(a, b, c)
LatticeConstants(
    ::BravaisLattice{T},
    a,
    b,
    c,
) where {T<:Union{Triclinic,Monoclinic,Orthorhombic}} = LatticeConstants(a, b, c)
LatticeConstants(::BravaisLattice{Tetragonal}, a, b, c) = LatticeConstants(a, a, c)
LatticeConstants(::BravaisLattice{Cubic}, a, b, c) = LatticeConstants(a, a, a)
LatticeConstants(::BravaisLattice{Hexagonal{3},Primitive}, a, b, c) =
    LatticeConstants(a, a, c)
LatticeConstants(::BravaisLattice{Hexagonal{3},RhombohedralCentered}, a, b, c) =
    LatticeConstants(a, a, a)

struct AxisAngles{T} <: FieldVector{3,T}
    α::T
    β::T
    γ::T
end
AxisAngles(α, β, γ, unit::Symbol = :deg) = AxisAngles(Degree(α), Degree(β), Degree(γ), Val(unit))
AxisAngles(::BravaisLattice{Triclinic}, α, β, γ) = AxisAngles(α, β, γ)
AxisAngles(::BravaisLattice{Monoclinic,Primitive}, α, β, γ; unit::Symbol = :deg, view::Int = 1) =
    view == 1 ? AxisAngles(90, 90, γ, unit) : AxisAngles(90, β, 90, unit)
AxisAngles(::BravaisLattice{Monoclinic,BaseCentered{:C}}, α, β, γ; args...) =
    AxisAngles(90, 90, γ; args...)
AxisAngles(::BravaisLattice{Monoclinic,BaseCentered{:B}}, α, β, γ; args...) =
    AxisAngles(90, β, 90; args...)
AxisAngles(
    ::BravaisLattice{T},
    α,
    β,
    γ;
    args...,
) where {T<:Union{Orthorhombic,Tetragonal,Cubic}} = AxisAngles(Degree(90), Degree(90), Degree(90); args...)
AxisAngles(::BravaisLattice{Hexagonal{3},Primitive}, α, β, γ; args...) =
    AxisAngles(90, 90, 120; args...)
AxisAngles(::BravaisLattice{Hexagonal{3},RhombohedralCentered}, α, β, γ; args...) =
    AxisAngles(α, α, α; args...)
AxisAngles(α::Degree, β::Degree, γ::Degree, ::Val{:deg}) = AxisAngles(α, β, γ)
AxisAngles(α::Degree, β::Degree, γ::Degree, ::Val{:rad}) = AxisAngles(Radian(α), Radian(β), Radian(γ))
AxisAngles(α::Degree, β::Degree, γ::Degree, ::Val{:cos}) = AxisAngles(Cosine(α), Cosine(β), Cosine(γ))
AxisAngles(α::Radian, β::Radian, γ::Radian) = AxisAngles(rad2deg(α), rad2deg(β), rad2deg(γ))
AxisAngles(α::Cosine, β::Cosine, γ::Cosine) = AxisAngles(_cos2deg(α), _cos2deg(β), _cos2deg(γ))

struct Degree{N} end
struct Radian{N} end
struct Cosine{N} end
_cos2deg(x) = rad2deg(acos(x))

struct CellParameters{S,T}
    x::LatticeConstants{S}
    y::AxisAngles{T}
end
CellParameters(a::S, b::S, c::S, α::T, β::T, γ::T; args...) where {S,T} =
    CellParameters{S,T}(LatticeConstants(a, b, c), AxisAngles(α, β, γ; args...))
CellParameters(x::BravaisLattice) = args -> CellParameters(x, args...)

struct Lattice{T} <: AbstractMatrix{T}
    m::SMatrix{3,3,T}
end
Lattice(m::AbstractMatrix) = Lattice(SMatrix{3,3}(m))
Lattice(v1::AbstractVector, v2::AbstractVector, v3::AbstractVector) =
    vcat(transpose.((v1, v2, v3))...)
function Lattice(p::CellParameters)
    # From https://github.com/LaurentRDC/crystals/blob/dbb3a92/crystals/lattice.py#L321-L354
    a, b, c, α, β, γ = p
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

struct Cell{
    L<:AbstractVecOrMat,
    P<:AbstractVecOrMat,
    N<:AbstractVector,
    M<:Union{AbstractVector,Nothing},
}
    lattice::L
    positions::P
    numbers::N
    magmoms::M
end
Cell(lattice, positions, numbers) = Cell(lattice, positions, numbers, nothing)

struct MetricTensor{T} <: AbstractMatrix{T}
    m::SHermitianCompact{3,T}
end
MetricTensor(m::AbstractMatrix) = MetricTensor(SHermitianCompact{3}(m))
function MetricTensor(v1::AbstractVector, v2::AbstractVector, v3::AbstractVector)
    vecs = (v1, v2, v3)
    return MetricTensor([dot(vecs[i], vecs[j]) for i in 1:3, j in 1:3])
end
function MetricTensor(a, b, c, α, β, γ)
    g12 = a * b * cos(γ)
    g13 = a * c * cos(β)
    g23 = b * c * cos(α)
    return MetricTensor(SHermitianCompact(SVector(a^2, g12, g13, b^2, g23, c^2)))
end
MetricTensor(p::CellParameters) = MetricTensor(p...)

struct MillerIndices{S<:AbstractSpace} <: AbstractVector{Int}
    v::SVector{3,Int}
    function MillerIndices{S}(x) where {S}
        return new(iszero(x) ? x : x .÷ gcd(x))
    end
end
MillerIndices{S}(i, j, k) where {S} = MillerIndices{S}([i, j, k])

struct MillerBravaisIndices{S<:AbstractSpace} <: AbstractVector{Int}
    v::SVector{4,Int}
    function MillerBravaisIndices{S}(x) where {S}
        return new(iszero(x) ? x : x .÷ gcd(x))
    end
end
MillerBravaisIndices{S}(i, j, k, l) where {S} = MillerBravaisIndices{S}([i, j, k, l])

# This is a helper type and should not be exported!
const INDICES = Union{MillerIndices,MillerBravaisIndices}

crystalsystem(::BravaisLattice{C}) where {C} = C()
function crystalsystem(p::CellParameters)
    a, b, c, α, β, γ = p
    if a == b == c
        if α == β == γ
            α == 90 ? Cubic() : Trigonal()
        else
            α == β == 90 && γ == 120 ? Hexagonal() : Triclinic()
        end
    else
        if α == β == γ == 90
            a == b || a == c || b == c ? Tetragonal() : Orthorhombic()
        else
            α == β == 90 || β == γ == 90 || α == γ == 90 ? Monoclinic() : Triclinic()
        end
    end
end # function whatsystem
function crystalsystem(lattice::AbstractMatrix)
    v1, v2, v3 = _splitlattice(lattice)
    a, b, c = norm(v1), norm(v2), norm(v3)
    γ = acos(dot(v1, v2) / a / b)
    β = acos(dot(v2, v3) / b / c)
    α = acos(dot(v1, v3) / a / c)
    return crystalsystem(CellParameters(a, b, c, α, β, γ))
end # function crystalsystem

# This is a helper function and should not be exported.
_splitlattice(m::AbstractMatrix) = collect(Iterators.partition(m', 3))

function _checkpositive(v)  # This is a helper function and should not be exported.
    v <= zero(v) && @warn "The volume of the cell is not positive! Check your input!"
    return v
end # function _checkpositive

"""
    cellvolume(param::CellParameters)

Calculates the cell volume from a set of `CellParameters`.
"""
function cellvolume(param::CellParameters)
    a, b, c, α, β, γ = param
    return a * b * c * sqrt(sin(α)^2 - cos(β)^2 - cos(γ)^2 + 2 * cos(α) * cos(β) * cos(γ))
end # function cellvolume
"""
    cellvolume(m::AbstractMatrix)

Calculates the cell volume from a general 3×3 matrix.
"""
function cellvolume(m::AbstractMatrix)
    @assert(size(m) == (3, 3), "The matrix must be of size 3×3!")
    return _checkpositive(det(m))
end # function cellvolume
"""
    cellvolume(g::MetricTensor)

Calculates the cell volume from a `MetricTensor`.
"""
cellvolume(g::MetricTensor) = sqrt(det(g.m))  # `sqrt` is always positive!

function reciprocalof(mat::AbstractMatrix, twopi::Bool = false)
    @assert size(mat) == (3, 3)
    volume = abs(det(mat))
    a1, a2, a3 = mat[1, :], mat[2, :], mat[3, :]
    factor = twopi ? 2 * SymPy.PI : 1
    return factor / volume * [cross(a2, a3) cross(a3, a1) cross(a1, a2)]
end # function reciprocalof

directioncosine(a::CrystalCoordinates, g::MetricTensor, b::CrystalCoordinates) =
    dot(a, g, b) / (norm(a, g) * norm(b, g))

directionangle(a::CrystalCoordinates, g::MetricTensor, b::CrystalCoordinates) =
    acos(directioncosine(a, g, b))

distance(a::CrystalCoordinates, g::MetricTensor, b::CrystalCoordinates) = norm(b - a, g)

interplanar_spacing(a::CrystalCoordinates, g::MetricTensor) = 1 / norm(a, g)

Base.size(::Union{MetricTensor,Lattice}) = (3, 3)
Base.size(::Union{MillerIndices}) = (3,)
Base.size(::Union{MillerBravaisIndices}) = (4,)

Base.getindex(A::Union{MetricTensor,Lattice}, I::Vararg{Int}) = getindex(A.m, I...)
Base.getindex(A::Union{MillerIndices,MillerBravaisIndices}, i::Int) = getindex(A.v, i)

Base.inv(g::MetricTensor) = MetricTensor(inv(SymPy.N(g.m)))

Base.convert(::Type{T}, x::T) where {T<:INDICES} = x
Base.convert(::Type{MillerIndices{T}}, mb::MillerBravaisIndices{T}) where {T<:RealSpace} =
    MillerIndices{T}(2 * mb[1] + mb[2], 2 * mb[2] + mb[1], mb[4])
Base.convert(
    ::Type{MillerIndices{T}},
    mb::MillerBravaisIndices{T},
) where {T<:ReciprocalSpace} = MillerIndices{T}(mb[1], mb[2], mb[4])
Base.convert(::Type{MillerBravaisIndices{T}}, m::MillerIndices{T}) where {T<:RealSpace} =
    MillerBravaisIndices{T}(2 * m[1] - m[2], 2 * m[2] - m[1], -(m[1] + m[2]), 3 * m[3])
Base.convert(
    ::Type{MillerBravaisIndices{T}},
    m::MillerIndices{T},
) where {T<:ReciprocalSpace} = MillerBravaisIndices{T}(m[1], m[2], -(m[1] + m[2]), m[3])

LinearAlgebra.dot(a::CrystalCoordinates, g::MetricTensor, b::CrystalCoordinates) =
    a' * g.m * b
LinearAlgebra.norm(a::CrystalCoordinates, g::MetricTensor) = sqrt(dot(a, g, a))

StaticArrays.similar_type(
    ::Type{<:CrystalCoordinates},  # Do not delete the `<:`!
    ::Type{T},
    size::Size{(3,)},
) where {T} = CrystalCoordinates{T}
StaticArrays.similar_type(
    ::Type{<:CartesianCoordinates},  # Do not delete the `<:`!
    ::Type{T},
    size::Size{(3,)},
) where {T} = CartesianCoordinates{T}
StaticArrays.similar_type(::Type{<:CellParameters}, ::Type{T}, size::Size{(6,)}) where {T} =
    CellParameters{T}
