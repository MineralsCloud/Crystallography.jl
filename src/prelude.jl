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
    makelattice,
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

crystalsystem(::BravaisLattice{C}) where {C} = C()

dimensionof(c::CrystalSystem) = first(supertype(typeof(c)).parameters)
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

struct CellParameters{T} <: FieldVector{6,T}
    a::T
    b::T
    c::T
    α::T
    β::T
    γ::T
    function CellParameters{T}(a, b, c, α, β, γ) where {T}
        @assert all(x > zero(T) for x in (a, b, c))
        return new(a, b, c, α, β, γ)
    end
end
CellParameters(a::T, b::T, c::T, α::T, β::T, γ::T) where {T} =
    CellParameters{T}(a, b, c, α, β, γ)
CellParameters(a, b, c, α, β, γ, angletype::Symbol = :deg) = CellParameters(a, b, c, α, β, γ, Val(angletype))
CellParameters(a, b, c, α, β, γ, ::Val{:deg}) = CellParameters(a, b, c, α, β, γ)
CellParameters(a, b, c, α, β, γ, ::Val{:rad}) = CellParameters(a, b, c, rad2deg(α), rad2deg(β), rad2deg(γ))
function CellParameters(a, b, c, α, β, γ, ::Val{:cos})
    v = (a, b, c, acos(α), acos(β), acos(γ))
    return CellParameters{Base.promote_typeof(v...)}(v...)
end
CellParameters(bravais::BravaisLattice) = args -> CellParameters(bravais, args...)
CellParameters(::BravaisLattice{Triclinic}, a, b, c, α, β, γ, args...) =
    CellParameters(a, b, c, α, β, γ)  # Triclinic
CellParameters(::BravaisLattice{Monoclinic,Primitive}, a, b, c, α, β, γ, args...) =
    CellParameters(a, b, c, SymPy.PI / 2, SymPy.PI / 2, γ)  # `α`, `β` are ignored.
CellParameters(::BravaisLattice{Monoclinic,Primitive}, a, b, c, α, β, γ, args...) =
    CellParameters(a, b, c, SymPy.PI / 2, β, SymPy.PI / 2)  # `α`, `γ` are ignored.
CellParameters(::BravaisLattice{Monoclinic,BaseCentered{:C}}, a, b, c, α, β, γ, args...) =
    CellParameters(a, b, c, SymPy.PI / 2, SymPy.PI / 2, γ)  # `α`, `β` are ignored.
CellParameters(::BravaisLattice{Monoclinic,BaseCentered{:B}}, a, b, c, α, β, γ, args...) =
    CellParameters(a, b, c, SymPy.PI / 2, β, SymPy.PI / 2)  # `α`, `γ` are ignored.
CellParameters(::BravaisLattice{Orthorhombic}, a, b, c, args...) =
    CellParameters(a, b, c, SymPy.PI / 2, SymPy.PI / 2, SymPy.PI / 2)  # `α`, `β`, `γ` are ignored.
CellParameters(::BravaisLattice{Tetragonal}, a, b, c, args...) =
    CellParameters(a, a, c, SymPy.PI / 2, SymPy.PI / 2, SymPy.PI / 2)  # `b` is ignored.
CellParameters(::BravaisLattice{Cubic}, a, args...) =
    CellParameters(a, a, a, SymPy.PI / 2, SymPy.PI / 2, SymPy.PI / 2)  # Only `a` matters.
CellParameters(::BravaisLattice{Hexagonal{3},Primitive}, a, b, c, args...) =
    CellParameters(a, a, c, SymPy.PI / 2, SymPy.PI / 2, 2 * SymPy.PI / 3)  # `b` is ignored.
CellParameters(::BravaisLattice{Hexagonal{3},RhombohedralCentered}, a, b, c, α, args...) =
    CellParameters(a, a, a, α, α, α)  # `b`, `c` are ignored.

struct Lattice{T} <: AbstractMatrix{T}
    m::SMatrix{3,3,T}
end
Lattice(m::AbstractMatrix) = Lattice(SMatrix{3,3}(m))
Lattice(v1::AbstractVector, v2::AbstractVector, v3::AbstractVector) =
    vcat(transpose.((v1, v2, v3))...)

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

function makelattice(p::CellParameters)
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
end # function makelattice

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
