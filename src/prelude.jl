using LinearAlgebra: Symmetric, cross, det, dot, norm

using AutoHashEquals: @auto_hash_equals
using CoordinateTransformations
using StaticArrays: FieldVector, Size
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
    CellParameters

export pearsonsymbol,
    arithmeticclass,
    centeringof,
    crystalsystem,
    dimensionof,
    axessetting,
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
function Hexagonal(N::Int = 3)
    @assert N ∈ (2, 3)
    return Hexagonal{N}()
end

abstract type Centering end
struct Primitive <: Centering end
struct BodyCentered <: Centering end
struct FaceCentered <: Centering end
struct RhombohedralCentered <: Centering end
struct BaseCentered{T} <: Centering end
function BaseCentered(T::Symbol)
    @assert T ∈ (:A, :B, :C)
    return BaseCentered{T}()
end # function BaseCentered

struct BravaisLattice{A<:CrystalSystem,B<:Centering,N} end
BravaisLattice(::A, ::B) where {A,B} = BravaisLattice{A,B,1}()
BravaisLattice(::A, ::B, N::Integer) where {A,B} = BravaisLattice{A,B,N}()
BravaisLattice(ibrav::Integer) = BravaisLattice(Val(ibrav))
BravaisLattice(::Val{1}) = BravaisLattice(Cubic(), Primitive())
BravaisLattice(::Val{2}) = BravaisLattice(Cubic(), FaceCentered())
BravaisLattice(::Val{3}) = BravaisLattice(Cubic(), BodyCentered())
BravaisLattice(::Val{4}) = BravaisLattice(Hexagonal(), Primitive())
BravaisLattice(::Val{5}) = BravaisLattice(Hexagonal(), RhombohedralCentered())
BravaisLattice(::Val{-5}) = BravaisLattice(Hexagonal(), RhombohedralCentered())
BravaisLattice(::Val{6}) = BravaisLattice(Tetragonal(), Primitive())
BravaisLattice(::Val{7}) = BravaisLattice(Tetragonal(), BodyCentered())
BravaisLattice(::Val{8}) = BravaisLattice(Orthorhombic(), Primitive())
BravaisLattice(::Val{9}) = BravaisLattice(Orthorhombic(), BaseCentered(:B))
BravaisLattice(::Val{-9}) = BravaisLattice(Orthorhombic(), BaseCentered(:C))
BravaisLattice(::Val{91}) = BravaisLattice(Orthorhombic(), BaseCentered(:A))  # New in QE 6.5
BravaisLattice(::Val{10}) = BravaisLattice(Orthorhombic(), FaceCentered())
BravaisLattice(::Val{11}) = BravaisLattice(Orthorhombic(), BodyCentered())
BravaisLattice(::Val{12}) = BravaisLattice(Monoclinic(), Primitive())
BravaisLattice(::Val{-12}) = BravaisLattice(Monoclinic(), Primitive())
BravaisLattice(::Val{13}) = BravaisLattice(Monoclinic(), BaseCentered(:B))
BravaisLattice(::Val{14}) = BravaisLattice(Triclinic(), Primitive())

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

axessetting(::BravaisLattice{C,T,N}) where {C,T,N} = N

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
function CellParameters(a, b, c, α, β, γ, angle_iscosine::Bool = false)
    v = angle_iscosine ? (a, b, c, acos(α), acos(β), acos(γ)) : (a, b, c, α, β, γ)
    return CellParameters{Base.promote_typeof(v...)}(v...)
end
CellParameters(bravais::BravaisLattice) = args -> CellParameters(bravais, args...)
CellParameters(::BravaisLattice{Triclinic}, a, b, c, α, β, γ, args...) =
    CellParameters(a, b, c, α, β, γ)  # Triclinic
CellParameters(::BravaisLattice{Monoclinic,Primitive,1}, a, b, c, α, β, γ, args...) =
    CellParameters(a, b, c, SymPy.PI / 2, SymPy.PI / 2, γ)  # `α`, `β` are ignored.
CellParameters(::BravaisLattice{Monoclinic,Primitive,2}, a, b, c, α, β, γ, args...) =
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

@auto_hash_equals struct MetricTensor{T<:AbstractMatrix}
    m::T
    function MetricTensor{T}(m) where {T}
        @assert(size(m) == (3, 3), "The metric tensor must be of size 3×3!")
        return new(m)
    end
end
MetricTensor(m::T) where {T} = MetricTensor{T}(m)
function MetricTensor(v1::AbstractVector, v2::AbstractVector, v3::AbstractVector)
    vecs = (v1, v2, v3)
    return MetricTensor(map(x -> dot(x...), Iterators.product(vecs, vecs)))
end
function MetricTensor(a, b, c, α, β, γ)
    g12 = a * b * cos(γ)
    g13 = a * c * cos(β)
    g23 = b * c * cos(α)
    return MetricTensor(Symmetric([
        a^2 g12 g13
        g12 b^2 g23
        g13 g23 c^2
    ]))
end
MetricTensor(p::CellParameters) = MetricTensor(p...)

directioncosine(a::CrystalCoordinates, g::MetricTensor, b::CrystalCoordinates) =
    dot(a, g, b) / (norm(a, g) * norm(b, g))

directionangle(a::CrystalCoordinates, g::MetricTensor, b::CrystalCoordinates) =
    acos(directioncosine(a, g, b))

distance(a::CrystalCoordinates, g::MetricTensor, b::CrystalCoordinates) = norm(b - a, g)

interplanar_spacing(a::CrystalCoordinates, g::MetricTensor) = 1 / norm(a, g)

function reciprocalof(mat::AbstractMatrix, twopi::Bool = false)
    @assert size(mat) == (3, 3)
    volume = abs(det(mat))
    a1, a2, a3 = mat[1, :], mat[2, :], mat[3, :]
    factor = twopi ? 2 * SymPy.PI : 1
    return factor / volume * [cross(a2, a3) cross(a3, a1) cross(a1, a2)]
end # function reciprocalof

struct MillerIndices{S<:AbstractSpace,T<:Integer}
    v::NTuple{3,T}
    function MillerIndices{S,T}(x) where {S,T}
        y = collect(x)
        return new(iszero(y) ? x : x .÷ gcd(y))
    end
end
MillerIndices{S}(x::NTuple{3,T}) where {S,T} = MillerIndices{S,T}(x)
MillerIndices{S}(i, j, k) where {S} = MillerIndices{S}((i, j, k))
MillerIndices{S}(x::AbstractVector) where {S} = MillerIndices{S}(Tuple(x))

struct MillerBravaisIndices{S<:AbstractSpace,T<:Integer}
    v::NTuple{4,T}
    function MillerBravaisIndices{S,T}(x) where {S,T}
        y = collect(x)
        return new(iszero(y) ? x : x .÷ gcd(y))
    end
end
MillerBravaisIndices{S}(x::NTuple{4,T}) where {S,T} = MillerBravaisIndices{S,T}(x)
MillerBravaisIndices{S}(i, j, k, l) where {S} = MillerBravaisIndices{S}((i, j, k, l))
MillerBravaisIndices{S}(x::AbstractVector) where {S} = MillerBravaisIndices{S}(Tuple(x))

CrystalCoordinates(m::MillerIndices) = CrystalCoordinates(m.i, m.j, m.k)
CrystalCoordinates(mb::MillerBravaisIndices{T}) where {T} =
    CrystalCoordinates(convert(MillerIndices{T}, mb))

function makelattice(b::BravaisLattice, params...; vecform::Bool = false, view::Int = 1)
    lattice = makelattice(b, CellParameters(b, params...))
    return vecform ? _splitlattice(lattice) : lattice
end # function makelattice
makelattice(::BravaisLattice{Cubic,Primitive}, cell::CellParameters) =
    cell[1] * [
        1 0 0
        0 1 0
        0 0 1
    ]
makelattice(::BravaisLattice{Cubic,FaceCentered}, cell::CellParameters) =
    cell[1] / 2 * [
        -1 0 1
        0 1 1
        -1 1 0
    ]
function makelattice(
    ::BravaisLattice{Cubic,BodyCentered},
    cell::CellParameters,
    view::Int = 1,
)
    if view == 1
        cell[1] / 2 * [
            1 1 1
            -1 1 1
            -1 -1 1
        ]
    elseif view == 2
        cell[1] / 2 * [
            -1 1 1
            1 -1 1
            1 1 -1
        ]
    else
        error("wrong `view` $view input!")
    end
end # function makelattice
makelattice(::BravaisLattice{Hexagonal{3},Primitive}, cell::CellParameters) =
    cell[1] * [
        1 0 0
        -1 / 2 √3 / 2 0
        0 0 cell[3] / cell[1]
    ]
function makelattice(
    ::BravaisLattice{Hexagonal{3},RhombohedralCentered},
    cell::CellParameters,
    view::Int = 1,
)
    if view == 1
        r = cos(cell[4])
        tx = sqrt((1 - r) / 2)
        ty = sqrt((1 - r) / 6)
        tz = sqrt((1 + 2r) / 3)
        cell[1] * [
            tx -ty tz
            0 2ty tz
            -tx -ty tz
        ]
    elseif view == 2
        ap = cell[1] / √3
        c = acos(cell[4])
        ty = sqrt((1 - c) / 6)
        tz = sqrt((1 + 2c) / 3)
        u = tz - 2 * √2 * ty
        v = tz + √2 * ty
        ap * [
            u v v
            v u v
            v v u
        ]
    else
        error("wrong `view` $view input!")
    end
end
makelattice(::BravaisLattice{Tetragonal,Primitive}, cell::CellParameters) =
    cell[1] * [
        1 0 0
        0 1 0
        0 0 cell[3] / cell[1]
    ]
function makelattice(::BravaisLattice{Tetragonal,BodyCentered}, cell::CellParameters)
    r = cell[3] / cell[1]
    return cell[1] / 2 * [
        1 -1 r
        1 1 r
        -1 -1 r
    ]
end
makelattice(::BravaisLattice{Orthorhombic,Primitive}, cell::CellParameters) = [
    cell[1] 0 0
    0 cell[2] 0
    0 0 cell[3]
]
# TODO: BravaisLattice{Orthorhombic,CCentered}
function makelattice(::BravaisLattice{Orthorhombic,FaceCentered}, cell::CellParameters)
    a, b, c = cell[1:3]
    return [
        a 0 c
        a b 0
        0 b c
    ] / 2
end
function makelattice(::BravaisLattice{Orthorhombic,BodyCentered}, cell::CellParameters)
    a, b, c = cell[1:3]
    return [
        a b c
        -a b c
        -a -b c
    ] / 2
end
function makelattice(::BravaisLattice{Monoclinic,Primitive}, cell::CellParameters)
    a, b, c = cell[1:3]
    return [
        a 0 0
        0 b 0
        c * cos(cell[5]) 0 c * sin(cell[5])
    ]
end
# TODO: BravaisLattice{Monoclinic,BCentered}
# TODO: BravaisLattice{Triclinic,Primitive}

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

Base.inv(g::MetricTensor) = MetricTensor(inv(SymPy.N(g.m)))

Base.convert(::Type{<:MillerIndices}, mb::MillerBravaisIndices) =
    error("unsupported conversion!")
Base.convert(::Type{<:MillerBravaisIndices}, m::MillerIndices) =
    error("unsupported conversion!")
Base.convert(T::Type{<:MillerIndices}, m::MillerIndices) =
    isa(m, T) ? m : error("unsupported conversion!")
Base.convert(T::Type{<:MillerBravaisIndices}, mb::MillerBravaisIndices) =
    isa(mb, T) ? mb : error("unsupported conversion!")
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

LinearAlgebra.dot(a::CrystalCoordinates, g::MetricTensor, b::CrystalCoordinates) =
    a' * g.m * b
LinearAlgebra.norm(a::CrystalCoordinates, g::MetricTensor) = sqrt(dot(a, g, a))
