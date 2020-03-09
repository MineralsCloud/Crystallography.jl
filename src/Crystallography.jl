module Crystallography

using LinearAlgebra: Diagonal, cross, det, dot, norm

using CoordinateTransformations: Transformation, IdentityTransformation
using StaticArrays: FieldVector, SVector, SMatrix, SHermitianCompact, Size
using SymPy

import LinearAlgebra
import StaticArrays
import CoordinateTransformations

export CrystalSystem,
    Triclinic,
    Monoclinic,
    Orthorhombic,
    Tetragonal,
    Cubic,
    Trigonal,
    Hexagonal,
    Centering,
    BaseCentering,
    Primitive,
    BodyCentering,
    FaceCentering,
    RhombohedralCentering,
    BaseCentering,
    BravaisLattice,
    BravaisLattice,
    PrimitiveTriclinic,
    PrimitiveMonoclinic,
    ACenteredMonoclinic,
    BCenteredMonoclinic,
    CCenteredMonoclinic,
    PrimitiveOrthorhombic,
    ACenteredOrthorhombic,
    BCenteredOrthorhombic,
    CCenteredOrthorhombic,
    BodyCenteredOrthorhombic,
    FaceCenteredOrthorhombic,
    PrimitiveTetragonal,
    BodyCenteredTetragonal,
    PrimitiveCubic,
    BodyCenteredCubic,
    FaceCenteredCubic,
    PrimitiveHexagonal,
    RCenteredHexagonal,
    RealSpace,
    ReciprocalSpace,
    Crystal,
    MetricTensor,
    Miller,
    MillerBravais,
    Cell,
    LatticeConstants,
    AxisAngles,
    CellParameters,
    Lattice,
    RealFromReciprocal,
    ReciprocalFromReal,
    CrystalFromCartesian,
    CartesianFromCrystal,
    CrystalFromCrystal
export pearsonsymbol,
    arithmeticclass,
    centering,
    crystalsystem,
    directioncosine,
    directionangle,
    distance,
    interplanar_spacing,
    cellvolume,
    reciprocal,
    @m_str,
    @mb_str

abstract type CrystalSystem end
struct Triclinic <: CrystalSystem end
struct Monoclinic <: CrystalSystem end
struct Orthorhombic <: CrystalSystem end
struct Tetragonal <: CrystalSystem end
struct Cubic <: CrystalSystem end
struct Trigonal <: CrystalSystem end
struct Hexagonal <: CrystalSystem end

abstract type Centering end
struct Primitive <: Centering end
struct BodyCentering <: Centering end
struct FaceCentering <: Centering end
struct RhombohedralCentering <: Centering end
struct BaseCentering{T} <: Centering end
BaseCentering(T::Symbol) = T ∈ (:A, :B, :C) ? BaseCentering{T}() :
    throw(ArgumentError("centering must be either :A, :B, or :C!"))

const BravaisLattice = Tuple{CrystalSystem,Centering}
const PrimitiveTriclinic = Tuple{Triclinic,Primitive}
const PrimitiveMonoclinic = Tuple{Monoclinic,Primitive}
const ACenteredMonoclinic = Tuple{Monoclinic,BaseCentering{:A}}
const BCenteredMonoclinic = Tuple{Monoclinic,BaseCentering{:B}}
const CCenteredMonoclinic = Tuple{Monoclinic,BaseCentering{:C}}
const PrimitiveOrthorhombic = Tuple{Orthorhombic,Primitive}
const ACenteredOrthorhombic = Tuple{Orthorhombic,BaseCentering{:A}}
const BCenteredOrthorhombic = Tuple{Orthorhombic,BaseCentering{:B}}
const CCenteredOrthorhombic = Tuple{Orthorhombic,BaseCentering{:C}}
const BodyCenteredOrthorhombic = Tuple{Orthorhombic,BodyCentering}
const FaceCenteredOrthorhombic = Tuple{Orthorhombic,FaceCentering}
const PrimitiveTetragonal = Tuple{Tetragonal,Primitive}
const BodyCenteredTetragonal = Tuple{Tetragonal,BodyCentering}
const PrimitiveCubic = Tuple{Cubic,Primitive}
const BodyCenteredCubic = Tuple{Cubic,BodyCentering}
const FaceCenteredCubic = Tuple{Cubic,FaceCentering}
const PrimitiveHexagonal = Tuple{Hexagonal,Primitive}
const RCenteredHexagonal = Tuple{Hexagonal,RhombohedralCentering}
(::Type{Tuple{A,B}})() where {A<:CrystalSystem,B<:Centering} = (A(), B())

pearsonsymbol(::Triclinic) = "a"
pearsonsymbol(::Monoclinic) = "m"
pearsonsymbol(::Orthorhombic) = "o"
pearsonsymbol(::Tetragonal) = "t"
pearsonsymbol(::Cubic) = "c"
pearsonsymbol(::Hexagonal) = "h"
pearsonsymbol(::Trigonal) = "h"
pearsonsymbol(::Primitive) = "P"
pearsonsymbol(::BaseCentering{T}) where {T} = string(T)
pearsonsymbol(::BodyCentering) = "I"
pearsonsymbol(::FaceCentering) = "F"
pearsonsymbol(::RhombohedralCentering) = "R"
pearsonsymbol(b::BravaisLattice) =
    pearsonsymbol(crystalsystem(b)) * pearsonsymbol(centering(b))

arithmeticclass(::Triclinic) = "-1"
arithmeticclass(::Monoclinic) = "2/m"
arithmeticclass(::Orthorhombic) = "mmm"
arithmeticclass(::Tetragonal) = "4/mmm"
arithmeticclass(::Hexagonal) = "6/mmm"
arithmeticclass(::Cubic) = "m-3m"
arithmeticclass(::Primitive) = "P"
arithmeticclass(::BaseCentering) = "S"
arithmeticclass(::BodyCentering) = "I"
arithmeticclass(::FaceCentering) = "F"
arithmeticclass(::RhombohedralCentering) = "R"
arithmeticclass(b::BravaisLattice) =
    arithmeticclass(crystalsystem(b)) * arithmeticclass(centering(b))
arithmeticclass(::RCenteredHexagonal) = "-3mR"

centering(b::BravaisLattice) = last(b)

crystalsystem(b::BravaisLattice) = first(b)

const TETRAGONAL = Union{PrimitiveTetragonal,BodyCenteredTetragonal}
const CUBIC = Union{PrimitiveCubic,BodyCenteredCubic,FaceCenteredCubic}
const ORTHORHOMBIC = Union{
    PrimitiveOrthorhombic,
    BCenteredOrthorhombic,
    CCenteredOrthorhombic,
    BodyCenteredOrthorhombic,
    FaceCenteredCubic,
}
const MONOCLINIC = Union{PrimitiveMonoclinic,BCenteredMonoclinic,CCenteredMonoclinic}

struct RealFromReciprocal <: Transformation end
struct ReciprocalFromReal <: Transformation end
struct CartesianFromCrystal <: Transformation end
struct CrystalFromCartesian <: Transformation end
struct CrystalFromCrystal <: Transformation end

abstract type AbstractSpace end
struct RealSpace <: AbstractSpace end
struct ReciprocalSpace <: AbstractSpace end

struct Crystal{T} <: FieldVector{3,T}
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
LatticeConstants(::Union{PrimitiveTriclinic,MONOCLINIC,ORTHORHOMBIC}, a, b, c) =
    LatticeConstants(a, b, c)
LatticeConstants(::TETRAGONAL, a, b, c) = LatticeConstants(a, a, c)
LatticeConstants(::CUBIC, a, b, c) = LatticeConstants(a, a, a)
LatticeConstants(::PrimitiveHexagonal, a, b, c) = LatticeConstants(a, a, c)
LatticeConstants(::RCenteredHexagonal, a, b, c) = LatticeConstants(a, a, a)

struct AxisAngles{T} <: FieldVector{3,T}
    α::T
    β::T
    γ::T
end
AxisAngles(::PrimitiveTriclinic, α, β, γ) = AxisAngles(α, β, γ)
AxisAngles(::PrimitiveMonoclinic, α, β, γ, view::Int = 1) =
    view == 1 ? AxisAngles(90, 90, γ) : AxisAngles(90, β, 90)
AxisAngles(::CCenteredMonoclinic, α, β, γ) = AxisAngles(90, 90, γ)
AxisAngles(::BCenteredMonoclinic, α, β, γ) = AxisAngles(90, β, 90)
AxisAngles(::T, α, β, γ) where {T<:Union{ORTHORHOMBIC,TETRAGONAL,CUBIC}} =
    AxisAngles(90, 90, 90)
AxisAngles(::PrimitiveHexagonal, α, β, γ) = AxisAngles(90, 90, 120)
AxisAngles(::RCenteredHexagonal, α, β, γ) = AxisAngles(α, α, α)

struct CellParameters{S,T}
    data::NamedTuple{(:a, :b, :c, :α, :β, :γ),Tuple{S,S,S,T,T,T}}
end
function CellParameters(a, b, c, α, β, γ)
    a, b, c = promote(a, b, c)
    α, β, γ = promote(α, β, γ)
    return CellParameters((a = a, b = b, c = c, α = α, β = β, γ = γ))
end
CellParameters(a::LatticeConstants, b::AxisAngles) = CellParameters(a..., b...)
CellParameters(x::BravaisLattice) = args -> CellParameters(x, args...)
CellParameters(x::BravaisLattice, a, b, c, α, β, γ) =
    CellParameters(LatticeConstants(x, a, b, c), AxisAngles(x, α, β, γ))

struct Lattice{T}
    data::SVector{3,SVector{3,T}}
end
Lattice(v1::AbstractVector, v2::AbstractVector, v3::AbstractVector) =
    Lattice(SVector(map(SVector{3}, (v1, v2, v3))))
function Lattice(m::AbstractMatrix, rowmajor::Bool = false)
    f = rowmajor ? transpose : identity
    return Lattice(Iterators.partition(f(m), 3)...)
end # function Lattice
function Lattice(a, b, c, α, β, γ)
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
Lattice(p::CellParameters) = Lattice(p...)

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
Cell(lattice::Lattice, positions, numbers, args...) =
    Cell(lattice.data, positions, numbers, args...)

struct MetricTensor{T} <: AbstractMatrix{T}
    data::SHermitianCompact{3,T}
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
function crystalsystem(lattice::Lattice)
    v1, v2, v3 = lattice.data
    a, b, c = norm(v1), norm(v2), norm(v3)
    γ = acos(dot(v1, v2) / a / b)
    β = acos(dot(v2, v3) / b / c)
    α = acos(dot(v1, v3) / a / c)
    return crystalsystem(CellParameters(a, b, c, α, β, γ))
end # function crystalsystem

"""
    cellvolume(a, b, c, α, β, γ)
    cellvolume(p::CellParameters)

Calculates the cell volume from 6 cell parameters.
"""
cellvolume(a, b, c, α, β, γ) =
    a * b * c * sqrt(sin(α)^2 - cos(β)^2 - cos(γ)^2 + 2 * cos(α) * cos(β) * cos(γ))
cellvolume(p::CellParameters) = cellvolume(p...)
"""
    cellvolume(l::Lattice)
    cellvolume(c::Cell)

Calculates the cell volume from a `Lattice` or a `Cell`.
"""
cellvolume(l::Lattice) = abs(det(convert(Matrix{eltype(l)}, l)))
cellvolume(c::Cell) = cellvolume(c.lattice)
"""
    cellvolume(g::MetricTensor)

Calculates the cell volume from a `MetricTensor`.
"""
cellvolume(g::MetricTensor) = sqrt(det(g.data))  # `sqrt` is always positive!

function reciprocal(lattice::Lattice, twopi::Bool = false)
    volume = cellvolume(lattice)
    a1, a2, a3 = lattice.data
    factor = twopi ? 2 * SymPy.PI : 1
    return factor / volume * [cross(a2, a3) cross(a3, a1) cross(a1, a2)]
end # function reciprocal

directioncosine(a::Crystal, g::MetricTensor, b::Crystal) =
    dot(a, g, b) / (norm(a, g) * norm(b, g))

directionangle(a::Crystal, g::MetricTensor, b::Crystal) = acos(directioncosine(a, g, b))

distance(a::Crystal, g::MetricTensor, b::Crystal) = norm(b - a, g)

interplanar_spacing(a::Crystal, g::MetricTensor) = 1 / norm(a, g)

"""
    supercell(cell::Lattice, expansion::AbstractMatrix{<:Integer})

Allow the supercell to be a tilted extension of `cell`.
"""
function supercell(cell::Lattice, expansion::AbstractMatrix{<:Integer})
    @assert(det(expansion) != 0, "matrix `expansion` cannot be a singular integer matrix!")
    return expansion * cell
end # function supercell
"""
    supercell(cell::Lattice, expansion::AbstractVector{<:Integer})

Return a supercell based on `cell` and expansion coefficients.
"""
function supercell(cell::Lattice, expansion::AbstractVector{<:Integer})
    @assert length(expansion) == 3
    return supercell(cell, Diagonal(expansion))
end # function supercell

Base.size(::Union{MetricTensor,Lattice}) = (3, 3)
Base.size(::Miller) = (3,)
Base.size(::MillerBravais) = (4,)
Base.size(::CellParameters) = (6,)

Base.getindex(A::MetricTensor, I::Vararg{Int}) = getindex(A.data, I...)
Base.getindex(
    A::Union{Miller,MillerBravais,CellParameters,Lattice,CellParameters},
    i::Int,
) = getindex(A.data, i)
Base.getindex(A::Lattice, i::Int, j::Int) = getindex(getindex(A.data, i), j)

Base.inv(g::MetricTensor) = MetricTensor(inv(SymPy.N(g.data)))
Base.inv(::RealFromReciprocal) = ReciprocalFromReal()
Base.inv(::ReciprocalFromReal) = RealFromReciprocal()
Base.inv(::CartesianFromCrystal) = CrystalFromCartesian()
Base.inv(::CrystalFromCartesian) = CartesianFromCrystal()

CoordinateTransformations.compose(::RealFromReciprocal, ::ReciprocalFromReal) =
    IdentityTransformation()
CoordinateTransformations.compose(::ReciprocalFromReal, ::RealFromReciprocal) =
    IdentityTransformation()
CoordinateTransformations.compose(::CrystalFromCartesian, ::CartesianFromCrystal) =
    IdentityTransformation()
CoordinateTransformations.compose(::CartesianFromCrystal, ::CrystalFromCartesian) =
    IdentityTransformation()

(::CrystalFromCartesian)(to::Lattice) = convert(Matrix{eltype(to)}, to)
(::CartesianFromCrystal)(from::Lattice) = convert(Matrix{eltype(from)}, from)'
(::CrystalFromCrystal)(to::Lattice, from::Lattice) =
    convert(Matrix{eltype(to)}, to) * convert(Matrix{eltype(from)}, from)

Base.convert(::Type{Crystal}, lattice::Lattice, v::AbstractVector) =
    Crystal(CrystalFromCartesian()(lattice) * v)
Base.convert(::Type{T}, lattice::Lattice, v::Crystal) where {T<:AbstractVector} =
    T(CartesianFromCrystal()(lattice) * collect(v))
Base.convert(::Type{Crystal}, from::Lattice, to::Lattice, v::Crystal) =
    CrystalFromCrystal()(to, from) * v
Base.convert(::Type{Matrix{T}}, lattice::Lattice{T}) where {T} = hcat(lattice.data...)
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
    a2, b2, c2, ab, ac, bc = data[1, 1], data[2, 2], data[3, 3], data[1, 2], data[1, 3], data[2, 3]
    a, b, c = map(sqrt, (a2, b2, c2))
    γ, β, α = acos(ab / (a * b)), acos(ac / (a * c)), acos(bc / (b * c))
    return CellParameters(a, b, c, α, β, γ)
end # function Base.convert
Base.convert(::Type{Lattice}, g::MetricTensor) = Lattice(convert(CellParameters, g))

Base.iterate(c::CellParameters, args...) = iterate(c.data, args...)

Base.eltype(::Lattice{T}) where {T} = T

Base.firstindex(::CellParameters) = 1

Base.lastindex(::CellParameters) = 6

Base.getproperty(p::CellParameters, name::Symbol) =
    name ∈ (:a, :b, :c, :α, :β, :γ) ? getfield(p.data, name) : getfield(p, name)

LinearAlgebra.dot(a::Crystal, g::MetricTensor, b::Crystal) = a' * g.data * b
LinearAlgebra.norm(a::Crystal, g::MetricTensor) = sqrt(dot(a, g, a))

StaticArrays.similar_type(
    ::Type{<:Crystal},  # Do not delete the `<:`!
    ::Type{T},
    size::Size{(3,)},
) where {T} = Crystal{T}

include("Symmetry.jl")

end # module
