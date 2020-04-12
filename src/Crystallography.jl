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
    CrystalCoord,
    MetricTensor,
    Miller,
    MillerBravais,
    AtomicPosition,
    Cell,
    latticeconst,
    axisangles,
    CellParameters,
    Lattice,
    RealFromReciprocal,
    ReciprocalFromReal,
    CrystalFromCartesian,
    CartesianFromCrystal,
    CrystalFromCrystal
export centering,
    crystalsystem,
    directioncosine,
    directionangle,
    distance,
    interplanar_spacing,
    cellvolume,
    reciprocal,
    eachatom,
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
BaseCentering(T::Symbol) = T ∈ (:A, :B, :C) ? BaseCentering{T}() : throw(ArgumentError("centering must be either :A, :B, or :C!"))

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

function latticeconst(a, b, c)
    x = Base.vect(a, b, c)  # Do a `promote_typeof`, better than just a `Tuple`
    # See the difference between `eltype(1u"m", 2.0u"m", 3.0, 4)` & `eltype(Base.vect(1u"m", 2.0u"m", 3.0, 4))`.
    z = zero(eltype(x))
    all(x .> z) ? x : error("lattice constants must all be positive!")
end # function latticeconst
latticeconst(::Union{PrimitiveTriclinic,MONOCLINIC,ORTHORHOMBIC}, a, b, c) = latticeconst(a, b, c)
latticeconst(::Union{TETRAGONAL,PrimitiveHexagonal}, a, c) = latticeconst(a, a, c)
latticeconst(::Union{CUBIC,RCenteredHexagonal}, a) = latticeconst(a, a, a)

axisangles(::PrimitiveTriclinic, α, β, γ) = axisangles(α, β, γ)
axisangles(::PrimitiveMonoclinic, α, θ, view::Int = 1) =
    view == 1 ? axisangles(90, 90, θ) : axisangles(90, θ, 90)
axisangles(::CCenteredMonoclinic, γ) = axisangles(90, 90, γ)
axisangles(::BCenteredMonoclinic, β) = axisangles(90, β, 90)
axisangles(::Union{ORTHORHOMBIC,TETRAGONAL,CUBIC}) = axisangles(90, 90, 90)
axisangles(::PrimitiveHexagonal) = axisangles(90, 90, 120)
axisangles(::RCenteredHexagonal, α) = axisangles(α, α, α)

const CellParameters = NamedTuple{(:a, :b, :c, :α, :β, :γ)}  # Use as type
CellParameters(a, b, c, α, β, γ) = (a = a, b = b, c = c, α = α, β = β, γ = γ)  # Use as constructor

struct Lattice{T}
    data::SMatrix{3,3,T}
end
Lattice(m::AbstractMatrix) = Lattice(SMatrix{3,3}(m))
Lattice(a::AbstractVector, b::AbstractVector, c::AbstractVector) = Lattice(transpose(hcat(a, b, c)))
Lattice(v::AbstractVector{<:AbstractVector}) = Lattice(v...)
function Lattice(p::CellParameters)
    a, b, c, α, β, γ = p
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

struct AtomicPosition{S,T}
    atom::S
    pos::SVector{3,T}
end
AtomicPosition(atom, pos::AbstractVector) = AtomicPosition(atom, SVector{3}(pos))

struct Cell{N,L,S,T}
    atompos::SVector{N,AtomicPosition{S,T}}
    lattice::Lattice{L}
end
function Cell(atoms::AbstractVector, positions::AbstractVector{<:AbstractVector}, lattice::AbstractVecOrMat)
    if length(positions) == length(atoms)
        N = length(positions)
        return Cell(SVector{N}([AtomicPosition(atom, position) for (atom, position) in zip(atoms, positions)]), Lattice(lattice))
    else
        throw(DimensionMismatch("the number of positions should equal the number of atoms!"))
    end
end # function Cell

# This is an internal type and should not be exported!
struct AtomicIterator{T}
    data::T
end

eachatom(atompos::AbstractVector{<:AtomicPosition}) = AtomicIterator(atompos)
eachatom(cell::Cell) = AtomicIterator(cell.atompos)

centering(b::BravaisLattice) = last(b)

crystalsystem(b::BravaisLattice) = first(b)
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

function reciprocal(lattice::Lattice, twopi::Bool = false)
    volume = cellvolume(lattice)
    a1, a2, a3 = lattice.data
    factor = twopi ? 2 * SymPy.PI : 1
    return factor / volume * [cross(a2, a3) cross(a3, a1) cross(a1, a2)]
end # function reciprocal

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

Base.size(::Lattice) = (3, 3)
Base.length(::Lattice) = 9  # Number of elements
Base.getindex(A::Lattice, i::Integer, j::Integer) = getindex(A.data, i, j)
Base.eltype(::Lattice{T}) where {T} = T
Base.convert(::Type{Matrix{T}}, lattice::Lattice{T}) where {T} = Matrix{T}(transpose(lattice.data))  # Use with care!

Base.length(iter::AtomicIterator) = length(iter.data)
Base.size(iter::AtomicIterator) = (length(iter.data),)
Base.iterate(iter::AtomicIterator{<:AbstractVector{<:AtomicPosition}}, i = 1) = i > length(iter) ? nothing : (iter.data[i], i + 1)
Base.eltype(::AtomicIterator{<:AbstractVector{T}}) where {T<:AtomicPosition} = T

include("transform.jl")
include("geometry.jl")
include("Symmetry.jl")

end # module
