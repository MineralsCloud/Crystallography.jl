module Crystallography

using LinearAlgebra: Diagonal, det, dot, norm
using Spglib: Cell
using StaticArrays: SVector, SMatrix
using Unitful: AbstractQuantity, ustrip, unit

import Spglib: basis_vectors

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
    Bravais,
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
    Cell,
    Lattice
export centering, crystalsystem, cellvolume, basis_vectors

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
const ACentering = BaseCentering{:A}
const BCentering = BaseCentering{:B}
const CCentering = BaseCentering{:C}

struct Bravais{A<:CrystalSystem,B<:Centering}
    obverse::Bool
end
Bravais(a::CrystalSystem, b::Centering, obverse::Bool = true) =
    Bravais{typeof(a),typeof(b)}(obverse)

const PrimitiveTriclinic = Bravais{Triclinic,Primitive}
const PrimitiveMonoclinic = Bravais{Monoclinic,Primitive}
const ACenteredMonoclinic = Bravais{Monoclinic,ACentering}
const BCenteredMonoclinic = Bravais{Monoclinic,BCentering}
const CCenteredMonoclinic = Bravais{Monoclinic,CCentering}
const PrimitiveOrthorhombic = Bravais{Orthorhombic,Primitive}
const ACenteredOrthorhombic = Bravais{Orthorhombic,ACentering}
const BCenteredOrthorhombic = Bravais{Orthorhombic,BCentering}
const CCenteredOrthorhombic = Bravais{Orthorhombic,CCentering}
const BodyCenteredOrthorhombic = Bravais{Orthorhombic,BodyCentering}
const FaceCenteredOrthorhombic = Bravais{Orthorhombic,FaceCentering}
const PrimitiveTetragonal = Bravais{Tetragonal,Primitive}
const BodyCenteredTetragonal = Bravais{Tetragonal,BodyCentering}
const PrimitiveCubic = Bravais{Cubic,Primitive}

const BodyCenteredCubic = Bravais{Cubic,BodyCentering}
const FaceCenteredCubic = Bravais{Cubic,FaceCentering}
const PrimitiveHexagonal = Bravais{Hexagonal,Primitive}
const RCenteredHexagonal = Bravais{Hexagonal,RhombohedralCentering}

const TETRAGONAL = Union{PrimitiveTetragonal,BodyCenteredTetragonal}
const CUBIC = Union{PrimitiveCubic,BodyCenteredCubic,FaceCenteredCubic}
const ORTHORHOMBIC = Union{
    PrimitiveOrthorhombic,
    BCenteredOrthorhombic,
    CCenteredOrthorhombic,
    BodyCenteredOrthorhombic,
    FaceCenteredOrthorhombic,
}
const MONOCLINIC = Union{PrimitiveMonoclinic,BCenteredMonoclinic,CCenteredMonoclinic}

struct Lattice{T}
    data::SMatrix{3,3,T}
end
Lattice(m::AbstractMatrix) = Lattice(SMatrix{3,3}(m))
Lattice(a::AbstractVector, b::AbstractVector, c::AbstractVector) = Lattice(hcat(a, b, c))
Lattice(x::Lattice) = x

function basis_vectors(lattice::Lattice)
    data = lattice.data
    return data[:, 1], data[:, 2], data[:, 3]
end

centering(::Bravais{A,B}) where {A,B} = B()

crystalsystem(::Bravais{A,B}) where {A,B} = A()
function crystalsystem(a, b, c, α, β, γ)
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
end
function crystalsystem(lattice::Lattice)
    v1, v2, v3 = basis_vectors(lattice)
    a, b, c = norm(v1), norm(v2), norm(v3)
    γ = acos(dot(v1, v2) / a / b)
    β = acos(dot(v2, v3) / b / c)
    α = acos(dot(v1, v3) / a / c)
    return crystalsystem(a, b, c, α, β, γ)
end

"""
    supercell(cell::Lattice, expansion::AbstractMatrix{<:Integer})

Allow the supercell to be a tilted extension of `cell`.
"""
function supercell(cell::Lattice, expansion::AbstractMatrix{<:Integer})
    @assert(det(expansion) != 0, "matrix `expansion` cannot be a singular integer matrix!")
    return expansion * cell
end
"""
    supercell(cell::Lattice, expansion::AbstractVector{<:Integer})

Return a supercell based on `cell` and expansion coefficients.
"""
function supercell(cell::Lattice, expansion::AbstractVector{<:Integer})
    @assert length(expansion) == 3
    return supercell(cell, Diagonal(expansion))
end

"""
    cellvolume(p::CellParameters)

Calculates the cell volume from 6 cell parameters.
"""
cellvolume(a, b, c, α, β, γ) =
    a * b * c * sqrt(sin(α)^2 - cos(β)^2 - cos(γ)^2 + 2 * cos(α) * cos(β) * cos(γ))
"""
    cellvolume(l::Lattice)
    cellvolume(c::Cell)

Calculates the cell volume from a `Lattice` or a `Cell`.
"""
cellvolume(lattice::Lattice) = abs(det(lattice.data))
cellvolume(cell::Cell) = cellvolume(cell.lattice)

Base.size(::Lattice) = (3, 3)
Base.length(::Lattice) = 9  # Number of elements
Base.getindex(A::Lattice, i::Integer, j::Integer) = getindex(A.data, i, j)
Base.eltype(::Lattice{T}) where {T} = T

include("Arithmetics.jl")
include("Symmetry.jl")
include("transform.jl")

end
