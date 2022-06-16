using LinearAlgebra: Diagonal, I

import Spglib: Cell

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
export centering, crystalsystem, basis_vectors, cellparameters

"Represent one of the seven crystal systems."
abstract type CrystalSystem end
"""
    Triclinic()

Represent the triclinic system.
"""
struct Triclinic <: CrystalSystem end
"""
    Monoclinic()

Represent the monoclinic system.
"""
struct Monoclinic <: CrystalSystem end
"""
    Orthorhombic()

Represent the orthorhombic system.
"""
struct Orthorhombic <: CrystalSystem end
"""
    Tetragonal()

Represent the tetragonal system.
"""
struct Tetragonal <: CrystalSystem end
"""
    Cubic()

Represent the cubic system.
"""
struct Cubic <: CrystalSystem end
"""
    Trigonal()

Represent the trigonal system.
"""
struct Trigonal <: CrystalSystem end
"""
    Hexagonal()

Represent the hexagonal system.
"""
struct Hexagonal <: CrystalSystem end

"Represent the centering types."
abstract type Centering end
"""
    Primitive()

Represent no centering.
"""
struct Primitive <: Centering end
"""
    BodyCentering()

Represent the body-centering.
"""
struct BodyCentering <: Centering end
"""
    FaceCentering()

Represent the face-centering.
"""
struct FaceCentering <: Centering end
"""
    RhombohedralCentering()

Represent the rhombohedral-centering of the hexagonal system.
"""
struct RhombohedralCentering <: Centering end
"""
    BaseCentering{:A}()
    BaseCentering{:B}()
    BaseCentering{:C}()

Represent the base-centering.
"""
struct BaseCentering{T} <: Centering end
const ACentering = BaseCentering{:A}
const BCentering = BaseCentering{:B}
const CCentering = BaseCentering{:C}

"""
    Bravais(a::CrystalSystem, b::Centering, obverse::Bool=true)

Represent a Bravais lattice type.
"""
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

"Represent the real lattices and the reciprocal lattices."
abstract type AbstractLattice{T} end

struct Lattice{T} <: AbstractLattice{T}
    data::SMatrix{3,3,T,9}
end
"""
    Lattice(mat::AbstractMatrix)

Construct a `Lattice` from a matrix.

!!! note
    The basis vectors of the matrix are stored as columns.
"""
Lattice(mat::AbstractMatrix) = Lattice(SMatrix{3,3}(mat))
"""
    Lattice(𝐚::AbstractVector, 𝐛::AbstractVector, 𝐜::AbstractVector)

Construct a `Lattice` from three basis vectors.
"""
Lattice(𝐚::AbstractVector, 𝐛::AbstractVector, 𝐜::AbstractVector) = Lattice(hcat(𝐚, 𝐛, 𝐜))
"""
    Lattice(cell::Cell)

Get the lattice of a `Cell`.
"""
Lattice(cell::Cell) = Lattice(cell.lattice)
"""
    Lattice(a, b, c, α, β, γ)

Construct a `Lattice` from the six cell parameters.

The convention we used here is that edge vector 𝐚 in the positive x-axis direction,
edge vector 𝐛 in the x-y plane with positive y-axis component,
and edge vector 𝐜 with positive z-axis component in the Cartesian-system.
See [Wikipedia](https://en.wikipedia.org/w/index.php?title=Fractional_coordinates&oldid=961675499#In_crystallography).
"""
function Lattice(a, b, c, α, β, γ)
    Ω = cellvolume(a, b, c, α, β, γ)
    sinγ, cosγ, cosα, cosβ = sind(γ), cosd(γ), cosd(α), cosd(β)
    return Lattice(
        [a, zero(a), zero(a)],
        [b * cosγ, b * sinγ, zero(b)],
        [c * cosβ, c * (cosα - cosβ * cosγ) / sinγ, Ω / (a * b * sinγ)],
    )
end
@functor Lattice

Cell(lattice::Lattice, positions, types, magmoms = zeros(length(types))) =
    Cell(lattice.data, positions, types, magmoms)

"""
    basis_vectors(lattice::Lattice)

Get the three basis vectors from a `lattice`.
"""
basis_vectors(lattice::Lattice) = lattice[:, 1], lattice[:, 2], lattice[:, 3]

"""
    centering(bravais::Bravais)

Get the centering type of a Bravais type.
"""
centering(::Bravais{A,B}) where {A,B} = B()

"""
    crystalsystem(bravais::Bravais)

Get the crystal system of a Bravais type.
"""
crystalsystem(::Bravais{A,B}) where {A,B} = A()
"""
    crystalsystem(a, b, c, α, β, γ)

Guess the crystal system from the six cell parameters.
"""
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
"""
    crystalsystem(lattice::Lattice)

Get the crystal system of a `lattice`.
"""
crystalsystem(lattice::Lattice) = crystalsystem(cellparameters(lattice)...)

"""
    cellparameters(lattice::Lattice)

Get the six cell parameters from a `lattice`.
"""
function cellparameters(lattice::Lattice)
    𝐚, 𝐛, 𝐜 = basis_vectors(lattice)
    a, b, c = norm(𝐚), norm(𝐛), norm(𝐜)
    γ, β, α =
        acosd(dot(𝐚, 𝐛) / (a * b)), acosd(dot(𝐚, 𝐜) / (a * c)), acosd(dot(𝐛, 𝐜) / (b * c))
    return a, b, c, α, β, γ
end

# See https://en.wikipedia.org/wiki/Supercell_(crystal)
"""
    supercell(lattice::Lattice, scaling_factors::AbstractMatrix{<:Integer})
    supercell(lattice::Lattice, scaling_factors::AbstractVector{<:Integer})
    supercell(lattice::Lattice, scaling_factor::Integer)

Create a supercell from `lattice`.
"""
function supercell(lattice::Lattice, scaling_factors::AbstractMatrix{<:Integer})
    @assert det(scaling_factors) >= 1
    return Lattice(lattice.data * scaling_factors)
end
supercell(lattice::Lattice, scaling_factors::AbstractVector{<:Integer}) =
    supercell(lattice, Diagonal(scaling_factors))
supercell(lattice::Lattice, scaling_factor::Integer) =
    supercell(lattice, scaling_factor * I)

Base.iterate(lattice::AbstractLattice) = iterate(lattice.data)
Base.iterate(lattice::AbstractLattice, state) = iterate(lattice.data, state)

Base.eltype(::AbstractLattice{T}) where {T} = T

Base.length(::AbstractLattice) = 9

Base.size(::AbstractLattice) = (3, 3)
Base.size(::AbstractLattice, dim::Integer) = dim <= 2 ? 3 : 1

Base.IteratorSize(::Type{<:AbstractLattice}) = Base.HasShape{2}()

Base.axes(lattice::AbstractLattice, dim::Integer) = axes(lattice.data, dim)

Base.getindex(lattice::AbstractLattice, i) = getindex(lattice.data, i)
Base.getindex(lattice::AbstractLattice, I::Vararg) = getindex(lattice.data, I...)

Base.firstindex(::AbstractLattice) = 1

Base.lastindex(::AbstractLattice) = 9

for op in (:+, :-)
    @eval Base.broadcast(::typeof($op), lattice::AbstractLattice, number::Number) =
        Lattice(broadcast($op, lattice.data, number))
    @eval Base.broadcast(::typeof($op), number::Number, lattice::AbstractLattice) =
        broadcast($op, lattice, number)
end
for op in (:*, :/, ://)
    @eval Base.$op(lattice::AbstractLattice, number::Number) =
        Lattice(($op)(lattice.data, number))
    @eval Base.$op(number::Number, lattice::AbstractLattice) = ($op)(lattice, number)
end

function Base.show(io::IO, x::AbstractLattice)
    if get(io, :compact, false) || get(io, :typeinfo, nothing) == typeof(x)
        Base.show_default(IOContext(io, :limit => true), x)  # From https://github.com/mauro3/Parameters.jl/blob/ecbf8df/src/Parameters.jl#L556
    else
        println(io, string(typeof(x)))
        for row in eachrow(x.data)
            println(io, ' ', join(row, "  "))
        end
    end
end
