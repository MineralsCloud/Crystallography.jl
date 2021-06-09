using LinearAlgebra: Diagonal
using Spglib: Cell

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

struct Lattice{T}
    data::SMatrix{3,3,T,9}
end
Lattice(mat::AbstractMatrix) = Lattice(SMatrix{3,3}(mat))
Lattice(𝐚::AbstractVector, 𝐛::AbstractVector, 𝐜::AbstractVector) = Lattice(hcat(𝐚, 𝐛, 𝐜))
Lattice(lattice::Lattice) = lattice
Lattice(cell::Cell) = Lattice(cell.lattice)
function Lattice(a, b, c, α, β, γ)
    # From https://github.com/LaurentRDC/crystals/blob/dbb3a92/crystals/lattice.py#L321-L354
    v = cellvolume(1, 1, 1, α, β, γ)
    # reciprocal lattice
    a_recip = sind(α) / (a * v)
    csg = (cosd(α) * cosd(β) - cosd(γ)) / (sind(α) * sind(β))
    sg = sqrt(1 - csg^2)
    a1 = [1 / a_recip, -csg / sg / a_recip, cosd(β) * a]
    a2 = [0, b * sind(α), b * cosd(α)]
    a3 = [0, 0, c]
    return Lattice(a1, a2, a3)
end
@functor Lattice

basis_vectors(lattice::Lattice) = lattice[:, 1], lattice[:, 2], lattice[:, 3]

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
crystalsystem(lattice::Lattice) = crystalsystem(cellparameters(lattice))

function cellparameters(lattice::Lattice)
    𝐚, 𝐛, 𝐜 = basis_vectors(lattice)
    a, b, c = norm(𝐚), norm(𝐛), norm(𝐜)
    γ, β, α =
        acosd(dot(𝐚, 𝐛) / (a * b)), acosd(dot(𝐚, 𝐜) / (a * c)), acosd(dot(𝐛, 𝐜) / (b * c))
    return a, b, c, α, β, γ
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
function supercell(cell::Cell, expansion) end

Base.iterate(lattice::Lattice) = iterate(lattice.data)
Base.iterate(lattice::Lattice, state) = iterate(lattice.data, state)

Base.eltype(::Lattice{T}) where {T} = T

Base.length(::Lattice) = 9

Base.size(::Lattice) = (3, 3)
Base.size(::Lattice, dim::Integer) = dim <= 2 ? 3 : 1

Base.IteratorSize(::Type{<:Lattice}) = Base.HasShape{2}()

Base.axes(lattice::Lattice, dim::Integer) = axes(lattice.data, dim)

Base.getindex(lattice::Lattice, i) = getindex(lattice.data, i)
Base.getindex(lattice::Lattice, I::Vararg) = getindex(lattice.data, I...)

Base.firstindex(::Lattice) = 1

Base.lastindex(::Lattice) = 9

for op in (:+, :-)
    @eval Base.broadcast(::typeof($op), lattice::Lattice, number::Number) =
        Lattice(broadcast($op, lattice.data, number))
    @eval Base.broadcast(::typeof($op), number::Number, lattice::Lattice) =
        broadcast($op, lattice, number)
end
for op in (:*, :/, ://)
    @eval Base.$op(lattice::Lattice, number::Number) = Lattice(($op)(lattice.data, number))
    @eval Base.$op(number::Number, lattice::Lattice) = ($op)(lattice, number)
end
