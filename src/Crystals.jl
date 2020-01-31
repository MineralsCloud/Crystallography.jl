"""
# module Crystals



# Examples

```jldoctest
julia>
```
"""
module Crystals

using CoordinateTransformations
using StaticArrays

using Crystallography

export Cell, CellParameters, makelattice

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
end
function CellParameters(a, b, c, α, β, γ, angle_iscosine::Bool = false)
    v = angle_iscosine ? (a, b, c, acos(α), acos(β), acos(γ)) : (a, b, c, α, β, γ)
    return CellParameters{Base.promote_typeof(v...)}(v...)
end
CellParameters(::BravaisLattice{Triclinic}, args...) = CellParameters(args...)  # Triclinic
CellParameters(::BravaisLattice{Monoclinic}, a, b, c, γ) =
    CellParameters(a, b, c, π / 2, π / 2, γ)
CellParameters(::BravaisLattice{Orthorhombic}, a, b, c) =
    CellParameters(a, b, c, π / 2, π / 2, π / 2)
CellParameters(::BravaisLattice{Tetragonal}, a, c) =
    CellParameters(a, a, c, π / 2, π / 2, π / 2)
CellParameters(::BravaisLattice{Cubic}, a) = CellParameters(a, a, a, π / 2, π / 2, π / 2)
CellParameters(::BravaisLattice{Hexagonal}, a, c) =
    CellParameters(a, a, c, π / 2, π / 2, 2π / 3)
CellParameters(::BravaisLattice{Trigonal}, a, c) =
    CellParameters(a, a, c, π / 2, π / 2, 2π / 3)
CellParameters(::BravaisLattice{Hexagonal,RhombohedralCentered}, a, α) =
    CellParameters(a, a, a, α, α, α)

makelattice(ibrav::Integer, cell::CellParameters) = makelattice(BravaisLattice(ibrav), cell)
makelattice(::BravaisLattice{Cubic,Primitive}, cell::CellParameters) =
    cell.a * [
        1 0 0
        0 1 0
        0 0 1
    ]
makelattice(::BravaisLattice{Cubic,BodyCentered}, cell::CellParameters) =
    cell.a / 2 * [
        1 1 1
        -1 1 1
        -1 -1 1
    ]
makelattice(::BravaisLattice{Cubic,FaceCentered}, cell::CellParameters) =
    cell.a / 2 * [
        -1 0 1
        0 1 1
        -1 1 0
    ]
makelattice(::BravaisLattice{Hexagonal,Primitive}, cell::CellParameters) =
    cell.a * [
        1 0 0
        -1 / 2 √3 / 2 0
        0 0 cell.c / cell.a
    ]
function makelattice(::BravaisLattice{Hexagonal,RhombohedralCentered}, cell::CellParameters)
    r = cos(cell.α)
    tx = sqrt((1 - r) / 2)
    ty = sqrt((1 - r) / 6)
    tz = sqrt((1 + 2r) / 3)
    return cell.a * [
        tx -ty tz
        0 2ty tz
        -tx -ty tz
    ]
end
makelattice(::BravaisLattice{Tetragonal,Primitive}, cell::CellParameters) =
    cell.a * [
        1 0 0
        0 1 0
        0 0 cell.c / cell.a
    ]
function makelattice(::BravaisLattice{Tetragonal,BodyCentered}, cell::CellParameters)
    r = cell.c / cell.a
    return cell.a / 2 * [
        1 -1 r
        1 1 r
        -1 -1 r
    ]
end
makelattice(::BravaisLattice{Orthorhombic,Primitive}, cell::CellParameters) = [
    cell.a 0 0
    0 cell.b 0
    0 0 cell.c
]
# TODO: BravaisLattice{Orthorhombic,CCentered}
function makelattice(::BravaisLattice{Orthorhombic,FaceCentered}, cell::CellParameters)
    a, b, c = cell.a, cell.b, cell.c
    return [
        a 0 c
        a b 0
        0 b c
    ] / 2
end
function makelattice(::BravaisLattice{Orthorhombic,BodyCentered}, cell::CellParameters)
    a, b, c = cell.a, cell.b, cell.c
    return [
        a b c
        -a b c
        -a -b c
    ] / 2
end
function makelattice(::BravaisLattice{Monoclinic,Primitive}, cell::CellParameters)
    a, b, c = cell.a, cell.b, cell.c
    return [
        a 0 0
        0 b 0
        c * cos(cell.β) 0 c * sin(cell.β)
    ]
end
# TODO: BravaisLattice{Monoclinic,BCentered}
# TODO: BravaisLattice{Triclinic,Primitive}

StaticArrays.similar_type(::Type{<:CellParameters}, ::Type{T}, size::Size{(6,)}) where {T} =
    CellParameters{T}

end
