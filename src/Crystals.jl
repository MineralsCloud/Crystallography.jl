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

export Cell, CellParameters, makelattice, transform

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
CellParameters(::BravaisLattice{C,Triclinic}, args...) where {C} = CellParameters(args...)  # Triclinic
CellParameters(::BravaisLattice{C,Monoclinic}, a, b, c, γ) where {C} =
    CellParameters(a, b, c, π / 2, π / 2, γ)
CellParameters(::BravaisLattice{C,Orthorhombic}, a, b, c) where {C} =
    CellParameters(a, b, c, π / 2, π / 2, π / 2)
CellParameters(::BravaisLattice{C,Tetragonal}, a, c) where {C} =
    CellParameters(a, a, c, π / 2, π / 2, π / 2)
CellParameters(::BravaisLattice{C,Cubic}, a) where {C} =
    CellParameters(a, a, a, π / 2, π / 2, π / 2)
CellParameters(::BravaisLattice{C,Hexagonal}, a, c) where {C} =
    CellParameters(a, a, c, π / 2, π / 2, 2π / 3)
CellParameters(::BravaisLattice{C,Trigonal}, a, c) where {C} =
    CellParameters(a, a, c, π / 2, π / 2, 2π / 3)
CellParameters(::BravaisLattice{RhombohedralCentered,Hexagonal}, a, α) =
    CellParameters(a, a, a, α, α, α)

makelattice(::BravaisLattice{Primitive,Cubic}, cell::CellParameters) =
    cell.a * [
        1 0 0
        0 1 0
        0 0 1
    ]
makelattice(::BravaisLattice{BodyCentered,Cubic}, cell::CellParameters) =
    cell.a / 2 * [
        1 1 1
        -1 1 1
        -1 -1 1
    ]
makelattice(::BravaisLattice{FaceCentered,Cubic}, cell::CellParameters) =
    cell.a / 2 * [
        -1 0 1
        0 1 1
        -1 1 0
    ]
makelattice(::BravaisLattice{Primitive,Hexagonal}, cell::CellParameters) =
    cell.a * [
        1 0 0
        -1 / 2 √3 / 2 0
        0 0 cell.c / cell.a
    ]
function makelattice(::BravaisLattice{RhombohedralCentered,Hexagonal}, cell::CellParameters)
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
makelattice(::BravaisLattice{Primitive,Tetragonal}, cell::CellParameters) =
    cell.a * [
        1 0 0
        0 1 0
        0 0 cell.c / cell.a
    ]
function makelattice(::BravaisLattice{BodyCentered,Tetragonal}, cell::CellParameters)
    r = cell.c / cell.a
    return cell.a / 2 * [
        1 -1 r
        1 1 r
        -1 -1 r
    ]
end
makelattice(::BravaisLattice{Primitive,Orthorhombic}, cell::CellParameters) = [
    cell.a 0 0
    0 cell.b 0
    0 0 cell.c
]
# TODO: BravaisLattice{CCentered,Orthorhombic}
function makelattice(::BravaisLattice{FaceCentered,Orthorhombic}, cell::CellParameters)
    a, b, c = cell.a, cell.b, cell.c
    return [
        a 0 c
        a b 0
        0 b c
    ] / 2
end
function makelattice(::BravaisLattice{BodyCentered,Orthorhombic}, cell::CellParameters)
    a, b, c = cell.a, cell.b, cell.c
    return [
        a b c
        -a b c
        -a -b c
    ] / 2
end
function makelattice(::BravaisLattice{Primitive,Monoclinic}, cell::CellParameters)
    a, b, c = cell.a, cell.b, cell.c
    return [
        a 0 0
        0 b 0
        c * cos(cell.β) 0 c * sin(cell.β)
    ]
end
# TODO: BravaisLattice{BCentered,Monoclinic}
# TODO: BravaisLattice{Primitive,Triclinic}

transform(m::LinearMap, a::CrystalCoordinates) = CrystalCoordinates(m.linear * collect(a))

StaticArrays.similar_type(::Type{<:CellParameters}, ::Type{T}, size::Size{(6,)}) where {T} =
    CellParameters{T}

end
