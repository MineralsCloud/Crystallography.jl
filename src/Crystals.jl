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

export Cell, CellParameters, transformation, transform

struct Cell{L<:AbstractVecOrMat,P<:AbstractVecOrMat,N<:AbstractVector,M<:Union{AbstractVector,Nothing}}
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
function CellParameters(a, b, c, α, β, γ)
    T = Base.promote_typeof(a, b, c, α, β, γ)
    return CellParameters{T}(a, b, c, α, β, γ)
end
CellParameters(::Triclinic, args...) = CellParameters(args...)
CellParameters(::Monoclinic, a, b, c, γ) = CellParameters(a, b, c, π / 2, π / 2, γ)
CellParameters(::Orthorhombic, a, b, c) = CellParameters(a, b, c, π / 2, π / 2, π / 2)
CellParameters(::Tetragonal, a, c) = CellParameters(a, a, c, π / 2, π / 2, π / 2)
CellParameters(::Cubic, a) = CellParameters(a, a, a, π / 2, π / 2, π / 2)
CellParameters(::Hexagonal, a, c) = CellParameters(a, a, c, π / 2, π / 2, 2π / 3)
CellParameters(::Trigonal, a, c) = CellParameters(a, a, c, π / 2, π / 2, 2π / 3)
CellParameters(::BravaisLattice{RhombohedralCentered,Hexagonal}, a, α) =
    CellParameters(a, a, a, α, α, α)

transformation(
    ::Type{BravaisLattice{Primitive,Cubic}},
    ::Type{RealSpace},
    cell::CellParameters,
) = cell.a * [
    1 0 0
    0 1 0
    0 0 1
] |> LinearMap
transformation(
    ::Type{BravaisLattice{BodyCentered,Cubic}},
    ::Type{RealSpace},
    cell::CellParameters,
) = cell.a / 2 * [
    1 1 1
    -1 1 1
    -1 -1 1
] |> LinearMap
transformation(
    ::Type{BravaisLattice{FaceCentered,Cubic}},
    ::Type{RealSpace},
    cell::CellParameters,
) = cell.a / 2 * [
    -1 0 1
    0 1 1
    -1 1 0
] |> LinearMap
transformation(
    ::Type{BravaisLattice{Primitive,Hexagonal}},
    ::Type{RealSpace},
    cell::CellParameters,
) = cell.a * [
    1 0 0
    -1 / 2 √3 / 2 0
    0 0 cell.c / cell.a
] |> LinearMap
function transformation(
    ::Type{BravaisLattice{RhombohedralCentered,Hexagonal}},
    ::Type{RealSpace},
    cell::CellParameters,
)
    r = cos(cell.α)
    tx = sqrt((1 - r) / 2)
    ty = sqrt((1 - r) / 6)
    tz = sqrt((1 + 2r) / 3)
    return cell.a * [
        tx -ty tz
        0 2ty tz
        -tx -ty tz
    ] |> LinearMap
end
transformation(
    ::Type{BravaisLattice{Primitive,Tetragonal}},
    ::Type{RealSpace},
    cell::CellParameters,
) = cell.a * [
    1 0 0
    0 1 0
    0 0 cell.c / cell.a
] |> LinearMap
function transformation(
    ::Type{BravaisLattice{BodyCentered,Tetragonal}},
    ::Type{RealSpace},
    cell::CellParameters,
)
    r = cell.c / cell.a
    return cell.a / 2 * [
        1 -1 r
        1 1 r
        -1 -1 r
    ] |> LinearMap
end
transformation(
    ::Type{BravaisLattice{Primitive,Orthorhombic}},
    ::Type{RealSpace},
    cell::CellParameters,
) = [
    cell.a 0 0
    0 cell.b 0
    0 0 cell.c
] |> LinearMap
# TODO: BravaisLattice{CCentered,Orthorhombic}
function transformation(
    ::Type{BravaisLattice{FaceCentered,Orthorhombic}},
    ::Type{RealSpace},
    cell::CellParameters,
)
    a, b, c = cell.a, cell.b, cell.c
    return [
        a 0 c
        a b 0
        0 b c
    ] / 2 |> LinearMap
end
function transformation(
    ::Type{BravaisLattice{BodyCentered,Orthorhombic}},
    ::Type{RealSpace},
    cell::CellParameters,
)
    a, b, c = cell.a, cell.b, cell.c
    return [
        a b c
        -a b c
        -a -b c
    ] / 2 |> LinearMap
end
function transformation(
    ::Type{BravaisLattice{Primitive,Monoclinic}},
    ::Type{RealSpace},
    cell::CellParameters,
)
    a, b, c = cell.a, cell.b, cell.c
    return [
        a 0 0
        0 b 0
        c * cos(cell.β) 0 c * sin(cell.β)
    ] |> LinearMap
end
# TODO: BravaisLattice{BCentered,Monoclinic}
# TODO: BravaisLattice{Primitive,Triclinic}
transformation(
    ::Type{T},
    ::Type{ReciprocalSpace},
    cell::CellParameters,
) where {T<:BravaisLattice} = inv(transformation(T, RealSpace, cell))

transform(m::LinearMap, a::CrystalCoordinates{T}) where {T} =
    CrystalCoordinates{T}(m.linear * collect(a))

StaticArrays.similar_type(::Type{<:CellParameters}, ::Type{T}, size::Size{(6,)}) where {T} =
    CellParameters{T}

end
