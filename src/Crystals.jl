"""
# module Crystals



# Examples

```jldoctest
julia>
```
"""
module Crystals

using LinearAlgebra: det

using CoordinateTransformations
using StaticArrays

using Crystallography
using Crystallography.Directions: MetricTensor

export Cell, CellParameters, makelattice, cellvolume

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
CellParameters(bravais::BravaisLattice) = args -> CellParameters(bravais, args...)
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
makelattice(::BravaisLattice{Hexagonal,Primitive}, cell::CellParameters) =
    cell[1] * [
        1 0 0
        -1 / 2 √3 / 2 0
        0 0 cell[3] / cell[1]
    ]
function makelattice(
    ::BravaisLattice{Hexagonal,RhombohedralCentered},
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
    return _checkpositive(
        a * b * c * sqrt(1 - cos(α)^2 - cos(β)^2 - cos(γ)^2 + 2 * cos(α) * cos(β) * cos(γ)),
    )
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

StaticArrays.similar_type(::Type{<:CellParameters}, ::Type{T}, size::Size{(6,)}) where {T} =
    CellParameters{T}

end