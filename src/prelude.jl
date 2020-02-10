using StaticArrays: FieldVector, Size

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
    pearsonsymbol,
    arithmeticclass,
    centeringof,
    crystalsystem,
    dimensionof

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
BravaisLattice(::A, ::B, N::Integer) where {A,B} = BravaisLattice{A,B,supertype(A).parameters[1]}()
BravaisLattice(::A, ::B) where {A,B} = BravaisLattice{A,B,3}()
function BravaisLattice(ibrav::Integer)
    return if ibrav == 1
        BravaisLattice(Cubic(), Primitive())
    elseif ibrav == 2
        BravaisLattice(Cubic(), FaceCentered())
    elseif ibrav == 3
        BravaisLattice(Cubic(), BodyCentered())
    elseif ibrav == -3
        BravaisLattice(Cubic(), BodyCentered(), 2)
    elseif ibrav == 4
        BravaisLattice(Hexagonal(), Primitive())
    elseif ibrav == 5
        BravaisLattice(Hexagonal(), RhombohedralCentered())
    elseif ibrav == -5
        BravaisLattice(Hexagonal(), RhombohedralCentered(), 2)
    elseif ibrav == 6
        BravaisLattice(Tetragonal(), Primitive())
    elseif ibrav == 7
        BravaisLattice(Tetragonal(), BodyCentered())
    elseif ibrav == 8
        BravaisLattice(Orthorhombic(), Primitive())
    elseif ibrav == 9
        BravaisLattice(Orthorhombic(), BaseCentered(:B))
    elseif ibrav == -9
        BravaisLattice(Orthorhombic(), BaseCentered(:C))
    elseif ibrav == 91  # In QE 6.5
        BravaisLattice(Orthorhombic(), BaseCentered(:C))
    elseif ibrav == 10
        BravaisLattice(Orthorhombic(), FaceCentered())
    elseif ibrav == 11
        BravaisLattice(Orthorhombic(), BodyCentered())
    elseif ibrav == 12
        BravaisLattice(Monoclinic(), Primitive())
    elseif ibrav == -12
        BravaisLattice(Monoclinic(), Primitive())
    elseif ibrav == 13
        BravaisLattice(Monoclinic(), BaseCentered(:B))
    elseif ibrav == 14
        BravaisLattice(Triclinic(), Primitive())
    else
        error("undefined lattice!")
    end
end # function BravaisLattice
function BravaisLattice(; symbol::Bool = false)
    x = (
        BravaisLattice(Triclinic(), Primitive()),
        BravaisLattice(Monoclinic(), Primitive()),
        BravaisLattice(Monoclinic(), BaseCentered(:B)),
        BravaisLattice(Orthorhombic(), Primitive()),
        BravaisLattice(Orthorhombic(), BaseCentered(:C)),
        BravaisLattice(Orthorhombic(), BodyCentered()),
        BravaisLattice(Orthorhombic(), FaceCentered()),
        BravaisLattice(Tetragonal(), Primitive()),
        BravaisLattice(Tetragonal(), BodyCentered()),
        BravaisLattice(Cubic(), Primitive()),
        BravaisLattice(Cubic(), BodyCentered()),
        BravaisLattice(Cubic(), FaceCentered()),
        BravaisLattice(Hexagonal(), Primitive()),
        BravaisLattice(Hexagonal(), RhombohedralCentered()),
    )
    return symbol ? pearsonsymbol.(x) : x
end  # function BravaisLattice

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
arithmeticclass(::BravaisLattice{A,B}) where {A,B} = arithmeticclass(A()) * arithmeticclass(B())
arithmeticclass(::BravaisLattice{Hexagonal{3},RhombohedralCentered}) = "-3mR"
arithmeticclass(::BravaisLattice{Oblique}) = "2p"
arithmeticclass(::BravaisLattice{Rectangular,Primitive}) = "2mmp"
arithmeticclass(::BravaisLattice{Rectangular,<:BaseCentered}) = "2mmc"
arithmeticclass(::BravaisLattice{Square}) = "4mmp"
arithmeticclass(::BravaisLattice{Hexagonal{2}}) = "6mmh"

centeringof(::BravaisLattice{C,T}) where {C,T} = T()

crystalsystem(::BravaisLattice{C}) where {C} = C()

dimensionof(::BravaisLattice{C,T,N}) where {C,T,N} = N

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
