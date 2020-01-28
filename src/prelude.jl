using StaticArrays: FieldVector, Size

import StaticArrays

export AbstractSpace,
    RealSpace,
    ReciprocalSpace,
    CrystalCoordinates,
    CartesianCoordinates,
    CrystalSystem,
    Triclinic,
    Monoclinic,
    Orthorhombic,
    Tetragonal,
    Cubic,
    Trigonal,
    Hexagonal,
    CenteringType,
    BaseCentered,
    Primitive,
    ACentered,
    BCentered,
    CCentered,
    BodyCentered,
    FaceCentered,
    RhombohedralCentered,
    BravaisLattice,
    nomenclature,
    allbravaislattices,
    centeringtype,
    crystalsystem

abstract type AbstractSpace end
struct RealSpace <: AbstractSpace end
struct ReciprocalSpace <: AbstractSpace end

abstract type AbstractCoordinates{T} <: FieldVector{3,T} end

struct CrystalCoordinates{S<:AbstractSpace,T} <: AbstractCoordinates{T}
    x::T
    y::T
    z::T
end
CrystalCoordinates{S}(x::T, y::T, z::T) where {S,T} = CrystalCoordinates{S,T}(x, y, z)

struct CartesianCoordinates{T} <: AbstractCoordinates{T}
    x::T
    y::T
    z::T
end

StaticArrays.similar_type(
    ::Type{<:CrystalCoordinates{S}},
    ::Type{T},
    size::Size{(3,)},
) where {S,T} = CrystalCoordinates{S,T}
StaticArrays.similar_type(
    ::Type{<:CartesianCoordinates},
    ::Type{T},
    size::Size{(3,)},
) where {T} = CartesianCoordinates{T}

Base.inv(::RealSpace) = ReciprocalSpace()
Base.inv(::ReciprocalSpace) = RealSpace()

abstract type CrystalSystem end
struct Triclinic <: CrystalSystem end
struct Monoclinic <: CrystalSystem end
struct Orthorhombic <: CrystalSystem end
struct Tetragonal <: CrystalSystem end
struct Cubic <: CrystalSystem end
struct Trigonal <: CrystalSystem end
struct Hexagonal <: CrystalSystem end

abstract type CenteringType end
abstract type BaseCentered <: CenteringType end
struct Primitive <: CenteringType end
struct ACentered <: BaseCentered end
struct BCentered <: BaseCentered end
struct CCentered <: BaseCentered end
struct BodyCentered <: CenteringType end
struct FaceCentered <: CenteringType end
struct RhombohedralCentered <: CenteringType end

struct BravaisLattice{B<:CenteringType,C<:CrystalSystem} end
BravaisLattice(::B, ::C) where {B,C} = BravaisLattice{B,C}()

nomenclature(::Triclinic) = "a"
nomenclature(::Monoclinic) = "m"
nomenclature(::Orthorhombic) = "o"
nomenclature(::Tetragonal) = "t"
nomenclature(::Cubic) = "c"
nomenclature(::Hexagonal) = "h"
nomenclature(::Trigonal) = "h"
nomenclature(::Primitive) = "P"
nomenclature(::ACentered) = "A"
nomenclature(::BCentered) = "B"
nomenclature(::CCentered) = "C"
nomenclature(::BodyCentered) = "I"
nomenclature(::FaceCentered) = "F"
nomenclature(::RhombohedralCentered) = "R"
nomenclature(::BravaisLattice{B,C}) where {B,C} = nomenclature(C()) * nomenclature(B())

function allbravaislattices(; symbol::Bool = false)
    x = (
        BravaisLattice(Primitive(), Triclinic()),
        BravaisLattice(Primitive(), Monoclinic()),
        BravaisLattice(BCentered(), Monoclinic()),
        BravaisLattice(Primitive(), Orthorhombic()),
        BravaisLattice(CCentered(), Orthorhombic()),
        BravaisLattice(BodyCentered(), Orthorhombic()),
        BravaisLattice(FaceCentered(), Orthorhombic()),
        BravaisLattice(Primitive(), Tetragonal()),
        BravaisLattice(BodyCentered(), Tetragonal()),
        BravaisLattice(Primitive(), Cubic()),
        BravaisLattice(BodyCentered(), Cubic()),
        BravaisLattice(FaceCentered(), Cubic()),
        BravaisLattice(Primitive(), Hexagonal()),
        BravaisLattice(RhombohedralCentered(), Hexagonal()),
    )
    return symbol ? map(nomenclature, x) : x
end  # function allbravaislattices

centeringtype(::BravaisLattice{C}) where {C} = C()

crystalsystem(::BravaisLattice{C,T}) where {C,T} = T()

Base.show(io::IO, t::CrystalSystem) = show(io, lowercase(string(t)))
Base.show(io::IO, t::CenteringType) = show(io, lowercase(string(t)))
Base.show(io::IO, t::BravaisLattice{C,T}) where {C,T} =
    show(io, lowercase(string(C, ' ', T)))
