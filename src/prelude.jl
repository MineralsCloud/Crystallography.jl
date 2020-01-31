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

abstract type CrystalSystem end
struct Triclinic <: CrystalSystem end
struct Monoclinic <: CrystalSystem end
struct Orthorhombic <: CrystalSystem end
struct Tetragonal <: CrystalSystem end
struct Cubic <: CrystalSystem end
struct Trigonal <: CrystalSystem end
struct Hexagonal <: CrystalSystem end

abstract type CenteringType end
struct Primitive <: CenteringType end
struct BodyCentered <: CenteringType end
struct FaceCentered <: CenteringType end
struct RhombohedralCentered <: CenteringType end
struct BaseCentered{T} <: CenteringType end
function BaseCentered(T::Symbol)
    @assert T ∈ (:A, :B, :C)
    return BaseCentered{T}()
end # function BaseCentered

struct BravaisLattice{B<:CenteringType,C<:CrystalSystem} end
BravaisLattice(::B, ::C) where {B,C} = BravaisLattice{B,C}()
function BravaisLattice(ibrav::Integer)
    return if ibrav == 1
        BravaisLattice(Primitive(), Cubic())
    elseif ibrav == 2
        BravaisLattice(FaceCentered(), Cubic())
    elseif ibrav ∈ (3, -3)
        BravaisLattice(BodyCentered(), Cubic())
    elseif ibrav == 4
        BravaisLattice(Primitive(), Hexagonal())
    elseif ibrav ∈ (5, -5)
        BravaisLattice(RhombohedralCentered(), Hexagonal())
    elseif ibrav == 6
        BravaisLattice(Primitive(), Tetragonal())
    elseif ibrav == 7
        BravaisLattice(BodyCentered(), Tetragonal())
    elseif ibrav == 8
        BravaisLattice(Primitive(), Orthorhombic())
    elseif ibrav == 9
        BravaisLattice(BaseCentered(:A), Orthorhombic())
    elseif ibrav == -9
        BravaisLattice(BaseCentered(:C), Orthorhombic())
    elseif ibrav == 10
        BravaisLattice(FaceCentered(), Orthorhombic())
    elseif ibrav == 11
        BravaisLattice(BodyCentered(), Orthorhombic())
    elseif ibrav ∈ (12, -12)
        BravaisLattice(Primitive(), Monoclinic())
    elseif ibrav == 13
        BravaisLattice(BaseCentered(:A), Monoclinic())
    elseif ibrav == 14
        BravaisLattice(Primitive(), Triclinic())
    else
        error("undefined lattice!")
    end
end # function BravaisLattice

nomenclature(::Triclinic) = "a"
nomenclature(::Monoclinic) = "m"
nomenclature(::Orthorhombic) = "o"
nomenclature(::Tetragonal) = "t"
nomenclature(::Cubic) = "c"
nomenclature(::Hexagonal) = "h"
nomenclature(::Trigonal) = "h"
nomenclature(::Primitive) = "P"
nomenclature(::BaseCentered{T}) where {T} = string(T)
nomenclature(::BodyCentered) = "I"
nomenclature(::FaceCentered) = "F"
nomenclature(::RhombohedralCentered) = "R"
nomenclature(::BravaisLattice{B,C}) where {B,C} = nomenclature(C()) * nomenclature(B())

function allbravaislattices(; symbol::Bool = false)
    x = (
        BravaisLattice(Primitive(), Triclinic()),
        BravaisLattice(Primitive(), Monoclinic()),
        BravaisLattice(BaseCentered(:B), Monoclinic()),
        BravaisLattice(Primitive(), Orthorhombic()),
        BravaisLattice(BaseCentered(:C), Orthorhombic()),
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
