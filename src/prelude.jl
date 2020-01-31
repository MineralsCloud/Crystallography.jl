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
    Centering,
    BaseCentered,
    Primitive,
    BodyCentered,
    FaceCentered,
    RhombohedralCentered,
    BaseCentered,
    BravaisLattice,
    pearsonsymbol,
    centeringof,
    crystalsystem,
    viewsetting

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
struct BodyCentered <: Centering end
struct FaceCentered <: Centering end
struct RhombohedralCentered <: Centering end
struct BaseCentered{T} <: Centering end
function BaseCentered(T::Symbol)
    @assert T ∈ (:A, :B, :C)
    return BaseCentered{T}()
end # function BaseCentered

struct BravaisLattice{A<:CrystalSystem,B<:Centering,N} end
BravaisLattice(::A, ::B, N::Integer) where {A,B} = BravaisLattice{A,B,N}()
BravaisLattice(::A, ::B) where {A,B} = BravaisLattice{A,B,1}()
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

centeringof(::BravaisLattice{C,T}) where {C,T} = T()

crystalsystem(::BravaisLattice{C}) where {C} = C()

viewsetting(::BravaisLattice{C,T,N}) where {C,T,N} = N

Base.show(io::IO, t::CrystalSystem) = show(io, lowercase(string(t)))
Base.show(io::IO, t::Centering) = show(io, lowercase(string(t)))
Base.show(io::IO, t::BravaisLattice{C,T}) where {C,T} =
    show(io, lowercase(string(T, ' ', C)))

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
