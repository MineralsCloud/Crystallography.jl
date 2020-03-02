module Crystallography

export CrystalSystem,
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
    BaseCentering,
    Primitive,
    BodyCentering,
    FaceCentering,
    RhombohedralCentering,
    BaseCentering,
    BravaisLattice,
    BravaisLattice,
    PrimitiveTriclinic,
    PrimitiveMonoclinic,
    BCenteredMonoclinic,
    CCenteredMonoclinic,
    PrimitiveOrthorhombic,
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
    RCenteredHexagonal
export pearsonsymbol,
    arithmeticclass,
    centeringof,
    dimensionof

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
Hexagonal(N::Int = 3) =
    N ∈ (2, 3) ? Hexagonal{N}() : throw(ArgumentError("hexagonal must be 2D or 3D!"))

abstract type Centering end
struct Primitive <: Centering end
struct BodyCentering <: Centering end
struct FaceCentering <: Centering end
struct RhombohedralCentering <: Centering end
struct BaseCentering{T} <: Centering end
BaseCentering(T::Symbol) = T ∈ (:A, :B, :C) ? BaseCentering{T}() :
    throw(ArgumentError("centering must be either :A, :B, or :C!"))

const BravaisLattice = Tuple{CrystalSystem,Centering}
const PrimitiveTriclinic = Tuple{Triclinic,Primitive}
const PrimitiveMonoclinic = Tuple{Monoclinic,Primitive}
const BCenteredMonoclinic = Tuple{Monoclinic,BaseCentering{:B}}
const CCenteredMonoclinic = Tuple{Monoclinic,BaseCentering{:C}}
const PrimitiveOrthorhombic = Tuple{Orthorhombic,Primitive}
const BCenteredOrthorhombic = Tuple{Orthorhombic,BaseCentering{:B}}
const CCenteredOrthorhombic = Tuple{Orthorhombic,BaseCentering{:C}}
const BodyCenteredOrthorhombic = Tuple{Orthorhombic,BodyCentering}
const FaceCenteredOrthorhombic = Tuple{Orthorhombic,FaceCentering}
const PrimitiveTetragonal = Tuple{Tetragonal,Primitive}
const BodyCenteredTetragonal = Tuple{Tetragonal,BodyCentering}
const PrimitiveCubic = Tuple{Cubic,Primitive}
const BodyCenteredCubic = Tuple{Cubic,BodyCentering}
const FaceCenteredCubic = Tuple{Cubic,FaceCentering}
const PrimitiveHexagonal = Tuple{Hexagonal,Primitive}
const RCenteredHexagonal = Tuple{Hexagonal{3},RhombohedralCentering}
const PrimitiveOblique = Tuple{Oblique,Primitive}
const PrimitiveRectangular = Tuple{Rectangular,Primitive}
const BaseCenteredRectangular = Tuple{Rectangular,BaseCentering}
const PrimitiveSquare = Tuple{Square,Primitive}
const PrimitiveHexagonal2D = Tuple{Hexagonal{2},Primitive}

const TetragonalBravais = Union{PrimitiveTetragonal,BodyCenteredTetragonal}
const CubicBravais = Union{PrimitiveCubic,BodyCenteredCubic,FaceCenteredCubic}
const OrthorhombicBravais = Union{
    PrimitiveOrthorhombic,
    BCenteredOrthorhombic,
    CCenteredOrthorhombic,
    BodyCenteredOrthorhombic,
    FaceCenteredCubic,
}
const MonoclinicBravais = Union{PrimitiveMonoclinic,BCenteredMonoclinic,CCenteredMonoclinic}

pearsonsymbol(::Triclinic) = "a"
pearsonsymbol(::Monoclinic) = "m"
pearsonsymbol(::Orthorhombic) = "o"
pearsonsymbol(::Tetragonal) = "t"
pearsonsymbol(::Cubic) = "c"
pearsonsymbol(::Hexagonal) = "h"
pearsonsymbol(::Trigonal) = "h"
pearsonsymbol(::Primitive) = "P"
pearsonsymbol(::BaseCentering{T}) where {T} = string(T)
pearsonsymbol(::BodyCentering) = "I"
pearsonsymbol(::FaceCentering) = "F"
pearsonsymbol(::RhombohedralCentering) = "R"
pearsonsymbol(b::BravaisLattice) = pearsonsymbol(first(b)) * pearsonsymbol(last(b))

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
arithmeticclass(::BaseCentering) = "S"
arithmeticclass(::BodyCentering) = "I"
arithmeticclass(::FaceCentering) = "F"
arithmeticclass(::RhombohedralCentering) = "R"
arithmeticclass(b::BravaisLattice) = arithmeticclass(first(b)) * arithmeticclass(last(b))
arithmeticclass(::RCenteredHexagonal) = "-3mR"
arithmeticclass(::PrimitiveOblique) = "2p"
arithmeticclass(::PrimitiveRectangular) = "2mmp"
arithmeticclass(::BaseCenteredRectangular) = "2mmc"
arithmeticclass(::PrimitiveSquare) = "4mmp"
arithmeticclass(::PrimitiveHexagonal2D) = "6mmh"

centeringof(b::BravaisLattice) = last(b)

dimensionof(::CrystalSystem{N}) where {N} = N
dimensionof(b::BravaisLattice) = dimensionof(first(b))

include("Crystals.jl")
include("Symmetry.jl")

end # module
