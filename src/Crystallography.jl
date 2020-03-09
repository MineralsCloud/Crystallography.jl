module Crystallography

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
    BravaisLattice,
    BravaisLattice,
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
    RCenteredHexagonal
export pearsonsymbol, arithmeticclass, centeringof, crystalsystem

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
BaseCentering(T::Symbol) = T ∈ (:A, :B, :C) ? BaseCentering{T}() :
    throw(ArgumentError("centering must be either :A, :B, or :C!"))

const BravaisLattice = Tuple{CrystalSystem,Centering}
const PrimitiveTriclinic = Tuple{Triclinic,Primitive}
const PrimitiveMonoclinic = Tuple{Monoclinic,Primitive}
const ACenteredMonoclinic = Tuple{Monoclinic,BaseCentering{:A}}
const BCenteredMonoclinic = Tuple{Monoclinic,BaseCentering{:B}}
const CCenteredMonoclinic = Tuple{Monoclinic,BaseCentering{:C}}
const PrimitiveOrthorhombic = Tuple{Orthorhombic,Primitive}
const ACenteredOrthorhombic = Tuple{Orthorhombic,BaseCentering{:A}}
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
const RCenteredHexagonal = Tuple{Hexagonal,RhombohedralCentering}

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
pearsonsymbol(b::BravaisLattice) = pearsonsymbol(crystalsystem(b)) * pearsonsymbol(centeringof(b))

arithmeticclass(::Triclinic) = "-1"
arithmeticclass(::Monoclinic) = "2/m"
arithmeticclass(::Orthorhombic) = "mmm"
arithmeticclass(::Tetragonal) = "4/mmm"
arithmeticclass(::Hexagonal) = "6/mmm"
arithmeticclass(::Cubic) = "m-3m"
arithmeticclass(::Primitive) = "P"
arithmeticclass(::BaseCentering) = "S"
arithmeticclass(::BodyCentering) = "I"
arithmeticclass(::FaceCentering) = "F"
arithmeticclass(::RhombohedralCentering) = "R"
arithmeticclass(b::BravaisLattice) = arithmeticclass(crystalsystem(b)) * arithmeticclass(centeringof(b))
arithmeticclass(::RCenteredHexagonal) = "-3mR"

centeringof(b::BravaisLattice) = last(b)

crystalsystem(b::BravaisLattice) = first(b)

include("Crystals.jl")
include("Symmetry.jl")

end # module
