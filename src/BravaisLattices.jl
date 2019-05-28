"""
# module BravaisLattices



# Examples

```jldoctest
julia>
```
"""
module BravaisLattices

export CrystalSystem,
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
    allbravaislattices

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

abstract type BravaisLattice{B <: CenteringType,C <: CrystalSystem} end

nomenclature(::Type{Triclinic}) = "a"
nomenclature(::Type{Monoclinic}) = "m"
nomenclature(::Type{Orthorhombic}) = "o"
nomenclature(::Type{Tetragonal}) = "t"
nomenclature(::Type{Cubic}) = "c"
nomenclature(::Type{Hexagonal}) = "h"
nomenclature(::Type{Trigonal}) = "h"
nomenclature(::Type{Primitive}) = "P"
nomenclature(::Type{ACentered}) = "A"
nomenclature(::Type{BCentered}) = "B"
nomenclature(::Type{CCentered}) = "C"
nomenclature(::Type{BodyCentered}) = "I"
nomenclature(::Type{FaceCentered}) = "F"
nomenclature(::Type{RhombohedralCentered}) = "R"
nomenclature(::Type{BravaisLattice{B,C}}) where {B,C} = nomenclature(C) * nomenclature(B)

function allbravaislattices(; if_nomenclature = false)
    x = (BravaisLattice{Primitive,Triclinic}, BravaisLattice{Primitive,Monoclinic}, BravaisLattice{BCentered,Monoclinic},
    BravaisLattice{Primitive,Orthorhombic}, BravaisLattice{CCentered,Orthorhombic}, BravaisLattice{BodyCentered,Orthorhombic}, BravaisLattice{FaceCentered,Orthorhombic},
    BravaisLattice{Primitive,Tetragonal}, BravaisLattice{BodyCentered,Tetragonal},
    BravaisLattice{Primitive,Cubic}, BravaisLattice{BodyCentered,Cubic}, BravaisLattice{FaceCentered,Cubic},
    BravaisLattice{Primitive,Hexagonal}, BravaisLattice{RhombohedralCentered,Hexagonal})
    if_nomenclature ? x : map(nomenclature, x)
end  # function allbravaislattices

end
