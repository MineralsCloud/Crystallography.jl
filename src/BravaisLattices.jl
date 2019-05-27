"""
# module BravaisLattices



# Examples

```jldoctest
julia>
```
"""
module BravaisLattices

export CrystalSystem,
    TriclinicSystem,
    MonoclinicSystem,
    OrthorhombicSystem,
    TetragonalSystem,
    CubicSystem,
    TrigonalSystem,
    HexagonalSystem,
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
    nomenclature

abstract type CrystalSystem end
struct TriclinicSystem <: CrystalSystem end
struct MonoclinicSystem <: CrystalSystem end
struct OrthorhombicSystem <: CrystalSystem end
struct TetragonalSystem <: CrystalSystem end
struct CubicSystem <: CrystalSystem end
struct TrigonalSystem <: CrystalSystem end
struct HexagonalSystem <: CrystalSystem end

abstract type CenteringType end
abstract type BaseCentered <: CenteringType end
struct Primitive <: CenteringType end
struct ACentered <: BaseCentered end
struct BCentered <: BaseCentered end
struct CCentered <: BaseCentered end
struct BodyCentered <: CenteringType end
struct FaceCentered <: CenteringType end
struct RhombohedralCentered <: CenteringType end

abstract type BravaisLattice{B,C} end

nomenclature(::Type{TriclinicSystem}) = "a"
nomenclature(::Type{MonoclinicSystem}) = "m"
nomenclature(::Type{OrthorhombicSystem}) = "o"
nomenclature(::Type{TetragonalSystem}) = "t"
nomenclature(::Type{CubicSystem}) = "c"
nomenclature(::Type{HexagonalSystem}) = "h"
nomenclature(::Type{TrigonalSystem}) = "h"
nomenclature(::Type{Primitive}) = "P"
nomenclature(::Type{ACentered}) = "A"
nomenclature(::Type{BCentered}) = "B"
nomenclature(::Type{CCentered}) = "C"
nomenclature(::Type{BodyCentered}) = "I"
nomenclature(::Type{FaceCentered}) = "F"
nomenclature(::Type{RhombohedralCentered}) = "R"
nomenclature(::Type{BravaisLattice{B,C}}) where {B <: CenteringType, C <: CrystalSystem} = nomenclature(C) * nomenclature(B)

end
