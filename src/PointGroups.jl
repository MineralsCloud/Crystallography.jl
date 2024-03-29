module PointGroups

using MLStyle: @match

using CrystallographyBase: CrystalSystem

export PointGroup
export pointgroups,
    hermann_mauguin, international, schönflies, schoenflies, orderof, laueclasses

abstract type PointGroup end
abstract type CyclicGroup <: PointGroup end
struct C1 <: CyclicGroup end
struct C2 <: CyclicGroup end
struct C3 <: CyclicGroup end
struct C4 <: CyclicGroup end
struct C6 <: CyclicGroup end
struct C2v <: CyclicGroup end
struct C3v <: CyclicGroup end
struct C4v <: CyclicGroup end
struct C6v <: CyclicGroup end
struct C1h <: CyclicGroup end
struct C2h <: CyclicGroup end
struct C3h <: CyclicGroup end
struct C4h <: CyclicGroup end
struct C6h <: CyclicGroup end
abstract type SpiegelGroup <: PointGroup end
struct S2 <: SpiegelGroup end
struct S4 <: SpiegelGroup end
struct S6 <: SpiegelGroup end
abstract type DihedralGroup <: PointGroup end
struct D2 <: DihedralGroup end
struct D3 <: DihedralGroup end
struct D4 <: DihedralGroup end
struct D6 <: DihedralGroup end
struct D2h <: DihedralGroup end
struct D3h <: DihedralGroup end
struct D4h <: DihedralGroup end
struct D6h <: DihedralGroup end
struct D2d <: DihedralGroup end
struct D3d <: DihedralGroup end
abstract type TetrahedralGroup <: PointGroup end
struct T <: TetrahedralGroup end
struct Th <: TetrahedralGroup end
struct Td <: TetrahedralGroup end
abstract type OctahedralGroup <: PointGroup end
struct O <: OctahedralGroup end
struct Oh <: OctahedralGroup end
const C1v = C1h
const Cs = C1h
const Ci = S2
const C3i = S6
const D1 = C2
const V = D2
const D1h = C2v
const D1d = C2h
const V2 = D2d
const Vh = D2h
function PointGroup(str::AbstractString)
    str = lowercase(filter(!isspace, str))
    @match str begin
        # International notation
        "1" => C1()
        "1̄" => S2()
        "2" => C2()
        "m" => C1h()
        "3" => C3()
        "4" => C4()
        "4̄" => S4()
        "2/m" => C2h()
        "222" => D2()
        "mm2" => C2v()
        "3̄" => S6()
        "6" => C6()
        "6̄" => C3h()
        "32" => D3()
        "3m" => C3v()
        "mmm" => D2h()
        "4/m" => C4h()
        "422" => D4()
        "4mm" => C4v()
        "4̄2m" => D2d()
        "6/m" => C6h()
        "23" => T()
        "3̄m" => D3d()
        "622" => D6()
        "6mm" => C6v()
        "6̄m2" => D3h()
        "4/mmm" => D4h()
        "6/mmm" => D6h()
        "m3̄" => Th()
        "432" => O()
        "4̄3m" => Td()
        "m3̄m" => Oh()
        # Schoenflies notation
        "c1" => C1()
        "c2" => C2()
        "c3" => C3()
        "c4" => C4()
        "c6" => C6()
        "c2v" => C2v()
        "c3v" => C3v()
        "c4v" => C4v()
        "c6v" => C6v()
        "c1h" => C1h()
        "c2h" => C2h()
        "c3h" => C3h()
        "c4h" => C4h()
        "c6h" => C6h()
        "s2" => S2()
        "s4" => S4()
        "s6" => S6()
        "d2" => D2()
        "d3" => D3()
        "d4" => D4()
        "d6" => D6()
        "d2h" => D2h()
        "d3h" => D3h()
        "d4h" => D4h()
        "d6h" => D6h()
        "d2d" => D2d()
        "d3d" => D3d()
        "t" => T()
        "th" => Th()
        "td" => Td()
        "o" => O()
        "oh" => Oh()
        "c1v" => C1h()
        "cs" => C1h()
        "ci" => S2()
        "c3i" => S6()
        "d1" => C2()
        "v" => D2()
        "d1h" => C2v()
        "d1d" => C2h()
        "v2" => D2d()
        "vh" => D2h()
        _ => throw(ArgumentError("\"$str\" is not a valid point group symbol!"))
    end
end

struct LaueClass{T<:PointGroup} end

function laueclasses(system::CrystalSystem)
    @match system begin
        CrystalSystem(:Triclinic) => (LaueClass{Ci}(),)
        CrystalSystem(:Monoclinic) => (LaueClass{C2h}(),)
        CrystalSystem(:Orthorhombic) => (LaueClass{D2h}(),)
        CrystalSystem(:Tetragonal) => (LaueClass{C4h}(), LaueClass{D4h}())
        CrystalSystem(:Trigonal) => (LaueClass{C3i}(), LaueClass{D3d}())
        CrystalSystem(:Hexagonal) => (LaueClass{C6h}(), LaueClass{D6h}())
        CrystalSystem(:Cubic) => (LaueClass{Th}(), LaueClass{Oh}())
        _ => error("Invalid crystal type")
    end
end

function pointgroups(crystal)
    @match crystal begin
        CrystalSystem(:Triclinic) => (C1(), S2())
        CrystalSystem(:Monoclinic) => (C2(), C1h(), C2h())
        CrystalSystem(:Orthorhombic) => (D2(), C2v(), D2h())
        CrystalSystem(:Tetragonal) => (C4(), S4(), C4h(), D4(), C4v(), D2d(), D4h())
        CrystalSystem(:Trigonal) => (C3(), S6(), D3(), C3v(), D3d())
        CrystalSystem(:Hexagonal) => (C6(), C3h(), C6h(), D6(), C6v(), D3h(), D6h())
        CrystalSystem(:Cubic) => (T(), Th(), O(), Td(), Oh())
        _ => error("Invalid crystal type")
    end
end

hermann_mauguin(::C1) = "1"
hermann_mauguin(::S2) = "1̄"
hermann_mauguin(::C2) = "2"
hermann_mauguin(::C1h) = "m"
hermann_mauguin(::C3) = "3"
hermann_mauguin(::C4) = "4"
hermann_mauguin(::S4) = "4̄"
hermann_mauguin(::C2h) = "2/m"
hermann_mauguin(::D2) = "222"
hermann_mauguin(::C2v) = "mm2"
hermann_mauguin(::S6) = "3̄"
hermann_mauguin(::C6) = "6"
hermann_mauguin(::C3h) = "6̄"
hermann_mauguin(::D3) = "32"
hermann_mauguin(::C3v) = "3m"
hermann_mauguin(::D2h) = "mmm"
hermann_mauguin(::C4h) = "4/m"
hermann_mauguin(::D4) = "422"
hermann_mauguin(::C4v) = "4mm"
hermann_mauguin(::D2d) = "4̄2m"
hermann_mauguin(::C6h) = "6/m"
hermann_mauguin(::T) = "23"
hermann_mauguin(::D3d) = "3̄m"
hermann_mauguin(::D6) = "622"
hermann_mauguin(::C6v) = "6mm"
hermann_mauguin(::D3h) = "6̄m2"
hermann_mauguin(::D4h) = "4/mmm"
hermann_mauguin(::D6h) = "6/mmm"
hermann_mauguin(::Th) = "m3̄"
hermann_mauguin(::O) = "432"
hermann_mauguin(::Td) = "4̄3m"
hermann_mauguin(::Oh) = "m3̄m"
const international = hermann_mauguin

function schönflies(g::PointGroup)
    if isconcretetype(typeof(g))
        return string(nameof(typeof(g)))
    else
        throw(ArgumentError("$g is not a specific point group!"))
    end
end
const schoenflies = schönflies

orderof(::C1) = 1
orderof(::S2) = 2
orderof(::C2) = 2
orderof(::C1h) = 2
orderof(::C3) = 3
orderof(::C4) = 4
orderof(::S4) = 4
orderof(::C2h) = 4
orderof(::D2) = 4
orderof(::C2v) = 4
orderof(::S6) = 6
orderof(::C6) = 6
orderof(::C3h) = 6
orderof(::D3) = 6
orderof(::C3v) = 6
orderof(::D2h) = 8
orderof(::C4h) = 8
orderof(::D4) = 8
orderof(::C4v) = 8
orderof(::D2d) = 8
orderof(::C6h) = 12
orderof(::T) = 12
orderof(::D3d) = 12
orderof(::D6) = 12
orderof(::C6v) = 12
orderof(::D3h) = 12
orderof(::D4h) = 16
orderof(::D6h) = 24
orderof(::Th) = 24
orderof(::O) = 24
orderof(::Td) = 24
orderof(::Oh) = 48

end
