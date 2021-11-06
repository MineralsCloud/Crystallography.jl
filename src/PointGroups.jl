module PointGroups

using ..Crystallography: Triclinic,
    Monoclinic,
    Orthorhombic,
    Tetragonal,
    Cubic,
    Trigonal,
    Hexagonal

export pointgroups, hermann_mauguin, international

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

pointgroups(::Triclinic) = (C1(), S2())
pointgroups(::Monoclinic) = (C2(), C1h(), C2h())
pointgroups(::Orthorhombic) = (D2(), C2v(), D2h())
pointgroups(::Tetragonal) = (C4(), S4(), C4h(), D4(), C4v(), D2d(), D4h())
pointgroups(::Trigonal) = (C3(), S6(), D3(), C3v(), D3d())
pointgroups(::Hexagonal) = (C6(), C3h(), C6h(), D6(), C6v(), D3h(), D6h())
pointgroups(::Cubic) = (T(), Th(), O(), Td(), Oh())

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

end
