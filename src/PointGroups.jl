module PointGroups

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



end
