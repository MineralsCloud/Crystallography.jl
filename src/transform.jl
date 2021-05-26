# Idea from https://spglib.github.io/spglib/definition.html#transformation-to-the-primitive-cell
struct PrimitiveFromStandardized{T<:Centering}
    transformation::SMatrix{3,3}
end
struct StandardizedFromPrimitive{T<:Centering}
    transformation::SMatrix{3,3}
end

(::PrimitiveFromStandardized{Primitive})(x) = x
(x::PrimitiveFromStandardized)(lattice::Lattice) = Lattice(lattice.data * x.transformation)
(x::PrimitiveFromStandardized)(coord::AbstractVector) = inv(x.transformation) * coord

(::StandardizedFromPrimitive{Primitive})(x) = x
(x::StandardizedFromPrimitive)(lattice::Lattice) = Lattice(lattice.data * x.transformation)
(x::StandardizedFromPrimitive)(coord::AbstractVector) = inv(x.transformation) * coord

const PRIM_STD_A = PrimitiveFromStandardized{ACentering}([
    1 0 0
    0 1//2 -1//2
    0 1//2 1//2
])
const PRIM_STD_C = PrimitiveFromStandardized{CCentering}([
    1//2 1//2 0
    -1//2 1//2 0
    0 0 1
])
const PRIM_STD_R = PrimitiveFromStandardized{RhombohedralCentering}(
    [
        2//3 -1//3 -1//3
        1//3 1//3 -2//3
        1//3 1//3 1//3
    ],
)
const PRIM_STD_I = PrimitiveFromStandardized{BodyCentering}(
    [
        -1//2 1//2 1//2
        1//2 -1//2 1//2
        1//2 1//2 -1//2
    ],
)
const PRIM_STD_F = PrimitiveFromStandardized{FaceCentering}([
    0 1//2 1//2
    1//2 0 1//2
    1//2 1//2 0
])
const STD_PRIM_A = StandardizedFromPrimitive{ACentering}([
    1 0 0
    0 1 1
    0 -1 1
])
const STD_PRIM_C = StandardizedFromPrimitive{CCentering}([
    1 -1 0
    1 1 0
    0 0 1
])
const STD_PRIM_R = StandardizedFromPrimitive{RhombohedralCentering}([
    1 0 1
    -1 1 1
    0 -1 1
])
const STD_PRIM_I = StandardizedFromPrimitive{BodyCentering}([
    0 1 1
    1 0 1
    1 1 0
])
const STD_PRIM_F = StandardizedFromPrimitive{FaceCentering}([
    -1 1 1
    1 -1 1
    1 1 -1
])

Base.inv(x::PrimitiveFromStandardized) = StandardizedFromPrimitive(inv(x.transformation))
Base.inv(x::StandardizedFromPrimitive) = PrimitiveFromStandardized(inv(x.transformation))
