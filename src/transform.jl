using CoordinateTransformations: IdentityTransformation
using LinearAlgebra: I

export FractionalFromCartesian,
    CartesianFromFractional, FractionalToCartesian, CartesianToFractional

struct CartesianFromFractional
    tf::SMatrix{3,3}
end
struct FractionalFromCartesian
    tf::SMatrix{3,3}
end
# This requires the a-vector is parallel to the Cartesian x-axis.
# See https://en.wikipedia.org/wiki/Fractional_coordinates
CartesianFromFractional(lattice::Lattice) = CartesianFromFractional(lattice.data)
function CartesianFromFractional(a, b, c, α, β, γ)
    Ω = cellvolume(a, b, c, α, β, γ)
    b_sinγ, b_cosγ = b .* sincos(γ)
    return CartesianFromFractional(
        [
            a b_cosγ c*cos(β)
            0 b_sinγ c*_auxiliary(α, β, γ)
            0 0 Ω/(a*b_sinγ)
        ],
    )
end
FractionalFromCartesian(lattice::Lattice) = FractionalFromCartesian(inv(lattice.data))
function FractionalFromCartesian(a, b, c, α, β, γ)
    Ω = cellvolume(a, b, c, α, β, γ)
    b_sinγ = b * sin(γ)
    return FractionalFromCartesian(
        [
            1/a -cot(γ)/a -b*c*_auxiliary(β, α, γ)/Ω
            0 1/b_sinγ -a*c*_auxiliary(α, β, γ)/Ω
            0 0 a*b_sinγ/Ω
        ],
    )
end
const FractionalToCartesian = CartesianFromFractional
const CartesianToFractional = FractionalFromCartesian

# This is a helper function and should not be exported!
_auxiliary(α, β, γ) = (cos(α) - cos(β) * cos(γ)) / sin(γ)

(x::Union{CartesianFromFractional,FractionalFromCartesian})(v) = x.tf * collect(v)

Base.inv(x::FractionalFromCartesian) = CartesianFromFractional(inv(x.tf))
Base.inv(x::CartesianFromFractional) = FractionalFromCartesian(inv(x.tf))
Base.:∘(x::FractionalFromCartesian, y::CartesianFromFractional) =
    x.tf * y.tf ≈ I ? IdentityTransformation() : error("undefined!")

# Idea from https://spglib.github.io/spglib/definition.html#transformation-to-the-primitive-cell
struct StandardizedFromPrimitive{T<:Centering}
    transformation::SMatrix{3,3}
end
struct PrimitiveFromStandardized{T<:Centering}
    transformation::SMatrix{3,3}
end

(::StandardizedFromPrimitive{Primitive})(x) = x
(x::StandardizedFromPrimitive)(lattice::Lattice) = Lattice(lattice.data * x.transformation)
(x::StandardizedFromPrimitive)(coord::AbstractVector) = inv(x.transformation) * coord

(::PrimitiveFromStandardized{Primitive})(x) = x
(x::PrimitiveFromStandardized)(lattice::Lattice) = Lattice(lattice.data * x.transformation)
(x::PrimitiveFromStandardized)(coord::AbstractVector) = inv(x.transformation) * coord

const PRIM_STD_A = StandardizedFromPrimitive{ACentering}([
    1 0 0
    0 1//2 -1//2
    0 1//2 1//2
])
const PRIM_STD_C = StandardizedFromPrimitive{CCentering}([
    1//2 1//2 0
    -1//2 1//2 0
    0 0 1
])
const PRIM_STD_R = StandardizedFromPrimitive{RhombohedralCentering}(
    [
        2//3 -1//3 -1//3
        1//3 1//3 -2//3
        1//3 1//3 1//3
    ],
)
const PRIM_STD_I = StandardizedFromPrimitive{BodyCentering}(
    [
        -1//2 1//2 1//2
        1//2 -1//2 1//2
        1//2 1//2 -1//2
    ],
)
const PRIM_STD_F = StandardizedFromPrimitive{FaceCentering}([
    0 1//2 1//2
    1//2 0 1//2
    1//2 1//2 0
])
const STD_PRIM_A = PrimitiveFromStandardized{ACentering}([
    1 0 0
    0 1 1
    0 -1 1
])
const STD_PRIM_C = PrimitiveFromStandardized{CCentering}([
    1 -1 0
    1 1 0
    0 0 1
])
const STD_PRIM_R = PrimitiveFromStandardized{RhombohedralCentering}([
    1 0 1
    -1 1 1
    0 -1 1
])
const STD_PRIM_I = PrimitiveFromStandardized{BodyCentering}([
    0 1 1
    1 0 1
    1 1 0
])
const STD_PRIM_F = PrimitiveFromStandardized{FaceCentering}([
    -1 1 1
    1 -1 1
    1 1 -1
])

Base.inv(x::StandardizedFromPrimitive) = PrimitiveFromStandardized(inv(x.transformation))
Base.inv(x::PrimitiveFromStandardized) = StandardizedFromPrimitive(inv(x.transformation))
