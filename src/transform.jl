using CoordinateTransformations: IdentityTransformation, Transformation
using LinearAlgebra: I

export FractionalFromCartesian,
    CartesianFromFractional,
    FractionalToCartesian,
    CartesianToFractional,
    StandardizedFromPrimitive,
    PrimitiveFromStandardized,
    PrimitiveToStandardized,
    StandardizedToPrimitive

struct CartesianFromFractional{T} <: Transformation
    tf::SMatrix{3,3,T,9}
end
struct FractionalFromCartesian{T} <: Transformation
    tf::SMatrix{3,3,T,9}
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
Base.:∘(x::CartesianFromFractional, y::FractionalFromCartesian) = ∘(y, x)
Base.:∘(x::FractionalFromCartesian, y::CartesianFromFractional) =
    x.tf * y.tf ≈ I ? IdentityTransformation() : error("undefined!")

# Idea from https://spglib.github.io/spglib/definition.html#transformation-to-the-primitive-cell
struct StandardizedFromPrimitive{T} <: Transformation
    tf::SMatrix{3,3,T,9}
end
struct PrimitiveFromStandardized{T} <: Transformation
    tf::SMatrix{3,3,T,9}
end
PrimitiveFromStandardized(tf::AbstractMatrix) = PrimitiveFromStandardized{eltype(tf)}(tf)
StandardizedFromPrimitive(tf::AbstractMatrix) = StandardizedFromPrimitive{eltype(tf)}(tf)
const PrimitiveToStandardized = StandardizedFromPrimitive
const StandardizedToPrimitive = PrimitiveFromStandardized

(::StandardizedFromPrimitive{Primitive})(v) = v
(::PrimitiveFromStandardized{Primitive})(v) = v
(x::Union{StandardizedFromPrimitive,PrimitiveFromStandardized})(v) = x.tf * collect(v)

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

Base.inv(x::StandardizedFromPrimitive) = PrimitiveFromStandardized(inv(x.tf))
Base.inv(x::PrimitiveFromStandardized) = StandardizedFromPrimitive(inv(x.tf))
Base.:∘(x::PrimitiveFromStandardized, y::StandardizedFromPrimitive) = ∘(y, x)
Base.:∘(x::StandardizedFromPrimitive, y::PrimitiveFromStandardized) =
    x.tf * y.tf ≈ I ? IdentityTransformation() : error("undefined!")

function Base.show(
    io::IO,
    x::Union{
        CartesianFromFractional,
        FractionalFromCartesian,
        StandardizedFromPrimitive,
        PrimitiveFromStandardized,
    },
)
    if get(io, :compact, false)
        print(io, string(typeof(x)), '(')
        print(io, x.tf, ')')
    else
        println(io, string(nameof(typeof(x))))
        for row in eachrow(x.tf)
            print(io, " ")
            println(io, row)
        end
    end
end
