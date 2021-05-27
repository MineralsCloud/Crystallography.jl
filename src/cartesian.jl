using CoordinateTransformations: IdentityTransformation
using LinearAlgebra: I
using StaticArrays: FieldVector

export FractionalFromCartesian,
    CartesianFromFractional, FractionalToCartesian, CartesianToFractional

abstract type Coordinates{T} <: FieldVector{3,T} end
struct Cartesian{T} <: Coordinates{T}
    x::T
    y::T
    z::T
end
struct Fractional{T} <: Coordinates{T}
    x::T
    y::T
    z::T
end

struct CartesianFromFractional
    m::SMatrix{3,3}
end
struct FractionalFromCartesian
    m::SMatrix{3,3}
end
for T in (:CartesianFromFractional, :FractionalFromCartesian)
    eval(quote
        $T(lattice::Lattice) = $T(convert(Matrix{eltype(lattice)}, lattice))
    end)
end
# This requires the a-vector is parallel to the Cartesian x-axis.
# See https://en.wikipedia.org/wiki/Fractional_coordinates
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

(x::CartesianFromFractional)(v::AbstractVector) = x.m * v
(x::FractionalFromCartesian)(v::AbstractVector) = inv(x.m) * v

Base.inv(x::FractionalFromCartesian) = CartesianFromFractional(inv(x.m))
Base.inv(x::CartesianFromFractional) = FractionalFromCartesian(inv(x.m))
Base.:∘(x::FractionalFromCartesian, y::CartesianFromFractional) =
    x.m * y.m ≈ I ? IdentityTransformation() : error("undefined!")
