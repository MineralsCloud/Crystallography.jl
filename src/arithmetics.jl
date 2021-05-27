export RealSpace, ReciprocalSpace, CrystalFromCartesian, CartesianFromCrystal

abstract type AbstractSpace end
struct RealSpace <: AbstractSpace end
struct ReciprocalSpace <: AbstractSpace end

"""
    cellvolume(a, b, c, α, β, γ)

Calculates the cell volume from 6 cell parameters.
"""
cellvolume(a, b, c, α, β, γ) =
    a * b * c * sqrt(sin(α)^2 - cos(β)^2 - cos(γ)^2 + 2 * cos(α) * cos(β) * cos(γ))
"""
    cellvolume(l::Lattice)
    cellvolume(c::Cell)

Calculates the cell volume from a `Lattice` or a `Cell`.
"""
cellvolume(lattice::Lattice) = abs(det(lattice.data))
cellvolume(cell::Cell) = cellvolume(cell.lattice)

include("miller.jl")

struct CartesianFromCrystal
    m::SMatrix{3,3}
end
struct CrystalFromCartesian
    m::SMatrix{3,3}
end
for T in (:CartesianFromCrystal, :CrystalFromCartesian)
    eval(quote
        $T(lattice::Lattice) = $T(convert(Matrix{eltype(lattice)}, lattice))
    end)
end
function CartesianFromCrystal(a, b, c, α, β, γ)
    v = cellvolume(a, b, c, α, β, γ)
    x, y = sin(γ), cos(γ)
    return CartesianFromCrystal(
        [
            a b*y -c/x^2*(_F(γ, α, β)+_F(β, α, γ)*y)
            0 b*x -c*_F(β, α, γ)/x
            0 0 v/(a*b*x)
        ],
    )
end
function CrystalFromCartesian(a, b, c, α, β, γ)  # This is wrong
    v = cellvolume(a, b, c, α, β, γ)
    x = sin(γ)
    return CrystalFromCartesian(
        [
            1/a -1/(a*tan(γ)) b*c*_F(γ, α, β)/(v*x)
            0 1/(b*x) a*c*_F(β, γ, α)/(v*x)
            0 0 a*b*x/v
        ],
    )
end

# This is a helper function and should not be exported!
_F(α, β, γ) = cos(α) * cos(β) - cos(γ)

(x::CartesianFromCrystal)(v::AbstractVector) = x.m * v
(x::CrystalFromCartesian)(v::AbstractVector) = inv(x.m) * v

include("metric.jl")
include("reciprocal.jl")

Base.inv(x::CrystalFromCartesian) = CartesianFromCrystal(inv(x.m))
Base.inv(x::CartesianFromCrystal) = CrystalFromCartesian(inv(x.m))
Base.:∘(x::CrystalFromCartesian, y::CartesianFromCrystal) =
    x.m * y.m ≈ I ? IdentityTransformation() : error("undefined!")
