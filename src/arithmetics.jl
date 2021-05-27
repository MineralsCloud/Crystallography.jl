export RealSpace, ReciprocalSpace, CrystalFromCartesian, CartesianFromCrystal

abstract type AbstractSpace end
struct RealSpace <: AbstractSpace end
struct ReciprocalSpace <: AbstractSpace end

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

function Lattice(a, b, c, α, β, γ)
    # From https://github.com/LaurentRDC/crystals/blob/dbb3a92/crystals/lattice.py#L321-L354
    v = cellvolume(1, 1, 1, α, β, γ)
    # reciprocal lattice
    a_recip = sin(α) / (a * v)
    csg = (cos(α) * cos(β) - cos(γ)) / (sin(α) * sin(β))
    sg = sqrt(1 - csg^2)
    a1 = [1 / a_recip, -csg / sg / a_recip, cos(β) * a]
    a2 = [0, b * sin(α), b * cos(α)]
    a3 = [0, 0, c]
    return Lattice(a1, a2, a3)
end

function reciprocal(lattice::Lattice)
    volume = cellvolume(lattice)
    a1, a2, a3 = basis_vectors(lattice)
    return 1 / volume * [cross(a2, a3) cross(a3, a1) cross(a1, a2)]
end

include("metric.jl")

Base.inv(x::CrystalFromCartesian) = CartesianFromCrystal(inv(x.m))
Base.inv(x::CartesianFromCrystal) = CrystalFromCartesian(inv(x.m))
Base.:∘(x::CrystalFromCartesian, y::CartesianFromCrystal) =
    x.m * y.m ≈ I ? IdentityTransformation() : error("undefined!")
