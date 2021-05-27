export RealSpace, ReciprocalSpace, MetricTensor, CrystalFromCartesian, CartesianFromCrystal
export directioncosine,
    directionangle, distance, interplanar_spacing, reciprocal, cellparameters

abstract type AbstractSpace end
struct RealSpace <: AbstractSpace end
struct ReciprocalSpace <: AbstractSpace end

struct MetricTensor{T} <: AbstractMatrix{T}
    data::SHermitianCompact{3,T}
end
MetricTensor(m::AbstractMatrix) = MetricTensor(SHermitianCompact{3}(m))
function MetricTensor(v1::AbstractVector, v2::AbstractVector, v3::AbstractVector)
    vecs = (v1, v2, v3)
    return MetricTensor([dot(vecs[i], vecs[j]) for i in 1:3, j in 1:3])
end
function MetricTensor(a, b, c, α, β, γ)
    g12 = a * b * cos(γ)
    g13 = a * c * cos(β)
    g23 = b * c * cos(α)
    return MetricTensor(SHermitianCompact(SVector(a^2, g12, g13, b^2, g23, c^2)))
end

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

"""
    cellvolume(g::MetricTensor)

Calculates the cell volume from a `MetricTensor`.
"""
cellvolume(g::MetricTensor) = sqrt(det(g.data))  # `sqrt` is always positive!

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
Lattice(g::MetricTensor) = Lattice(cellparameters(g))

function cellparameters(g::MetricTensor)
    data = g.data
    a2, b2, c2, ab, ac, bc =
        data[1, 1], data[2, 2], data[3, 3], data[1, 2], data[1, 3], data[2, 3]
    a, b, c = map(sqrt, (a2, b2, c2))
    γ, β, α = acos(ab / (a * b)), acos(ac / (a * c)), acos(bc / (b * c))
    return a, b, c, α, β, γ
end

function reciprocal(lattice::Lattice)
    volume = cellvolume(lattice)
    a1, a2, a3 = basis_vectors(lattice)
    return 1 / volume * [cross(a2, a3) cross(a3, a1) cross(a1, a2)]
end

directioncosine(a::AbstractVector, g::MetricTensor, b::AbstractVector) =
    dot(a, g, b) / (norm(a, g) * norm(b, g))

directionangle(a::AbstractVector, g::MetricTensor, b::AbstractVector) =
    acos(directioncosine(a, g, b))

distance(a::AbstractVector, g::MetricTensor, b::AbstractVector) = norm(b - a, g)

interplanar_spacing(a::AbstractVector, g::MetricTensor) = 1 / norm(a, g)

Base.size(::MetricTensor) = (3, 3)

Base.getindex(A::MetricTensor, I::Vararg{Int}) = getindex(A.data, I...)

Base.inv(g::MetricTensor) = MetricTensor(inv(g.data))
Base.inv(x::CrystalFromCartesian) = CartesianFromCrystal(inv(x.m))
Base.inv(x::CartesianFromCrystal) = CrystalFromCartesian(inv(x.m))
Base.:∘(x::CrystalFromCartesian, y::CartesianFromCrystal) =
    x.m * y.m ≈ I ? IdentityTransformation() : error("undefined!")

dot(a::AbstractVector, g::MetricTensor, b::AbstractVector) = a' * g.data * b
norm(a::AbstractVector, g::MetricTensor) = sqrt(dot(a, g, a))
