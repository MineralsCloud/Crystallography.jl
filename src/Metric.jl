"""
# module Metric



# Examples

```jldoctest
julia>
```
"""
module Metric

using LinearAlgebra: cross, det, dot

using Crystallography:
    AbstractSpace,
    RealSpace,
    ReciprocalSpace,
    CrystalCoordinates,
    BravaisLattice,
    Triclinic,
    Monoclinic,
    Orthorhombic,
    Tetragonal,
    Cubic,
    Hexagonal,
    Trigonal,
    RhombohedralCentered

import LinearAlgebra

export MetricTensor, directioncosine, directionangle, distance, interplanar_spacing

struct MetricTensor{T}
    m::T
    function MetricTensor{T}(m) where {T}
        @assert(size(m) == (3, 3), "The metric tensor must be of size 3×3!")
        return new(m)
    end
end
MetricTensor(m::T) where {T} = MetricTensor{T}(m)
function MetricTensor(v1::AbstractVector, v2::AbstractVector, v3::AbstractVector)
    vecs = (v1, v2, v3)
    return MetricTensor(map(x -> dot(x...), Iterators.product(vecs, vecs)))
end
function MetricTensor(a::Real, b::Real, c::Real, α::Real, β::Real, γ::Real)
    g12 = a * b * cos(γ)
    g13 = a * c * cos(β)
    g23 = b * c * cos(α)
    return MetricTensor([
        a^2 g12 g13
        g12 b^2 g23
        g13 g23 c^2
    ])
end
# MetricTensor(::BravaisLattice, args...) = MetricTensor(args...)  # Triclinic
# MetricTensor{RealSpace}(::BravaisLattice{C,Monoclinic}, a, b, c, γ) where {C} =
#     MetricTensor{RealSpace}(a, b, c, π / 2, π / 2, γ)
# MetricTensor{RealSpace}(::BravaisLattice{C,Orthorhombic}, a, b, c) where {C} =
#     MetricTensor{RealSpace}(a, b, c, π / 2, π / 2, π / 2)
# MetricTensor{RealSpace}(::BravaisLattice{C,Tetragonal}, a, c) where {C} =
#     MetricTensor{RealSpace}(a, a, c, π / 2, π / 2, π / 2)
# MetricTensor{RealSpace}(::BravaisLattice{C,Cubic}, a) where {C} =
#     MetricTensor{RealSpace}(a, a, a, π / 2, π / 2, π / 2)
# MetricTensor{RealSpace}(::BravaisLattice{C,Hexagonal}, a, c) where {C} =
#     MetricTensor{RealSpace}(a, a, c, π / 2, π / 2, 2π / 3)
# MetricTensor{RealSpace}(::BravaisLattice{C,Trigonal}, a, c) where {C} =
#     MetricTensor{RealSpace}(BravaisLattice(Primitive(), Hexagonal()), a, c)
# MetricTensor{RealSpace}(::BravaisLattice{RhombohedralCentered,Hexagonal}, a, α) =
#     MetricTensor{RealSpace}(a, a, a, α, α, α)

function directioncosine(
    a::CrystalCoordinates{T},
    g::MetricTensor{T},
    b::CrystalCoordinates{T},
) where {T}
    return dot(a, g, b) / (length(a, g) * length(b, g))
end

directionangle(
    a::CrystalCoordinates{T},
    g::MetricTensor{T},
    b::CrystalCoordinates{T},
) where {T} = acos(directioncosine(a, g, b))

distance(a::CrystalCoordinates{T}, g::MetricTensor{T}, b::CrystalCoordinates{T}) where {T} =
    length(b - a, g)

interplanar_spacing(
    a::CrystalCoordinates{T},
    g::MetricTensor{T},
) where {T<:ReciprocalSpace} = 1 / length(a, g)

Base.inv(::Type{MetricTensor{T}}) where {T} = MetricTensor{inv(T)}
Base.inv(g::MetricTensor) = inv(typeof(g))(inv(g.m))

LinearAlgebra.dot(
    a::CrystalCoordinates{T},
    g::MetricTensor{T},
    b::CrystalCoordinates{T},
) where {T} = a' * g.m * b

Base.length(a::CrystalCoordinates{T}, g::MetricTensor{T}) where {T} = sqrt(dot(a, g, a))

function reciprocalof(mat::AbstractMatrix)
    @assert size(mat) == (3, 3)
    volume = abs(det(mat))
    a1, a2, a3 = mat[1, :], mat[2, :], mat[3, :]
    return 2π / volume * [cross(a2, a3) cross(a3, a1) cross(a1, a2)]
end # function reciprocalof

end
