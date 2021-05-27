using StaticArrays: SHermitianCompact

export MetricTensor
export directioncosine, directionangle, distance, interplanar_spacing, cellparameters

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

"""
    cellvolume(g::MetricTensor)

Calculates the cell volume from a `MetricTensor`.
"""
cellvolume(g::MetricTensor) = sqrt(det(g.data))  # `sqrt` is always positive!

Lattice(g::MetricTensor) = Lattice(cellparameters(g))

function cellparameters(g::MetricTensor)
    data = g.data
    a2, b2, c2, ab, ac, bc =
        data[1, 1], data[2, 2], data[3, 3], data[1, 2], data[1, 3], data[2, 3]
    a, b, c = map(sqrt, (a2, b2, c2))
    γ, β, α = acos(ab / (a * b)), acos(ac / (a * c)), acos(bc / (b * c))
    return a, b, c, α, β, γ
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

dot(a::AbstractVector, g::MetricTensor, b::AbstractVector) = a' * g.data * b
norm(a::AbstractVector, g::MetricTensor) = sqrt(dot(a, g, a))
