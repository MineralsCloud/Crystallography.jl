using StaticArrays: SHermitianCompact

export MetricTensor
export directioncosine, directionangle, distance, interplanar_spacing, cellparameters

struct MetricTensor{T} <: AbstractMatrix{T}
    data::SHermitianCompact{3,T,6}
end
MetricTensor(m::AbstractMatrix) = MetricTensor(SHermitianCompact{3}(m))
function MetricTensor(𝐚::AbstractVector, 𝐛::AbstractVector, 𝐜::AbstractVector)
    vecs = (𝐚, 𝐛, 𝐜)
    return MetricTensor([dot(vecs[i], vecs[j]) for i in 1:3, j in 1:3])
end
function MetricTensor(a, b, c, α, β, γ)
    g12 = a * b * cos(γ)
    g13 = a * c * cos(β)
    g23 = b * c * cos(α)
    return MetricTensor(SHermitianCompact(SVector(a^2, g12, g13, b^2, g23, c^2)))
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

directioncosine(𝐚::AbstractVector, g::MetricTensor, 𝐛::AbstractVector) =
    dot(𝐚, g, 𝐛) / (norm(𝐚, g) * norm(𝐛, g))

directionangle(𝐚::AbstractVector, g::MetricTensor, 𝐛::AbstractVector) =
    acos(directioncosine(𝐚, g, 𝐛))

distance(𝐚::AbstractVector, g::MetricTensor, 𝐛::AbstractVector) = norm(𝐛 - 𝐚, g)

interplanar_spacing(𝐚::AbstractVector, g::MetricTensor) = 1 / norm(𝐚, g)

Base.size(::MetricTensor) = (3, 3)

Base.getindex(A::MetricTensor, I::Vararg{Int}) = getindex(A.data, I...)

Base.inv(g::MetricTensor) = MetricTensor(inv(g.data))

dot(𝐚::AbstractVector, g::MetricTensor, 𝐛::AbstractVector) = 𝐚' * g.data * 𝐛
norm(𝐚::AbstractVector, g::MetricTensor) = sqrt(dot(𝐚, g, 𝐚))
