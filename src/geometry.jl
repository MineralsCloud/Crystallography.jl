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
MetricTensor(p::CellParameters) = MetricTensor(p...)

struct Miller{S<:AbstractSpace} <: AbstractVector{Int}
    data::SVector{3,Int}
    Miller{S}(v) where {S} = new(iszero(v) ? v : v .÷ gcd(v))
end
Miller{S}(i, j, k) where {S} = Miller{S}([i, j, k])

struct MillerBravais{S<:AbstractSpace} <: AbstractVector{Int}
    data::SVector{4,Int}
    function MillerBravais{S}(v) where {S}
        @assert(
            v[3] == -v[1] - v[2],
            "the 3rd index of `MillerBravais` should equal to the negation of the first two!"
        )
        return new(iszero(v) ? v : v .÷ gcd(v))
    end
end
MillerBravais{S}(i, j, k, l) where {S} = MillerBravais{S}([i, j, k, l])

# This is a helper type and should not be exported!
const INDICES = Union{Miller,MillerBravais}

# This is a helper function and should not be exported!
function _indices_str(r::Regex, s::AbstractString, ::Type{T}) where {T<:INDICES}
    m = match(r, strip(s))
    isnothing(m) && error("not a valid expression!")
    brackets = first(m.captures) * last(m.captures)
    x = (parse(Int, x) for x in m.captures[2:(end-1)])
    if brackets ∈ ("()", "{}")
        return T{ReciprocalSpace}(x...)
    elseif brackets ∈ ("[]", "<>")
        return T{RealSpace}(x...)
    else
        error("not a valid expression!")
    end
end # function _indices_str

macro m_str(s)
    r = r"([({[<])\s*([-+]?[0-9]+)[\s,]+([-+]?[0-9]+)[\s,]+([-+]?[0-9]+)[\s,]*([>\]})])"
    _indices_str(r, s, Miller)
end

macro mb_str(s)
    r = r"([({[<])\s*([-+]?[0-9]+)[\s,]+([-+]?[0-9]+)[\s,]+([-+]?[0-9]+)[\s,]+([-+]?[0-9]+)[\s,]*([>\]})])"
    _indices_str(r, s, MillerBravais)
end

"""
    cellvolume(g::MetricTensor)

Calculates the cell volume from a `MetricTensor`.
"""
cellvolume(g::MetricTensor) = sqrt(det(g.data))  # `sqrt` is always positive!

directioncosine(a::CrystalCoord, g::MetricTensor, b::CrystalCoord) =
    dot(a, g, b) / (norm(a, g) * norm(b, g))

directionangle(a::CrystalCoord, g::MetricTensor, b::CrystalCoord) = acos(directioncosine(a, g, b))

distance(a::CrystalCoord, g::MetricTensor, b::CrystalCoord) = norm(b - a, g)

interplanar_spacing(a::CrystalCoord, g::MetricTensor) = 1 / norm(a, g)

Base.size(::Union{MetricTensor,Lattice}) = (3, 3)
Base.size(::Miller) = (3,)
Base.size(::MillerBravais) = (4,)

Base.getindex(A::MetricTensor, I::Vararg{Int}) = getindex(A.data, I...)
Base.getindex(
    A::Union{Miller,MillerBravais,CellParameters,Lattice,CellParameters},
    i::Int,
) = getindex(A.data, i)

Base.inv(g::MetricTensor) = MetricTensor(inv(SymPy.N(g.data)))

Base.convert(::Type{T}, x::T) where {T<:INDICES} = x
Base.convert(::Type{Miller{T}}, mb::MillerBravais{T}) where {T<:RealSpace} =
    Miller{T}(2 * mb[1] + mb[2], 2 * mb[2] + mb[1], mb[4])
Base.convert(::Type{Miller{T}}, mb::MillerBravais{T}) where {T<:ReciprocalSpace} =
    Miller{T}(mb[1], mb[2], mb[4])
Base.convert(::Type{MillerBravais{T}}, m::Miller{T}) where {T<:RealSpace} =
    MillerBravais{T}(2 * m[1] - m[2], 2 * m[2] - m[1], -(m[1] + m[2]), 3 * m[3])
Base.convert(::Type{MillerBravais{T}}, m::Miller{T}) where {T<:ReciprocalSpace} =
    MillerBravais{T}(m[1], m[2], -(m[1] + m[2]), m[3])
function Base.convert(::Type{CellParameters}, g::MetricTensor)
    data = g.data
    a2, b2, c2, ab, ac, bc =
        data[1, 1], data[2, 2], data[3, 3], data[1, 2], data[1, 3], data[2, 3]
    a, b, c = map(sqrt, (a2, b2, c2))
    γ, β, α = acos(ab / (a * b)), acos(ac / (a * c)), acos(bc / (b * c))
    return CellParameters(a, b, c, α, β, γ)
end # function Base.convert
Base.convert(::Type{Lattice}, g::MetricTensor) = Lattice(convert(CellParameters, g))

LinearAlgebra.dot(a::CrystalCoord, g::MetricTensor, b::CrystalCoord) = a' * g.data * b
LinearAlgebra.norm(a::CrystalCoord, g::MetricTensor) = sqrt(dot(a, g, a))
