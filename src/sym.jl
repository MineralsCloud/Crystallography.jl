module Symmetry

using LinearAlgebra: I, diagm, det
using StaticArrays: MMatrix

import StaticArrays: similar_type

export SeitzOperator, istranslation, ispointsymmetry

struct SeitzOperator{T} <: AbstractMatrix{T}
    data::MMatrix{4,4,T,16}
end
SeitzOperator{T}(::UndefInitializer, dims) where {T} =
    SeitzOperator(MMatrix{4,4,T,16}(undef, dims))
function SeitzOperator(m::AbstractMatrix)
    @assert size(m) == (3, 3)
    x = diagm(ones(eltype(m), 4))
    x[1:3, 1:3] = m
    return SeitzOperator(x)
end
function SeitzOperator(t::AbstractVector)
    @assert length(t) == 3
    data = diagm(ones(eltype(t), 4))
    data[1:3, 4] = t
    return SeitzOperator(data)
end
SeitzOperator(m::AbstractMatrix, t::AbstractVector) =
    SeitzOperator(MMatrix{4,4}(vcat(hcat(m, t), [zeros(eltype(m), 3)... one(eltype(t))])))

function istranslation(op::SeitzOperator)
    m = op.data
    if m[1:3, 1:3] != I || !(iszero(m[4, 1:3]) && isone(m[4, 4]))
        return false
    end
    return true
end

function ispointsymmetry(op::SeitzOperator)
    m = op.data
    if iszero(m[4, 1:3]) &&
        iszero(m[1:3, 4]) &&
        isone(m[4, 4]) &&
        abs(det(m[1:3, 1:3])) == 1
        return true
    else
        false
    end
end

similar_type(::SeitzOperator, ::Type{T}) where {T} = similar_type(SeitzOperator, T)
similar_type(::Type{<:SeitzOperator}, ::Type{T}) where {T} = SeitzOperator{T}

Base.size(::SeitzOperator) = (4, 4)

Base.parent(op::SeitzOperator) = op.data

Base.getindex(op::SeitzOperator, i) = getindex(parent(op), i)

Base.setindex!(op::SeitzOperator, v, i) = setindex!(parent(op), v, i)

Base.IndexStyle(::Type{SeitzOperator{T}}) where {T} = IndexLinear()

Base.BroadcastStyle(::Type{<:SeitzOperator}) = Broadcast.ArrayStyle{SeitzOperator}()

Base.similar(
    bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{SeitzOperator}}, ::Type{T}
) where {T} = similar(SeitzOperator{T}, axes(bc))

end
