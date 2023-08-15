module Symmetry

using StaticArrays: MMatrix, Size

import StaticArrays: similar_type

export SeitzOperator

struct SeitzOperator{T} <: AbstractMatrix{T}
    data::MMatrix{4,4,T,16}
end
SeitzOperator(matrix::AbstractMatrix) = SeitzOperator(MMatrix{4,4}(matrix))
SeitzOperator{T}(::UndefInitializer, dims) where {T} =
    SeitzOperator(MMatrix{4,4,T,16}(undef, dims))

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
