module Symmetry

using LinearAlgebra: I, diagm, det
using StaticArrays: MMatrix

import StaticArrays: similar_type

export SeitzOperator, istranslation, ispointsymmetry, isidentity

struct SeitzOperator{T} <: AbstractMatrix{T}
    data::MMatrix{4,4,T,16}
end
SeitzOperator{T}(::UndefInitializer, dims) where {T} =
    SeitzOperator(MMatrix{4,4,T,16}(undef, dims))
function SeitzOperator(matrix::AbstractMatrix)
    @assert size(matrix) == (3, 3)
    data = diagm(ones(eltype(matrix), 4))
    data[1:3, 1:3] = matrix
    return SeitzOperator(data)
end
function SeitzOperator(vector::AbstractVector)
    @assert length(vector) == 3
    data = diagm(ones(eltype(vector), 4))
    data[1:3, 4] = vector
    return SeitzOperator(data)
end
SeitzOperator(𝐑::AbstractMatrix, 𝐭::AbstractVector) =
    SeitzOperator(MMatrix{4,4}(vcat(hcat(𝐑, 𝐭), [zeros(eltype(𝐑), 3)... one(eltype(𝐭))])))

function istranslation(op::SeitzOperator)
    if op[1:3, 1:3] != I || !(iszero(op[4, 1:3]) && isone(op[4, 4]))
        return false
    end
    return true
end

function ispointsymmetry(op::SeitzOperator)
    if iszero(op[4, 1:3]) &&
        iszero(op[1:3, 4]) &&
        isone(op[4, 4]) &&
        abs(det(op[1:3, 1:3])) == 1
        return true
    else
        false
    end
end

isidentity(op::SeitzOperator) = op == I

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

function Base.inv(op::SeitzOperator)
    𝐑, 𝐭 = getpointsymmetry(op), gettranslation(op)
    𝐑⁻¹ = inv(𝐑)
    return SeitzOperator(𝐑⁻¹, -𝐑⁻¹ * 𝐭)
end

end
