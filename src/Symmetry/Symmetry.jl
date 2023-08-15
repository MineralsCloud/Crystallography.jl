# module Symmetry

using LinearAlgebra: I, Diagonal, diagm, det
using StaticArrays: MMatrix

import StaticArrays: similar_type

export SeitzOperator, istranslation, ispointsymmetry, gettranslation, getpointsymmetry

"""
    SeitzOperator(𝐑::AbstractMatrix, 𝐭::AbstractVector)
    SeitzOperator(𝐑::AbstractMatrix)
    SeitzOperator(𝐭::AbstractVector)

Construct a Seitz operator from matrices, vectors, or both.

The operator is defined by the following equation:

```math
\\left[
\\begin{array}{ccc|c}
& & & \\\\
& \\mathbf{R} & & \\mathbf{t} \\\\
& & & \\\\
\\hline 0 & 0 & 0 & 1
\\end{array}
\\right],
```

where ``\\mathbf{R}`` is a point group operation and ``\\mathbf{t}`` is a translation.
"""
struct SeitzOperator{T} <: AbstractMatrix{T}
    data::MMatrix{4,4,T,16}
end
# See https://github.com/JuliaArrays/StaticArraysCore.jl/blob/v1.4.2/src/StaticArraysCore.jl#L195-L198
SeitzOperator{T}(::UndefInitializer) where {T} = SeitzOperator(MMatrix{4,4,T,16}(undef))
function SeitzOperator(𝐑::AbstractMatrix)
    @assert size(𝐑) == (3, 3)
    data = diagm(ones(eltype(𝐑), 4))
    data[1:3, 1:3] = 𝐑
    return SeitzOperator(MMatrix{4,4}(data))
end
function SeitzOperator(𝐭::AbstractVector)
    @assert length(𝐭) == 3
    data = diagm(ones(eltype(𝐭), 4))
    data[1:3, 4] = 𝐭
    return SeitzOperator(MMatrix{4,4}(data))
end
SeitzOperator(𝐑::AbstractMatrix, 𝐭::AbstractVector) =
    SeitzOperator(MMatrix{4,4}(vcat(hcat(𝐑, 𝐭), [zeros(eltype(𝐑), 3)... one(eltype(𝐭))])))

function (op::SeitzOperator)(𝐫::AbstractVector)
    𝐑, 𝐭 = getpointsymmetry(op), gettranslation(op)
    return 𝐑 * 𝐫 + 𝐭
end

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

gettranslation(op::SeitzOperator) = op[1:3, 4]

getpointsymmetry(op::SeitzOperator) = op[1:3, 1:3]

similar_type(::SeitzOperator, ::Type{T}) where {T} = similar_type(SeitzOperator, T)
similar_type(::Type{<:SeitzOperator}, ::Type{T}) where {T} = SeitzOperator{T}

Base.size(::SeitzOperator) = (4, 4)

Base.parent(op::SeitzOperator) = op.data

Base.getindex(op::SeitzOperator, i) = getindex(parent(op), i)

Base.setindex!(op::SeitzOperator, v, i) = setindex!(parent(op), v, i)

Base.IndexStyle(::Type{SeitzOperator{T}}) where {T} = IndexLinear()

function Base.inv(op::SeitzOperator)
    𝐑, 𝐭 = getpointsymmetry(op), gettranslation(op)
    𝐑⁻¹ = inv(𝐑)
    return SeitzOperator(𝐑⁻¹, -𝐑⁻¹ * 𝐭)
end

# See https://github.com/JuliaLang/julia/blob/v1.10.0-beta1/stdlib/LinearAlgebra/src/uniformscaling.jl#L130-L131
Base.one(::Type{SeitzOperator{T}}) where {T} =
    SeitzOperator(MMatrix{4,4}(Diagonal(fill(one(T), 4))))
Base.one(op::SeitzOperator) = one(typeof(op))

# See https://github.com/JuliaLang/julia/blob/v1.10.0-beta1/stdlib/LinearAlgebra/src/uniformscaling.jl#L132-L133
Base.oneunit(::Type{SeitzOperator{T}}) where {T} =
    SeitzOperator(MMatrix{4,4}(Diagonal(fill(oneunit(T), 4))))
Base.oneunit(op::SeitzOperator) = oneunit(typeof(op))

# See https://github.com/JuliaArrays/StaticArraysCore.jl/blob/v1.4.2/src/StaticArraysCore.jl#L397-L398
struct SeitzOperatorStyle <: Broadcast.AbstractArrayStyle{2} end
SeitzOperatorStyle(::Val{N}) where {N} = SeitzOperatorStyle()

Base.BroadcastStyle(::Type{<:SeitzOperator}) = SeitzOperatorStyle()

Base.similar(::Broadcast.Broadcasted{SeitzOperatorStyle}, ::Type{T}) where {T} =
    similar(SeitzOperator{T})
# Override https://github.com/JuliaLang/julia/blob/v1.10.0-beta1/base/abstractarray.jl#L874
Base.similar(::Type{SeitzOperator{T}}, dims::Dims) where {T} = SeitzOperator{T}(undef)
# Override https://github.com/JuliaLang/julia/blob/v1.10.0-beta1/base/abstractarray.jl#L839C1-L839C93
Base.similar(::SeitzOperator, ::Type{T}, ::Dims{N}) where {T,N} = SeitzOperator{T}(undef)

# end
