using LinearAlgebra: Diagonal

import StaticArrays: similar_type

# See https://juliaarrays.github.io/StaticArrays.jl/dev/pages/api/#StaticArraysCore.Size
# and http://docs.julialang.org/en/v1/base/base/#Base.Val
struct Size{x} end
Size(x) = Size{x}()

apply(::Size, ::SeitzOperator, 𝐫::AbstractVector) = throw(
    DimensionMismatch(
        "`SeitzOperator` can be only applied onto vectors of lengths 3 or 4!"
    ),
)
function apply(::Size{(3,)}, op::SeitzOperator, 𝐫::AbstractVector)
    𝐑, 𝐭 = getpointsymmetry(op), gettranslation(op)
    return 𝐑 * 𝐫 + 𝐭
end
function apply(::Size{(4,)}, op::SeitzOperator, 𝐫::AbstractVector)
    @assert 𝐫[end] == 1
    return op * 𝐫
end

# Much faster than writing `SeitzOperator(𝐑 * 𝐒, 𝐑 * 𝐮 + 𝐭)`
Base.:*(op₁::SeitzOperator, op₂::SeitzOperator) = SeitzOperator(parent(op₁) * parent(op₂))

similar_type(::SeitzOperator, ::Type{T}) where {T} = similar_type(SeitzOperator, T)
similar_type(::Type{<:SeitzOperator}, ::Type{T}) where {T} = SeitzOperator{T}

# See https://github.com/JuliaLang/julia/blob/v1.10.0-beta1/stdlib/LinearAlgebra/src/uniformscaling.jl#L130-L131
Base.one(::Type{SeitzOperator{T}}) where {T} =
    SeitzOperator(MMatrix{4,4}(Diagonal(fill(one(T), 4))))
Base.one(op::SeitzOperator) = one(typeof(op))

# See https://github.com/JuliaLang/julia/blob/v1.10.0-beta1/stdlib/LinearAlgebra/src/uniformscaling.jl#L132-L133
Base.oneunit(::Type{SeitzOperator{T}}) where {T} =
    SeitzOperator(MMatrix{4,4}(Diagonal(fill(oneunit(T), 4))))
Base.oneunit(op::SeitzOperator) = oneunit(typeof(op))

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

# See https://github.com/JuliaArrays/StaticArraysCore.jl/blob/v1.4.2/src/StaticArraysCore.jl#L397-L398
struct SeitzOperatorStyle <: Broadcast.AbstractArrayStyle{2} end
SeitzOperatorStyle(::Val{N}) where {N} = SeitzOperatorStyle()

Base.BroadcastStyle(::Type{<:SeitzOperator}) = SeitzOperatorStyle()

Base.similar(::Broadcast.Broadcasted{SeitzOperatorStyle}, ::Type{T}) where {T} =
    similar(SeitzOperator{T})
# Override https://github.com/JuliaLang/julia/blob/v1.10.0-beta1/base/abstractarray.jl#L874
Base.similar(::Type{SeitzOperator{T}}, dims::Dims) where {T} = SeitzOperator{T}(undef)
# Override https://github.com/JuliaLang/julia/blob/v1.10.0-beta1/base/abstractarray.jl#L827
Base.similar(::SeitzOperator, ::Type{T}) where {T} = SeitzOperator{T}(undef)
