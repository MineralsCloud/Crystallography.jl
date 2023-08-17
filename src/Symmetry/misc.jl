using LinearAlgebra: Diagonal

import Base: +, -, *, /

# See https://juliaarrays.github.io/StaticArrays.jl/dev/pages/api/#StaticArraysCore.Size
# and http://docs.julialang.org/en/v1/base/base/#Base.Val
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

# Array interface
Base.parent(op::SeitzOperator) = op.data

Base.getindex(op::SeitzOperator, i::Int) = getindex(parent(op), i)

Base.setindex!(op::SeitzOperator, v, i::Int) = setindex!(parent(op), v, i)

# Customizing broadcasting
# See https://github.com/JuliaArrays/StaticArraysCore.jl/blob/v1.4.2/src/StaticArraysCore.jl#L397-L398
# and https://github.com/JuliaLang/julia/blob/v1.10.0-beta1/stdlib/LinearAlgebra/src/structuredbroadcast.jl#L7-L14
struct SeitzOperatorStyle <: Broadcast.AbstractArrayStyle{2} end
SeitzOperatorStyle(::Val{2}) = SeitzOperatorStyle()
SeitzOperatorStyle(::Val{N}) where {N} = Broadcast.DefaultArrayStyle{N}()

Base.BroadcastStyle(::Type{<:SeitzOperator}) = SeitzOperatorStyle()

Base.similar(::Broadcast.Broadcasted{SeitzOperatorStyle}, ::Type{T}) where {T} =
    similar(SeitzOperator{T})
# Override https://github.com/JuliaArrays/StaticArrays.jl/blob/v1.6.2/src/abstractarray.jl#L129
function Base.similar(op::SeitzOperator, ::Type{T}, _size::Size) where {T}
    if _size == size(op)
        one(SeitzOperator{T})
    else
        return similar(Array(op), T, _size)
    end
end
# Override https://github.com/JuliaLang/julia/blob/v1.10.0-beta2/base/abstractarray.jl#L839
function Base.similar(op::SeitzOperator, ::Type{T}, dims::Dims) where {T}
    if dims == size(op)
        one(SeitzOperator{T})
    else
        return similar(Array(op), T, dims)
    end
end
function Base.similar(::Type{<:SeitzOperator}, ::Type{T}, s::Size) where {T}
    if s == (4, 4)
        one(SeitzOperator{T})
    else
        return Array{T}(undef, Tuple(s))
    end
end
function Base.similar(::Type{<:SeitzOperator}, ::Type{T}, dim, dims...) where {T}
    if (dim, dims...) == (4, 4)
        one(SeitzOperator{T})
    else
        return Array{T}(undef, dims)
    end
end

# See https://github.com/JuliaArrays/StaticArrays.jl/blob/v1.6.2/src/linalg.jl#L7-L25
@inline +(op::SeitzOperator) = op
@inline +(op₁::SeitzOperator, op₂::SeitzOperator) = op₁ .+ op₂
@inline +(A::AbstractArray, op::SeitzOperator) = A .+ op
@inline +(op::SeitzOperator, A::AbstractArray) = op .+ A

@inline -(op::SeitzOperator) = -1 .* op
@inline -(op₁::SeitzOperator, op₂::SeitzOperator) = op₁ .- op₂
@inline -(A::AbstractArray, op::SeitzOperator) = A .- op
@inline -(op::SeitzOperator, A::AbstractArray) = op .- A

@inline *(n::Number, op::SeitzOperator) = n .* op
@inline *(op::SeitzOperator, n::Number) = n .* op

@inline /(op::SeitzOperator, n::Number) = op ./ n
