using LinearAlgebra: Diagonal

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

Base.size(::SeitzOperator) = (4, 4)

Base.getindex(op::SeitzOperator, i::Int) = getindex(parent(op), i)
Base.getindex(op::SeitzOperator, I...) = getindex(parent(op), I...)

Base.setindex!(op::SeitzOperator, v, i::Int) = setindex!(parent(op), v, i)
Base.setindex!(op::SeitzOperator, X, I...) = setindex!(parent(op), X, I...)

Base.IndexStyle(::Type{<:SeitzOperator}) = IndexLinear()

# Customizing broadcasting
# See https://github.com/JuliaArrays/StaticArraysCore.jl/blob/v1.4.2/src/StaticArraysCore.jl#L397-L398
# and https://github.com/JuliaLang/julia/blob/v1.10.0-beta1/stdlib/LinearAlgebra/src/structuredbroadcast.jl#L7-L14
struct SeitzOperatorStyle <: Broadcast.AbstractArrayStyle{2} end
SeitzOperatorStyle(::Val{2}) = SeitzOperatorStyle()
SeitzOperatorStyle(::Val{N}) where {N} = Broadcast.DefaultArrayStyle{N}()

Base.BroadcastStyle(::Type{<:SeitzOperator}) = SeitzOperatorStyle()

Base.similar(::Broadcast.Broadcasted{SeitzOperatorStyle}, ::Type{T}) where {T} =
    similar(SeitzOperator{T}, 4, 4)
# Override https://github.com/JuliaLang/julia/blob/v1.10.0-beta2/base/abstractarray.jl#L839
function Base.similar(op::SeitzOperator, ::Type{T}, dims::Dims) where {T}
    if dims == size(op)
        one(SeitzOperator{T})
    else
        return throw(ArgumentError("invalid dimensions `$dims` for `SeitzOperator`!"))
    end
end
# Override https://github.com/JuliaLang/julia/blob/v1.10.0-beta1/base/abstractarray.jl#L874
function Base.similar(::Type{SeitzOperator{T}}, dims::Dims) where {T}
    if dims == (4, 4)
        one(SeitzOperator{T})
    else
        return throw(ArgumentError("invalid dimensions `$dims` for `SeitzOperator`!"))
    end
end
