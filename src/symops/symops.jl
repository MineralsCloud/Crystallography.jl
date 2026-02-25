using LinearAlgebra: I, diagm, det
using StaticArrays: MMatrix, Size

import CrystallographyBase: shift

export SeitzOperator,
    shift, is_translation, is_point_symmetry, get_translation, get_point_symmetry

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
mutable struct SeitzOperator{T} <: AbstractMatrix{T}
    data::MMatrix{4,4,T,16}
end
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
function SeitzOperator(𝐑::AbstractMatrix, 𝐭::AbstractVector)
    @assert isone(abs(det(𝐑)))
    T = promote_type(eltype(𝐑), eltype(𝐭))
    data = zero(MMatrix{4,4,T,16})
    data[1:3, 1:3] = 𝐑
    data[1:3, 4] = 𝐭
    data[4, 4] = oneunit(T)
    return SeitzOperator(data)
end

"""
    shift(op::SeitzOperator, 𝐱::AbstractVector)

Return a `SeitzOperator` that locates at `𝐱` from a `SeitzOperator` passing through the
origin.
"""
function shift(op::SeitzOperator, 𝐱::AbstractVector)
    @assert length(𝐱) == 3
    op′ = SeitzOperator(𝐱)
    return ConjugacyOperation(op′)(op)
end

(op::SeitzOperator)(𝐫::AbstractVector) = apply(Size(size(𝐫)), op, 𝐫)

function is_translation(op::SeitzOperator)
    if op[1:3, 1:3] != I || !iszero(op[4, 1:3]) || !isone(op[4, 4])
        return false
    end
    return true
end

function is_point_symmetry(op::SeitzOperator)
    if iszero(op[4, 1:3]) &&
        iszero(op[1:3, 4]) &&
        isone(op[4, 4]) &&
        abs(det(op[1:3, 1:3])) == 1
        return true
    else
        false
    end
end

get_translation(op::SeitzOperator) = op[1:3, 4]

get_point_symmetry(op::SeitzOperator) = op[1:3, 1:3]

include("interface.jl")
include("spglib.jl")
include("conjugacy.jl")
