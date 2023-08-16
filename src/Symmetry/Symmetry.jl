# module Symmetry

using LinearAlgebra: I, diagm, det
using StaticArrays: MMatrix

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
mutable struct SeitzOperator{T} <: AbstractMatrix{T}
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

(op::SeitzOperator)(𝐫::AbstractVector) = apply(Size(size(𝐫)), op, 𝐫)

function istranslation(op::SeitzOperator)
    if op[1:3, 1:3] != I || !iszero(op[4, 1:3]) || !isone(op[4, 4])
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

include("misc.jl")
include("spglib.jl")

# end
