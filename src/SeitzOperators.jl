"""
# module SeitzOperators



# Examples

```jldoctest
julia>
```
"""
module SeitzOperators

using LinearAlgebra

export SeitzOperator,
    IdentityOperator,
    TranslationOperator,
    PointSymmetryOperator

abstract type SeitzOperator{T} end

struct IdentityOperator{T} <: SeitzOperator{T}
    m::T
end
IdentityOperator() = IdentityOperator(ones(4))

struct TranslationOperator{T} <: SeitzOperator{T}
    m::T
    function TranslationOperator(vec::AbstractVector{T}) where {T}
        vcat(hcat(diagm(0 => ones(T, 3)), vec), lastrow(T))
    end
end

struct PointSymmetryOperator{T} <: SeitzOperator{T}
    m::T
    function PointSymmetryOperator(mat::AbstractMatrix{T}) where {T}
        vcat(hcat(mat, zeros(T, 3)), lastrow(T))
    end
end

lastrow(::Type{T}) where {T} = [zeros(T, 3)... ones(T, 1)]

end
