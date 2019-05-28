"""
# module SeitzOperators



# Examples

```jldoctest
julia>
```
"""
module SeitzOperators

using LinearAlgebra: diagm, I

using LuxurySparse: IMatrix

export SeitzOperator,
    IdentityOperator,
    TranslationOperator,
    PointSymmetryOperator

abstract type SeitzOperator{T <: AbstractMatrix} end

struct IdentityOperator{T} <: SeitzOperator{T}
    m::T
    function IdentityOperator{T}(m) where {T}
        size(m) == (4, 4) || throw(DimensionMismatch("The operator must be of size 4x4!"))
        m == I || error("The matrix is not an identity matrix!")
        new(m)
    end
end
IdentityOperator(m::T) where {T} = IdentityOperator{T}(m)

struct TranslationOperator{T} <: SeitzOperator{T}
    m::T
    function TranslationOperator{T}(m) where {T}
        size(m) == (4, 4) || throw(DimensionMismatch("The operator must be of size 4x4!"))
        new(m)
    end
end
TranslationOperator(m::T) where {T} = TranslationOperator{T}(m)
function TranslationOperator(vec::AbstractVector{T}) where {T}
    TranslationOperator(vcat(hcat(diagm(0 => ones(T, 3)), vec), lastrow(T)))
end

struct PointSymmetryOperator{T} <: SeitzOperator{T}
    m::T
    function PointSymmetryOperator{T}(m) where {T}
        size(m) == (4, 4) || throw(DimensionMismatch("The operator must be of size 4x4!"))
        new(m)
    end
end
PointSymmetryOperator(m::T) where {T} = PointSymmetryOperator{T}(m)
function PointSymmetryOperator(mat::AbstractMatrix{T}) where {T}
    PointSymmetryOperator(vcat(hcat(mat, zeros(T, 3)), lastrow(T)))
end

lastrow(::Type{T}) where {T} = [zeros(T, 3)... ones(T, 1)]

end
