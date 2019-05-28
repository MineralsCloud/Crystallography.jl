"""
# module SeitzOperators



# Examples

```jldoctest
julia>
```
"""
module SeitzOperators

using LinearAlgebra: diagm, I

using CoordinateTransformations: Translation, LinearMap
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
IdentityOperator(T::Type{<: Real}) = IdentityOperator(IMatrix{4,T}())

struct TranslationOperator{T} <: SeitzOperator{T}
    m::T
    function TranslationOperator{T}(m) where {T}
        size(m) == (4, 4) || throw(DimensionMismatch("The operator must be of size 4x4!"))
        new(m)
    end
end
TranslationOperator(m::T) where {T} = TranslationOperator{T}(m)
function TranslationOperator(t::AbstractVector{T}) where {T}
    TranslationOperator(vcat(hcat(diagm(0 => ones(T, 3)), t), lastrow(T)))
end
TranslationOperator(t::Translation) = TranslationOperator(t.translation)

struct PointSymmetryOperator{T} <: SeitzOperator{T}
    m::T
    function PointSymmetryOperator{T}(m) where {T}
        size(m) == (4, 4) || throw(DimensionMismatch("The operator must be of size 4x4!"))
        new(m)
    end
end
PointSymmetryOperator(m::T) where {T} = PointSymmetryOperator{T}(m)
function PointSymmetryOperator(sym::AbstractMatrix{T}) where {T}
    x = vcat(hcat(sym, zeros(T, 3)), lastrow(T))
    PointSymmetryOperator{typeof(x)}(x)
end
PointSymmetryOperator(sym::LinearMap) = PointSymmetryOperator(sym.linear)

lastrow(T::Type{<: Real}) = [zeros(T, 3)... ones(T, 1)]

end
