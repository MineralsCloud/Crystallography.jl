module Symmetry

using LinearAlgebra: I, diagm, det

using CoordinateTransformations: AffineMap, Translation, LinearMap
using LibSymspg: get_symmetry, get_spacegroup, ir_reciprocal_mesh
using StaticArrays: SVector, SMatrix, SDiagonal

using Crystallography: Crystal, Cell

export SeitzOperator
export getsymmetry,
    getspacegroup, irreciprocalmesh, isidentity, istranslation, ispointsymmetry

# These are helper methods and should not be exported!
_numbers(a::AbstractVector{<:Integer}) = a
function _numbers(a::AbstractVector{T}) where {T}
    d = Dict(value => key for (key, value) in pairs(unique(a)))
    return [d[v] for v in a]
end

function getsymmetry(cell::Cell, symprec::AbstractFloat = 1e-5; seitz::Bool = false)
    maps, translations = get_symmetry(
        cell.lattice,
        cell.positions,
        _numbers(cell.atoms),
        length(cell.atoms),
        symprec,
    )
    return if seitz
        (SeitzOperator(LinearMap(m), Translation(t)) for (m, t) in zip(maps, translations))
    else
        (AffineMap(m, t) for (m, t) in zip(maps, translations))
    end
end

function getspacegroup(cell::Cell, symprec::AbstractFloat = 1e-5)
    return get_spacegroup(
        cell.lattice,
        cell.positions,
        _numbers(cell.atoms),
        length(cell.atoms),
        symprec,
    )
end # function getspacegroup

function irreciprocalmesh(
    cell::Cell,
    mesh::AbstractVector{Int},
    symprec::AbstractFloat = 1e-5;
    is_shift::AbstractVector{Bool} = falses(3),
    is_time_reversal::Bool = false,
)
    return ir_reciprocal_mesh(
        mesh,
        collect(is_shift),
        is_time_reversal,
        cell.lattice,
        cell.positions,
        _numbers(cell.atoms),
        length(cell.atoms),
        symprec,
    )
end # function irreciprocalmesh

"""
    SeitzOperator(m::AbstractMatrix)

Construct a `SeitzOperator` from a 4×4 matrix.
"""
struct SeitzOperator{T}
    data::SMatrix{4,4,T}
end
SeitzOperator(m::AbstractMatrix) = SeitzOperator(SMatrix{4,4}(m))
function SeitzOperator(l::LinearMap)
    m = l.linear
    @assert size(m) == (3, 3)
    x = diagm(ones(eltype(m), 4))
    x[1:3, 1:3] = m
    return SeitzOperator(x)
end # function PointSymmetryOperator
function SeitzOperator(t::Translation)
    v = t.translation
    @assert length(v) == 3
    x = diagm(ones(eltype(v), 4))
    x[1:3, 4] = v
    return SeitzOperator(x)
end # function TranslationOperator
"""
    SeitzOperator(l::LinearMap, t::Translation)
    SeitzOperator(a::AffineMap)
    SeitzOperator(l::LinearMap)
    SeitzOperator(t::Translation)

Construct a `SeitzOperator` from `LinearMap`s, `Translation`s or `AffineMap`s by the
following equation:

```math
S = \\left(
\\begin{array}{ccc|c}
& & & \\\\
& R & & {\\boldsymbol \\tau} \\\\
& & & \\\\
\\hline 0 & 0 & 0 & 1
\\end{array}
\\right),
```

where ``R`` is a `LinearMap` and ``{\\boldsymbol \\tau}`` is a `Translation`.
"""
function SeitzOperator(l::LinearMap, t::Translation)
    m, v = l.linear, t.translation
    return SeitzOperator(vcat(hcat(m, v), [zeros(eltype(m), 3)... one(eltype(v))]))
end # function SeitzOperator
SeitzOperator(a::AffineMap) = SeitzOperator(LinearMap(a.linear), Translation(a.translation))
"""
    SeitzOperator(s::SeitzOperator, pos::AbstractVector)

Construct a `SeitzOperator` that locates at `pos` from a `SeitzOperator` passing through the
origin.
"""
function SeitzOperator(s::SeitzOperator, pos::AbstractVector)
    @assert length(pos) == 3
    t = SeitzOperator(Translation(pos))
    return t * s * inv(t)
end # function SeitzOperator

isidentity(op::SeitzOperator) = op.data == I

function istranslation(op::SeitzOperator)
    m = op.data
    if m[1:3, 1:3] != I || !(iszero(m[4, 1:3]) && isone(m[4, 4]))
        return false
    end
    return true
end # function istranslation

function ispointsymmetry(op::SeitzOperator)
    m = op.data
    if !(
        iszero(m[4, 1:3]) &&
        iszero(m[1:3, 4]) && isone(m[4, 4]) && abs(det(m[1:3, 1:3])) == 1
    )
        return false
    end
    return true
end # function ispointsymmetry

Base.getindex(A::SeitzOperator, I::Vararg{Int}) = getindex(A.data, I...)

Base.one(::Type{SeitzOperator{T}}) where {T} =
    SeitzOperator(SDiagonal(SVector{4}(ones(T, 4))))
Base.one(A::SeitzOperator) = one(typeof(A))

Base.inv(op::SeitzOperator) = SeitzOperator(Base.inv(op.data))

Base.:*(m::SeitzOperator, c::Crystal) = Crystal((m.data*[c; 1])[1:3])
Base.:*(a::SeitzOperator, b::SeitzOperator) = SeitzOperator(a.data * b.data)

function Base.convert(::Type{Translation}, op::SeitzOperator)
    @assert(istranslation(op), "operator is not a translation!")
    return Translation(collect(op.data[1:3, 4]))
end # function Base.convert
function Base.convert(::Type{LinearMap}, op::SeitzOperator)
    @assert(ispointsymmetry(op), "operator is not a point symmetry!")
    return LinearMap(collect(op.data[1:3, 1:3]))
end # function Base.convert

end # module Symmetry
