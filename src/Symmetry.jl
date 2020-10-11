module Symmetry

using CoordinateTransformations: AffineMap, Translation, LinearMap
using LibSymspg: get_symmetry, get_spacegroup, ir_reciprocal_mesh
using LinearAlgebra: I, diagm, det, tr
using StaticArrays: SVector, SMatrix, SDiagonal

using Crystallography

import LinearAlgebra

export SeitzOperator,
    symmetrytype,
    getsymmetry,
    getspacegroup,
    irreciprocalmesh,
    isidentity,
    istranslation,
    ispointsymmetry,
    genpath,
    pearsonsymbol,
    arithmeticclass

pearsonsymbol(::Triclinic) = "a"
pearsonsymbol(::Monoclinic) = "m"
pearsonsymbol(::Orthorhombic) = "o"
pearsonsymbol(::Tetragonal) = "t"
pearsonsymbol(::Cubic) = "c"
pearsonsymbol(::Hexagonal) = "h"
pearsonsymbol(::Trigonal) = "h"
pearsonsymbol(::Primitive) = "P"
pearsonsymbol(::BaseCentering{T}) where {T} = string(T)
pearsonsymbol(::BodyCentering) = "I"
pearsonsymbol(::FaceCentering) = "F"
pearsonsymbol(::RhombohedralCentering) = "R"
pearsonsymbol(b::Bravais) = pearsonsymbol(crystalsystem(b)) * pearsonsymbol(centering(b))

arithmeticclass(::Triclinic) = "-1"
arithmeticclass(::Monoclinic) = "2/m"
arithmeticclass(::Orthorhombic) = "mmm"
arithmeticclass(::Tetragonal) = "4/mmm"
arithmeticclass(::Hexagonal) = "6/mmm"
arithmeticclass(::Cubic) = "m-3m"
arithmeticclass(::Primitive) = "P"
arithmeticclass(::BaseCentering) = "S"
arithmeticclass(::BodyCentering) = "I"
arithmeticclass(::FaceCentering) = "F"
arithmeticclass(::RhombohedralCentering) = "R"
arithmeticclass(b::Bravais) =
    arithmeticclass(crystalsystem(b)) * arithmeticclass(centering(b))
arithmeticclass(::RCenteredHexagonal) = "-3mR"

abstract type PointSymmetry end
struct Identity <: PointSymmetry end
struct RotationAxis{N} <: PointSymmetry end
struct Inversion <: PointSymmetry end
struct RotoInversion{N} <: PointSymmetry end
const Mirror = RotoInversion{2}
RotationAxis(N::Int) =
    N ∈ (2, 3, 4, 6) ? RotationAxis{N}() :
    throw(ArgumentError("rotation axis must be either 2, 3, 4 or 6!"))
RotoInversion(N::Int) =
    N ∈ (2, 3, 4, 6) ? RotoInversion{N}() :
    throw(ArgumentError("rotoinversion axis must be either 2, 3, 4 or 6!"))

struct PointSymmetryPower{T<:PointSymmetry,N} end
PointSymmetryPower(X::PointSymmetry, N::Int) = PointSymmetryPower{typeof(X),N}()
function PointSymmetryPower(::RotationAxis{N}, M::Int) where {N}
    if M >= N
        if M == N
            Identity()
        else
            PointSymmetryPower(RotationAxis(N), M - N)  # Recursive call
        end
    else  # Until M < N
        PointSymmetryPower{RotationAxis{N},M}()
    end
end # function PointSymmetryPower

function symmetrytype(trace, determinant)
    return Dict(
        (3, 1) => Identity(),
        (-1, 1) => RotationAxis(2),
        (0, 1) => RotationAxis(3),
        (1, 1) => RotationAxis(4),
        (2, 1) => RotationAxis(6),
        (-3, -1) => Inversion(),
        (1, -1) => Mirror(),
        (0, -1) => RotoInversion(3),
        (-1, -1) => RotoInversion(4),
        (-2, -1) => RotoInversion(6),
    )[(trace, determinant)]
end # function symmetrytype
symmetrytype(op::PointSymmetry) = symmetrytype(tr(op), det(op))

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
(op::SeitzOperator)(v::AbstractVector) = (op.data*[v; 1])[1:3]

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
        iszero(m[1:3, 4]) &&
        isone(m[4, 4]) &&
        abs(det(m[1:3, 1:3])) == 1
    )
        return false
    end
    return true
end # function ispointsymmetry

"""
    genpath(nodes, densities)
    genpath(nodes, density::Integer, iscircular = false)

Generate a reciprocal space path from each node.

# Arguments
- `nodes::AbstractVector`: a vector of 3-element k-points.
- `densities::AbstractVector{<:Integer}`: number of segments between current node and the next node.
- `density::Integer`: assuming constant density between nodes.

# Examples
```jldoctest
julia> nodes = [
    [0.0, 0.0, 0.5],
    [0.0, 0.0, 0.0],
    [0.3333333333, 0.3333333333, 0.5],
    [0.3333333333, 0.3333333333, -0.5],
    [0.3333333333, 0.3333333333, 0.0],
    [0.5, 0.0, 0.5],
    [0.5, 0.0, 0.0]
];

julia> genpath(nodes, 100, true)  # Generate a circular path
700-element Array{Array{Float64,1},1}:
...

julia> genpath(nodes, 100, false)  # Generate a noncircular path
600-element Array{Array{Float64,1},1}:
...

julia> genpath(nodes, [10 * i for i in 1:6])  # Generate a noncircular path
210-element Array{Array{Float64,1},1}:
...
```
"""
function genpath(nodes, densities)
    if length(densities) == length(nodes) || length(densities) == length(nodes) - 1
        path = similar(nodes, sum(densities .- 1) + length(densities))
        s = 0
        for (i, (thisnode, nextnode, density)) in
            enumerate(zip(nodes, circshift(nodes, -1), densities))
            step = @. (nextnode - thisnode) / density
            for j = 1:density
                path[s+j] = @. thisnode + j * step
            end
            s += density
        end
        return path
    else
        throw(DimensionMismatch("the length of `densities` is either `length(nodes)` or `length(nodes) - 1`!"))
    end
end # function genpath
genpath(nodes, density::Integer, iscircular::Bool = false) =
    genpath(nodes, (density for _ = 1:(length(nodes)-(iscircular ? 0 : 1))))

Base.getindex(A::SeitzOperator, I::Vararg{Int}) = getindex(A.data, I...)

Base.one(::Type{SeitzOperator{T}}) where {T} =
    SeitzOperator(SDiagonal(SVector{4}(ones(T, 4))))
Base.one(A::SeitzOperator) = one(typeof(A))

Base.inv(op::SeitzOperator) = SeitzOperator(Base.inv(op.data))

Base.:*(::Identity, ::Identity) = Identity()
Base.:*(::Identity, x::Union{PointSymmetryPower,PointSymmetry}) = x
Base.:*(x::Union{PointSymmetryPower,PointSymmetry}, ::Identity) = x
Base.:*(::Inversion, ::Inversion) = Identity()
Base.:*(::RotationAxis{N}, ::Inversion) where {N} = RotoInversion(N)
Base.:*(i::Inversion, r::RotationAxis) = r * i
Base.:*(::RotationAxis{N}, ::RotationAxis{N}) where {N} =
    PointSymmetryPower(RotationAxis(N), 2)
Base.:*(::RotationAxis{2}, ::RotationAxis{2}) = Identity()
Base.:*(::PointSymmetryPower{RotationAxis{N},M}, ::RotationAxis{N}) where {N,M} =
    PointSymmetryPower(RotationAxis(N), M + 1)
Base.:∘(a::SeitzOperator, b::SeitzOperator) = SeitzOperator(a.data * b.data)

function Base.convert(::Type{Translation}, op::SeitzOperator)
    @assert(istranslation(op), "operator is not a translation!")
    return Translation(collect(op.data[1:3, 4]))
end # function Base.convert
function Base.convert(::Type{LinearMap}, op::SeitzOperator)
    @assert(ispointsymmetry(op), "operator is not a point symmetry!")
    return LinearMap(collect(op.data[1:3, 1:3]))
end # function Base.convert

LinearAlgebra.tr(::Identity) = 3
LinearAlgebra.tr(::RotationAxis) where {N} = N - 3  # 2, 3, 4
LinearAlgebra.tr(::RotationAxis{6}) = 2
LinearAlgebra.tr(::Inversion) = -3
LinearAlgebra.tr(::RotoInversion{N}) where {N} = -tr(RotationAxis(N))
LinearAlgebra.det(::Union{Identity,RotationAxis}) = 1
LinearAlgebra.det(::Inversion) = -1
LinearAlgebra.det(::RotoInversion) = -1

end # module Symmetry
