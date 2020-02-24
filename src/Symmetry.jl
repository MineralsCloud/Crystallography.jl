module Symmetry

using LinearAlgebra: diagm, I

using CoordinateTransformations: AffineMap, Translation, LinearMap
using LibSymspg: get_symmetry
using StaticArrays: StaticMatrix

using Crystallography: CrystalCoordinates, Cell

export SeitzOperator
export getsymmetry, isidentity, istranslation, ispointsymmetry

function getsymmetry(cell::Cell, symprec::AbstractFloat = 1e-5)
    maps, translations = get_symmetry(
        cell.lattice,
        cell.positions,
        cell.numbers,
        length(cell.numbers),
        symprec,
    )
    return (AffineMap(m, t) for (m, t) in zip(maps, translations))
end

struct SeitzOperator{T}
    m::T
end
SeitzOperator() = SeitzOperator(ones(Int, 4, 4))
function SeitzOperator(m::LinearMap)
    @assert size(m) == (3, 3)
    x = diagm(ones(eltype(m.linear), 4))
    x[1:3, 1:3] = m
    return SeitzOperator(x)
end # function PointSymmetryOperator
function SeitzOperator(t::Translation)
    @assert length(t) == 3
    x = diagm(ones(eltype(t.translation), 4))
    x[1:3, 4] = t
    return SeitzOperator(x)
end # function TranslationOperator
SeitzOperator(m::LinearMap, t::Translation) = SeitzOperator(t ∘ m ∘ inv(t))
SeitzOperator(a::AffineMap) = SeitzOperator(a.linear, a.translation)
function SeitzOperator(s::SeitzOperator, pos::AbstractVector)
    @assert length(pos) == 3
    t = SeitzOperator(Translation(pos))
    return t * s * inv(t)
end # function SeitzOperator

isidentity(op::SeitzOperator) = op.m == I

function istranslation(op::SeitzOperator)
    m = op.m
    if m[1:3, 1:3] != I || !(iszero(m[4, 1:3]) && isone(m[4, 4]))
        return false
    end
    return true
end # function istranslation

function ispointsymmetry(op::SeitzOperator)
    m = op.m
    if !(iszero(m[4, 1:3]) && iszero(m[1:3, 4]) && isone(m[4, 4]))
        return false
    end
    return true
end # function ispointsymmetry

Base.:*(m::SeitzOperator, c::CrystalCoordinates) = CrystalCoordinates((m.m*[c; 1])[1:3])
Base.:*(a::SeitzOperator, b::SeitzOperator) = SeitzOperator(a.m * b.m)

end # module Symmetry
