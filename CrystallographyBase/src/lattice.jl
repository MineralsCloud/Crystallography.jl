using EnumX: @enumx
using LinearAlgebra: Diagonal, I

export CrystalSystem, LatticeSystem, Bravais, Lattice
export latticesystem, basis_vectors, cellparameters

@enumx LatticeSystem begin
    Triclinic = 1
    Monoclinic = 2
    Orthorhombic = 3
    Tetragonal = 4
    Rhombohedral = 5
    Hexagonal = 6
    Cubic = 7
end

@enumx CrystalSystem begin
    Triclinic = 1
    Monoclinic = 2
    Orthorhombic = 3
    Tetragonal = 4
    Trigonal = 5
    Hexagonal = 6
    Cubic = 7
end

@enumx BravaisArithmeticClass begin
    PrimitiveTriclinic = 1
    PrimitiveMonoclinic = 2
    BaseCenteredMonoclinic = 3
    PrimitiveOrthorhombic = 4
    BaseCenteredOrthorhombic = 5
    BodyCenteredOrthorhombic = 6
    FaceCenteredOrthorhombic = 7
    PrimitiveTetragonal = 8
    BodyCenteredTetragonal = 9
    PrimitiveHexagonal = 10
    PrimitiveRhombohedral = 11
    RCentredHexagonal = PrimitiveRhombohedral
    PrimitiveCubic = 12
    BodyCenteredCubic = 13
    FaceCenteredCubic = 14
end
const Bravais = BravaisArithmeticClass

"Represent the real lattices and the reciprocal lattices."
abstract type AbstractLattice{T} end

struct Lattice{T} <: AbstractLattice{T}
    data::SMatrix{3,3,T,9}
end
"""
    Lattice(mat::AbstractMatrix)

Construct a `Lattice` from a matrix.

!!! note
    The basis vectors of the matrix are stored as columns.
"""
Lattice(mat::AbstractMatrix) = Lattice(SMatrix{3,3}(mat))
"""
    Lattice(𝐚::AbstractVector, 𝐛::AbstractVector, 𝐜::AbstractVector)

Construct a `Lattice` from three basis vectors.
"""
Lattice(𝐚::AbstractVector, 𝐛::AbstractVector, 𝐜::AbstractVector) = Lattice(hcat(𝐚, 𝐛, 𝐜))
Lattice(vecs::AbstractVector{<:AbstractVector}) = Lattice(reduce(hcat, vecs))
"""
    Lattice(a, b, c, α, β, γ)

Construct a `Lattice` from the six cell parameters.

The convention we used here is that edge vector 𝐚 in the positive x-axis direction,
edge vector 𝐛 in the x-y plane with positive y-axis component,
and edge vector 𝐜 with positive z-axis component in the Cartesian-system.
See [Wikipedia](https://en.wikipedia.org/w/index.php?title=Fractional_coordinates&oldid=961675499#In_crystallography).
"""
function Lattice(a, b, c, α, β, γ)
    Ω = cellvolume(a, b, c, α, β, γ)
    sinγ, cosγ, cosα, cosβ = sind(γ), cosd(γ), cosd(α), cosd(β)
    return Lattice(
        [a, zero(a), zero(a)],
        [b * cosγ, b * sinγ, zero(b)],
        [c * cosβ, c * (cosα - cosβ * cosγ) / sinγ, Ω / (a * b * sinγ)],
    )
end
@functor Lattice

# See https://github.com/korsbo/Latexify.jl/blob/5859690/src/Latexify.jl#L19-L27
const ANGLE_TOLERANCE = 1e-5
function angle_tolerance(v)
    global ANGLE_TOLERANCE = v
end

const LENGTH_TOLERANCE = 1e-5
function length_tolerance(v)
    global LENGTH_TOLERANCE = v
end

"""
    basis_vectors(lattice::Lattice)

Get the three basis vectors from a `lattice`.
"""
basis_vectors(lattice::Lattice) = lattice[:, 1], lattice[:, 2], lattice[:, 3]

"""
    latticesystem(bravais::Bravais)

Get the crystal system of a Bravais type.
"""
function latticesystem(bravais::Bravais.T)
    index = Int(bravais)
    if index == 1
        return LatticeSystem.Triclinic
    elseif 2 <= index <= 3
        return LatticeSystem.Monoclinic
    elseif 4 <= index <= 7
        return LatticeSystem.Orthorhombic
    elseif 8 <= index <= 9
        return LatticeSystem.Tetragonal
    elseif index == 10
        return LatticeSystem.Hexagonal
    elseif index == 11
        return LatticeSystem.Rhombohedral
    else  # 12 <= index <= 14
        return LatticeSystem.Cubic
    end
end
"""
    latticesystem(a, b, c, α, β, γ)

Guess the crystal system from the six cell parameters.
"""
# See https://github.com/LaurentRDC/crystals/blob/2d3a570/crystals/lattice.py#L396-L475
function latticesystem(a, b, c, α, β, γ)
    lengths, angles = Base.vect(a, b, c), Base.vect(α, β, γ)
    unilength = all(x ≅ a for x in lengths)
    uniangle = all(θ ≊ α for θ in angles)
    # Checking for monoclinic system is generalized
    # to the case where a, b, and c can be cycled,
    # i.e., a != c && β != 90 && α == γ == 90
    #   || b != c && α != 90 && β == γ == 90
    #   || a != b && γ != 90 && α == β != 90
    for (lengths′, angles′) in zip(cyclic_perm(lengths), cyclic_perm(angles))
        (a′, _, c′), (α′, β′, γ′) = lengths′, angles′
        if !(a′ ≅ c′) && α′ ≊ 90 && γ′ ≊ 90 && !(β′ ≊ 90)
            return LatticeSystem.Monoclinic
        end
    end
    if unilength && uniangle
        return α ≊ 90 ? LatticeSystem.Cubic : LatticeSystem.Rhombohedral
    end
    if unilength && !uniangle  # Technically, a hexagonal system could have all 3 lengths equal.
        if any(θ ≊ 120 for θ in angles) && sum(θ ≊ 90 for θ in angles) == 2
            return LatticeSystem.Hexagonal
        end
    end
    if bilengths(lengths)  # At this point, two lengths are equal at most.
        if uniangle && α ≊ 90
            return LatticeSystem.Tetragonal
        elseif any(θ ≊ 120 for θ in angles) && sum(θ ≊ 90 for θ in angles) == 2
            return LatticeSystem.Hexagonal
        end
    else  # At this point, all lengths are not equal.
        return uniangle && α ≊ 90 ? LatticeSystem.Orthorhombic : LatticeSystem.Triclinic
    end
end
"""
    crystalsystem(lattice::Lattice)

Get the crystal system of a `lattice`.
"""
latticesystem(lattice::Lattice; kwargs...) =
    latticesystem(cellparameters(lattice)...; kwargs...)
# Auxiliary functions
≊(θ, φ) = isapprox(θ, φ; rtol = ANGLE_TOLERANCE)
≅(x, y) = isapprox(x, y; rtol = LENGTH_TOLERANCE)
cyclic_perm(vec) = (circshift(vec, i) for i in 1:length(vec))  # See https://stackoverflow.com/a/43035441
function bilengths(iterable)  # If and only if two lengths are equal.
    for i in iterable
        if sum(isapprox(i, j, rtol = LENGTH_TOLERANCE) for j in iterable) == 2
            return true
        end
    end
    return false
end

"""
    cellparameters(lattice::Lattice)

Get the six cell parameters from a `lattice`.
"""
function cellparameters(lattice::Lattice)
    𝐚, 𝐛, 𝐜 = basis_vectors(lattice)
    a, b, c = norm(𝐚), norm(𝐛), norm(𝐜)
    γ, β, α =
        acosd(dot(𝐚, 𝐛) / (a * b)), acosd(dot(𝐚, 𝐜) / (a * c)), acosd(dot(𝐛, 𝐜) / (b * c))
    return a, b, c, α, β, γ
end

# See https://en.wikipedia.org/wiki/Supercell_(crystal)
"""
    supercell(lattice::Lattice, scaling_factors::AbstractMatrix{<:Integer})
    supercell(lattice::Lattice, scaling_factors::AbstractVector{<:Integer})
    supercell(lattice::Lattice, scaling_factor::Integer)

Create a supercell from `lattice`.
"""
function supercell(lattice::Lattice, scaling_factors::AbstractMatrix{<:Integer})
    if size(scaling_factors) != (3, 3)
        throw(ArgumentError("`scaling_factors` must be a 3×3 matrix!"))
    end
    @assert det(scaling_factors) >= 1
    return Lattice(lattice.data * scaling_factors)
end
supercell(lattice::Lattice, scaling_factors::AbstractVector{<:Integer}) =
    supercell(lattice, Diagonal(scaling_factors))
supercell(lattice::Lattice, scaling_factor::Integer) =
    supercell(lattice, scaling_factor * I)

Base.iterate(lattice::AbstractLattice) = iterate(lattice.data)
Base.iterate(lattice::AbstractLattice, state) = iterate(lattice.data, state)

Base.eltype(::AbstractLattice{T}) where {T} = T

Base.length(::AbstractLattice) = 9

Base.size(::AbstractLattice) = (3, 3)
Base.size(::AbstractLattice, dim::Integer) = dim <= 2 ? 3 : 1

Base.IteratorSize(::Type{<:AbstractLattice}) = Base.HasShape{2}()

Base.axes(lattice::AbstractLattice, dim::Integer) = axes(lattice.data, dim)

Base.getindex(lattice::AbstractLattice, i) = getindex(lattice.data, i)
Base.getindex(lattice::AbstractLattice, I::Vararg) = getindex(lattice.data, I...)

Base.firstindex(::AbstractLattice) = 1

Base.lastindex(::AbstractLattice) = 9

for op in (:+, :-)
    @eval Base.broadcast(::typeof($op), lattice::AbstractLattice, number::Number) =
        Lattice(broadcast($op, lattice.data, number))
    @eval Base.broadcast(::typeof($op), number::Number, lattice::AbstractLattice) =
        broadcast($op, lattice, number)
end
for op in (:*, :/, ://)
    @eval Base.$op(lattice::AbstractLattice, number::Number) =
        Lattice(($op)(lattice.data, number))
    @eval Base.$op(number::Number, lattice::AbstractLattice) = ($op)(lattice, number)
end

function Base.show(io::IO, x::AbstractLattice)
    if get(io, :compact, false) || get(io, :typeinfo, nothing) == typeof(x)
        Base.show_default(IOContext(io, :limit => true), x)  # From https://github.com/mauro3/Parameters.jl/blob/ecbf8df/src/Parameters.jl#L556
    else
        println(io, string(typeof(x)))
        for row in eachrow(x.data)
            println(io, ' ', join(row, "  "))
        end
    end
end
