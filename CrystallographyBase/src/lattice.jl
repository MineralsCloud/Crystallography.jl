using LinearAlgebra: Diagonal, I

export Lattice, latticevectors, latticesystem, latticeconstants, periodicity, supercell

"Represent the real lattices and the reciprocal lattices."
abstract type AbstractLattice{T} end
@struct_hash_equal_isequal_isapprox struct Lattice{T} <: AbstractLattice{T}
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
"""
    Lattice(latticevectors::AbstractVector{<:AbstractVector})

Construct a `Lattice` from a vector of three basis vectors.
"""
Lattice(latticevectors::AbstractVector{<:AbstractVector}) =
    Lattice(reduce(hcat, latticevectors))
"""
    Lattice(a, b, c, α, β, γ; axis = :a)

Construct a `Lattice` from the six cell parameters.

The default convention we used here is that edge vector 𝐚 in the positive x-axis direction,
edge vector 𝐛 in the x-y plane with a positive y-axis component,
and edge vector 𝐜 with a positive z-axis component in the Cartesian system.
See [Wikipedia](https://en.wikipedia.org/w/index.php?title=Fractional_coordinates&oldid=961675499#In_crystallography).
You can also choose `axis = :c`.
"""
function Lattice(a, b, c, α, β, γ; axis=:a)
    Ω = cellvolume(a, b, c, α, β, γ)
    if axis == :a  # See https://en.wikipedia.org/w/index.php?title=Fractional_coordinates&oldid=961675499#In_crystallography
        sinγ, cosγ, cosα, cosβ, 𝟎 = sind(γ), cosd(γ), cosd(α), cosd(β), zero(a)
        return Lattice(
            [a, 𝟎, 𝟎],
            [b * cosγ, b * sinγ, zero(b)],
            [c * cosβ, c * (cosα - cosβ * cosγ) / sinγ, Ω / (a * b * sinγ)],
        )
    elseif axis == :c  # See https://github.com/LaurentRDC/crystals/blob/2d3a570/crystals/lattice.py#L356-L391
        sinα, cosα, sinβ, cosβ, 𝟎 = sind(α), cosd(α), sind(β), cosd(β), zero(c)
        x = Ω / (b * c * sinα)
        cos′ = (cosα * cosβ - cosd(γ)) / (sinα * sinβ)
        sin′ = sqrt(1 - cos′^2)
        return Lattice(
            [x, -x * cos′ / sin′, a * cosβ], [zero(b), b * sinα, b * cosα], [𝟎, 𝟎, c]
        )
    else
        error("aligning `$axis` axis is not supported!")
    end
end
@functor Lattice

"""
    latticevectors(lattice::Lattice)

Get the three primitive vectors from a `lattice`.
"""
latticevectors(lattice::Lattice) = lattice[:, 1], lattice[:, 2], lattice[:, 3]

"""
    latticesystem(bravais::Bravais)

Get the lattice system of a Bravais lattice.
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
    latticesystem(a, b, c, α, β, γ; kwargs...)

Guess the lattice system from the six lattice constants.

# Arguments
- `angletol=1e-5`: the absolute tolerance of angles (`α`, `β`, `γ`).
- `lengthtol=1e-5`: the absolute tolerance of edges (`a`, `b`, `c`).
"""
# See https://github.com/LaurentRDC/crystals/blob/2d3a570/crystals/lattice.py#L396-L475
function latticesystem(a, b, c, α, β, γ; angletol=1e-5, lengthtol=1e-5)
    lengths, angles = Base.vect(a, b, c), Base.vect(α, β, γ)
    ≊(θ, φ) = isapprox(θ, φ; atol=angletol)
    ≅(x, y) = isapprox(x, y; atol=lengthtol)
    unilength = all(x ≅ a for x in lengths)
    uniangle = all(θ ≊ α for θ in angles)
    function bilengths(vec)  # If and only if two lengths are equal.
        for x in vec
            if sum(isapprox(x, y; atol=lengthtol) for y in vec) == 2
                return true
            end
        end
        return false
    end
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
    latticesystem(lattice::Lattice; angletol=1e-5, lengthtol=1e-5)

Get the lattice system of a `Lattice`.
"""
latticesystem(lattice::Lattice; kwargs...) =
    latticesystem(latticeconstants(lattice)...; kwargs...)
# Auxiliary functions
cyclic_perm(vec) = (circshift(vec, i) for i in 1:length(vec))  # See https://stackoverflow.com/a/43035441

"""
    latticeconstants(lattice::Lattice)

Get the six lattice constants from a `lattice`.
"""
function latticeconstants(lattice::AbstractLattice)  # Works for `ReciprocalLattice`s, too
    𝐚, 𝐛, 𝐜 = latticevectors(lattice)
    a, b, c = norm(𝐚), norm(𝐛), norm(𝐜)
    γ, β, α = acosd(dot(𝐚, 𝐛) / (a * b)),
    acosd(dot(𝐚, 𝐜) / (a * c)),
    acosd(dot(𝐛, 𝐜) / (b * c))
    return a, b, c, α, β, γ
end

# See https://github.com/LaurentRDC/crystals/blob/2d3a570/crystals/lattice.py#L161-L176
# Add the absolute value of the component of every lattice vector
# along the three euclidian vectors, which is effectively the sum of
# absolutes of rows.
"""
    periodicity(lattice::Lattice)

Get crystal periodicity in ``x``, ``y``, and ``z`` direction from the `Lattice`.
"""
periodicity(lattice::Lattice) = Tuple(sum(abs, lattice.data; dims=2))

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
supercell(lattice_or_cell, scaling_factors::AbstractVector{<:Integer}) =
    supercell(lattice_or_cell, Diagonal(scaling_factors))
# See https://stackoverflow.com/a/57270841
supercell(lattice_or_cell, scaling_factor::Integer) =
    supercell(lattice_or_cell, Matrix(scaling_factor * I, 3, 3))

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
        println(io, " (a, b, c, α, β, γ) = ", latticeconstants(x))
    end
end
