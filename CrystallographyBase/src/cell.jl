using LinearAlgebra: isdiag, diag
using StaticArrays: MVector

export Cell, natoms, atomtypes, eachatom

"""
    Cell(lattice, positions, atoms)

Create a new cell.

Argument `lattice` is a [`Lattice`](@ref) type.
Fractional atomic positions `positions` are given
by a vector of ``N`` vectors with floating point values, where ``N`` is the number of atoms.
Argument `atoms` is a list of ``N`` values, where the same kind of atoms
need to be the same type.
"""
@struct_hash_equal_isequal_isapprox struct Cell{L,P,T}
    lattice::Lattice{L}
    positions::Vector{MVector{3,P}}
    atoms::Vector{T}
end
function Cell(lattice, positions, atoms)
    if !(lattice isa Lattice)
        lattice = Lattice(lattice)
    end
    if positions isa AbstractVector
        P = eltype(Base.promote_typeof(positions...))
        positions = collect(map(MVector{3,P}, positions))
    else
        throw(ArgumentError("`positions` must be a `Vector` of `Vector`s!"))
    end
    L, T = eltype(lattice), eltype(atoms)
    return Cell{L,P,T}(lattice, positions, atoms)
end

"""
    supercell(cell::Cell, scaling_factors::AbstractMatrix{<:Integer})
    supercell(cell::Cell, scaling_factors::AbstractVector{<:Integer})
    supercell(cell::Cell, scaling_factor::Integer)

Create a supercell from `cell`.

!!! note
    Currently, only integral replications are supported.
"""
function supercell(cell::Cell, scaling_factors::AbstractMatrix{<:Integer})
    if size(scaling_factors) != (3, 3)
        throw(ArgumentError("`scaling_factors` must be a 3×3 matrix!"))
    end
    @assert isdiag(scaling_factors) "currently not supported!"
    @assert det(scaling_factors) >= 1
    new_atoms = eltype(cell.atoms)[]
    new_positions = eltype(cell.positions)[]
    l, m, n = diag(scaling_factors)
    𝐚, 𝐛, 𝐜 = eachcol(Matrix(I, 3, 3))
    for (atom, position) in eachatom(cell)
        for (i, j, k) in Iterators.product(0:(l - 1), 0:(m - 1), 0:(n - 1))
            push!(new_atoms, atom)
            new_position = position + i * 𝐚 + j * 𝐛 + k * 𝐜
            new_position ./= (l, m, n)  # Make them within the boundary of the cell
            push!(new_positions, new_position)
        end
    end
    new_lattice = supercell(cell.lattice, scaling_factors)
    return Cell(new_lattice, new_positions, new_atoms)
end

natoms(cell::Cell) = length(cell.atoms)

atomtypes(cell::Cell) = unique(cell.atoms)

"""
    Lattice(cell::Cell)

Get the lattice of a `Cell`.
"""
Lattice(cell::Cell) = cell.lattice

struct EachAtom{A,B}
    atoms::Vector{A}
    positions::Vector{B}
end
EachAtom(cell::Cell) = EachAtom(cell.atoms, cell.positions)

eachatom(cell::Cell) = EachAtom(cell)

# Similar to https://github.com/JuliaCollections/IterTools.jl/blob/0ecaa88/src/IterTools.jl#L1028-L1032
function Base.iterate(iter::EachAtom, state = 1)
    if state > length(iter)
        return nothing
    else
        return (iter.atoms[state], iter.positions[state]), state + 1
    end
end

Base.eltype(::EachAtom{A,B}) where {A,B} = Tuple{A,B}

Base.length(iter::EachAtom) = length(iter.atoms)

Base.IteratorSize(::Type{<:EachAtom}) = Base.HasLength()

function Base.show(io::IO, cell::Cell)
    if get(io, :compact, false) || get(io, :typeinfo, nothing) == typeof(cell)
        Base.show_default(IOContext(io, :limit => true), cell)  # From https://github.com/mauro3/Parameters.jl/blob/ecbf8df/src/Parameters.jl#L556
    else
        println(io, string(typeof(cell)))
        println(io, " lattice:")
        for row in eachrow(cell.lattice.data)
            println(io, "  ", join(row, "  "))
        end
        N = natoms(cell)
        println(io, " $N atomic positions:")
        for pos in cell.positions
            println(io, "  ", pos)
        end
        println(io, " $N atoms:")
        println(io, "  ", cell.atoms)
    end
end
