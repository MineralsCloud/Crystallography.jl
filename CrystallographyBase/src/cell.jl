using StaticArrays: MMatrix, MVector
using StructHelpers: @batteries

export Cell

struct Cell{L,P,T}
    lattice::Lattice{L}
    positions::Vector{MVector{3,P}}
    atoms::Vector{T}
end
function Cell(lattice, positions, atoms)
    if !(lattice isa Lattice)
        lattice = Lattice(lattice)
        L = eltype(lattice)
    end
    if positions isa AbstractVector
        P = eltype(Base.promote_typeof(positions...))
        positions = collect(map(MVector{3,P}, positions))
    else
        throw(ArgumentError("`positions` must be a `Vector` of `Vector`s!"))
    end
    T = eltype(atoms)
    return Cell{L,P,T}(lattice, positions, atoms)
end

@batteries Cell eq = true hash = true

"""
    Lattice(cell::Cell)

Get the lattice of a `Cell`.
"""
Lattice(cell::Cell) = cell.lattice

"""
basis_vectors(cell::Cell)
Return the three basis vectors from `cell`.
"""
function basis_vectors(cell::Cell)
    lattice = cell.lattice
    return lattice[:, 1], lattice[:, 2], lattice[:, 3]
end
