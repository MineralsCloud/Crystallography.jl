using StaticArrays: MMatrix, MVector
using StructHelpers: @batteries

export Cell

struct Cell{L,P,T}
    lattice::MMatrix{3,3,L,9}
    positions::Vector{MVector{3,P}}
    types::Vector{T}
end
function Cell(lattice, positions, types)
    if !(lattice isa AbstractMatrix)
        lattice = reduce(hcat, lattice)  # Use `reduce` can make it type stable
    end
    P = eltype(Base.promote_typeof(positions...))
    positions = collect(map(MVector{3,P}, positions))
    L, T = eltype(lattice), eltype(types)
    return Cell{L,P,T}(lattice, positions, types)
end
Cell(lattice::Lattice, positions, types, magmoms = zeros(length(types))) =
    Cell(lattice.data, positions, types, magmoms)
@batteries Cell eq = true hash = true

"""
    Lattice(cell::Cell)

Get the lattice of a `Cell`.
"""
Lattice(cell::Cell) = Lattice(cell.lattice)

"""
basis_vectors(cell::Cell)
Return the three basis vectors from `cell`.
"""
function basis_vectors(cell::Cell)
    lattice = cell.lattice
    return lattice[:, 1], lattice[:, 2], lattice[:, 3]
end
