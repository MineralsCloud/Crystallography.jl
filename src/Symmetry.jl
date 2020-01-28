module Symmetry

using CoordinateTransformations: AffineMap
using LibSymspg: get_symmetry

export Cell, getsymmetry

struct Cell{L<:AbstractVecOrMat,P<:AbstractVecOrMat,N<:AbstractVector,M<:Union{AbstractVector,Nothing}}
    lattice::L
    positions::P
    numbers::N
    magmoms::M
end
Cell(lattice, positions, numbers) = Cell(lattice, positions, numbers, nothing)

function getsymmetry(cell::Cell, symprec::AbstractFloat = 1e-5)
    maps, translations = get_symmetry(cell.lattice, cell.positions, cell.numbers, length(cell.numbers), symprec)
    return (AffineMap(m, t) for (m, t) in zip(maps, translations))
end

end # module Symmetry
