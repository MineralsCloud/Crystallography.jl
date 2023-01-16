export crystaldensity, atomicmass

"""
    crystaldensity(lattice::Lattice, atoms)
    crystaldensity(cell::Cell)

Calculate the density of a crystal structure.

Here, `atoms` is an iterable of atomic numbers, element names, symbols, or `PeriodicTable.Element`s.
You can extend the `atomicmass` method to work with custom types.
"""
function crystaldensity(lattice::Lattice, atoms)
    mass = sum(atomicmass, atoms)
    volume = cellvolume(lattice)
    return mass / volume
end
crystaldensity(cell::Cell) = crystaldensity(Lattice(cell), cell.atoms)
const density = crystaldensity  # For compatibility reason

atomicmass(element::Element) = element.atomic_mass
atomicmass(i::Union{AbstractString,Integer,Symbol}) = elements[i].atomic_mass
