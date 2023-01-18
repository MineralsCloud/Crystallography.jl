using Mendeleev: Element, elements

export cellvolume, crystaldensity, atomicmass

"""
    cellvolume(a, b, c, α, β, γ)

Calculate the cell volume from 6 cell parameters.
"""
cellvolume(a, b, c, α, β, γ) =
    a * b * c * sqrt(sind(α)^2 - cosd(β)^2 - cosd(γ)^2 + 2 * cosd(α) * cosd(β) * cosd(γ))
"""
    cellvolume(l::Lattice)
    cellvolume(c::Cell)

Calculate the cell volume from a `Lattice` or a `Cell`.
"""
cellvolume(lattice::AbstractLattice) = abs(_det(lattice.data))
cellvolume(cell::Cell) = cellvolume(Lattice(cell))
"""
    cellvolume(g::MetricTensor)

Calculate the cell volume from a `MetricTensor`.
"""
cellvolume(g::MetricTensor) = sqrt(_det(g.data))  # `sqrt` is always positive!

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
