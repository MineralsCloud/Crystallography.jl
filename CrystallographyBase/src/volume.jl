using PeriodicTable: Element, elements

export cellvolume, density

"""
    cellvolume(a, b, c, α, β, γ)

Calculates the cell volume from 6 cell parameters.
"""
cellvolume(a, b, c, α, β, γ) =
    a * b * c * sqrt(sind(α)^2 - cosd(β)^2 - cosd(γ)^2 + 2 * cosd(α) * cosd(β) * cosd(γ))
"""
    cellvolume(l::Lattice)
    cellvolume(c::Cell)

Calculates the cell volume from a `Lattice` or a `Cell`.
"""
cellvolume(lattice::AbstractLattice) = abs(det(lattice.data))
cellvolume(cell::Cell) = cellvolume(cell.lattice)
"""
    cellvolume(g::MetricTensor)

Calculates the cell volume from a `MetricTensor`.
"""
cellvolume(g::MetricTensor) = sqrt(det(g.data))  # `sqrt` is always positive!

density(volume::Number, mass::Number) = mass / volume
function density(lattice::Lattice, atoms)
    mass = sum(atomicmass, atoms)
    volume = cellvolume(lattice)
    return mass / volume
end
density(cell::Cell) = density(Lattice(cell.lattice), cell.types)

atomicmass(element::Element) = element.atomic_mass
atomicmass(i::Union{AbstractString,Integer,Symbol}) = elements[i].atomic_mass
