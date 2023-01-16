export cellvolume

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
cellvolume(lattice::AbstractLattice) = abs(det(lattice.data))
cellvolume(cell::Cell) = cellvolume(Lattice(cell))
"""
    cellvolume(g::MetricTensor)

Calculate the cell volume from a `MetricTensor`.
"""
cellvolume(g::MetricTensor) = sqrt(det(g.data))  # `sqrt` is always positive!
