module IO

using CrystallographyBase: Lattice, Cell, natoms, eachatom, latticeconstants

export writexyz

function writexyz(io::Base.IO, cell::Cell, comment = "")
    write(io, natoms(cell))
    write(io, comment)
    for (atom, positions) in eachatom(cell)
        write(io, string(atom)[1:2], join(positions, "  "))
    end
    return
end
function writexyz(filename, cell::Cell, comment = "")
    open(filename, "w") do io
        writexyz(io, cell, comment)
    end
end

function writecif(io::Base.IO, cell::Cell, comment = "")
    a, b, c, α, β, γ = latticeconstants(Lattice(cell))
    block = CifBlock()
    for (key, value) in zip(
        (
            "_cell_length_a",
            "_cell_length_b",
            "_cell_length_c",
            "_cell_angle_alpha",
            "_cell_angle_beta",
            "_cell_angle_gamma",
        ),
        latticeconstants(Lattice(cell)),
    )
        block[key] = value
    end
    for (key, value) in symmetry_items.items()
        block[key] = value
    end
end
function writecif(filename, cell::Cell, comment = "")
    open(filename, "w") do io
        writexyz(io, cell, comment)
    end
end

end
