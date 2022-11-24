module IO

using CrystalInfoFramework: CifBlock, add_to_loop
using CrystallographyBase: Lattice, Cell, natoms, eachatom, latticeconstants, latticevectors

export writexyz

function writexyz(io::Base.IO, cell::Cell, comment="")
    write(io, natoms(cell))
    write(io, comment)
    for (atom, positions) in eachatom(cell)
        write(io, string(atom)[1:2], join(positions, "  "))
    end
    return
end
function writexyz(filename, cell::Cell, comment="")
    open(filename, "w") do io
        writexyz(io, cell, comment)
    end
end

function writecif(io::Base.IO, cell::Cell)
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
    add_to_loop(block, "datanames", [
        "_atom_site_type_symbol",
        "_atom_site_fract_x",
        "_atom_site_fract_y",
        "_atom_site_fract_z",
    ])
    add_to_loop(block, "length_check", false)
    atoms = list(crystal.primitive().unitcell)
    symbols = [atm.symbol for atm in atoms]
    xf = [atm.coords_fractional[0] for atm in atoms]
    yf = [atm.coords_fractional[1] for atm in atoms]
    zf = [atm.coords_fractional[2] for atm in atoms]
    block["_atom_site_type_symbol"] = symbols
    block["_atom_site_fract_x"] = xf
    block["_atom_site_fract_y"] = yf
    block["_atom_site_fract_z"] = zf
end
function writecif(filename, cell::Cell)
    open(filename, "w") do io
        writecif(io, cell)
    end
end

function writevasp(io::Base.IO, cell::Cell; scaling_factor=1, cartesian=false, comment="")
    println(io, comment)
    println(io, scaling_factor)
    lattice = Lattice(cell) / scaling_factor
    for vector in latticevectors(lattice)
        println(io, join(vector, "  "))
    end


    if cartesian
        println(io, "Cartesian")
    else
        println(io, "Direct")
        scaling_factor = 1
    end
    for (elem, atoms) in grouped
        for atom in atoms
            positions = positions / scaling_factor
            println(io, join(positions, " "))
        end
    end
end

end
