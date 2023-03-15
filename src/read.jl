using CrystalInfoFramework: Cif, get_loop

export readcif

readcif(filename) = open(io -> readcif(io), filename)
function readcif(io::IO)
    str = read(io, String)
    cif = Cif(str)
    block = first(cif).second
    a, b, c, α, β, γ = map((
        "_cell_length_a",
        "_cell_length_b",
        "_cell_length_c",
        "_cell_angle_alpha",
        "_cell_angle_beta",
        "_cell_angle_gamma",
    )) do key
        value = only(block[key])
        parse(Float64, value)
    end
    lattice = Lattice(a, b, c, α, β, γ)
    positions = map(eachrow(get_loop(block, "_atom_site_type_symbol"))) do row
        map(
            Base.Fix1(parse, Float64),
            (
                row["_atom_site_fract_x"],
                row["_atom_site_fract_y"],
                row["_atom_site_fract_z"],
            ),
        )
    end
    atoms = map(eachrow(get_loop(block, "_atom_site_type_symbol"))) do row
        row["_atom_site_type_symbol"]
    end
    return Cell(lattice, positions, atoms)
end
