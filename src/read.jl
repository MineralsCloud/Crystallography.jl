using CrystalInfoFramework: Cif, get_loop

export readcif

readcif(filename, ::Type{T}) where {T} = open(io -> readcif(io, T), filename)
function readcif(io::IO, ::Type{T}) where {T<:Lattice}
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
        parse(_eltype(T), value)
    end
    return Lattice(a, b, c, α, β, γ)
end
function readcif(io::IO, ::Type{T}) where {T<:Cell}
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
        parse(_lattice_eltype(T), value)
    end
    lattice = Lattice(a, b, c, α, β, γ)
    positions = map(eachrow(get_loop(block, "_atom_site_type_symbol"))) do row
        map(
            Base.Fix1(parse, _position_eltype(T)),
            (
                row["_atom_site_fract_x"],
                row["_atom_site_fract_y"],
                row["_atom_site_fract_z"],
            ),
        )
    end
    atoms = map(eachrow(get_loop(block, "_atom_site_type_symbol"))) do row
        _atom_type(T)(row["_atom_site_type_symbol"])
    end
    return Cell(lattice, positions, atoms)
end

_eltype(::Type{Lattice{T}}) where {T} = T
_eltype(::Type{<:Lattice}) = Float64

_lattice_eltype(::Type{Cell{T}}) where {T} = T
_lattice_eltype(::Type{<:Cell}) = Float64

_position_eltype(::Type{Cell{T,S}}) where {T,S} = S
_position_eltype(::Type{<:Cell}) = Float64

_atom_type(::Type{Cell{T,S,R}}) where {T,S,R} = R
_atom_type(::Type{<:Cell}) = String
