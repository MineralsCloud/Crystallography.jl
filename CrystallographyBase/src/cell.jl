using StaticArrays: MVector
using StructHelpers: @batteries

export Cell, natoms, eachatom

struct Cell{L,P,T}
    lattice::Lattice{L}
    positions::Vector{MVector{3,P}}
    atoms::Vector{T}
end
function Cell(lattice, positions, atoms)
    if !(lattice isa Lattice)
        lattice = Lattice(lattice)
    end
    if positions isa AbstractVector
        P = eltype(Base.promote_typeof(positions...))
        positions = collect(map(MVector{3,P}, positions))
    else
        throw(ArgumentError("`positions` must be a `Vector` of `Vector`s!"))
    end
    L, T = eltype(lattice), eltype(atoms)
    return Cell{L,P,T}(lattice, positions, atoms)
end

@batteries Cell eq = true hash = true

natoms(cell::Cell) = length(cell.atoms)

"""
    Lattice(cell::Cell)

Get the lattice of a `Cell`.
"""
Lattice(cell::Cell) = cell.lattice

struct EachAtom{A,B}
    atoms::Vector{A}
    positions::Vector{B}
end
EachAtom(cell::Cell) = EachAtom(cell.atoms, cell.positions)

eachatom(cell::Cell) = EachAtom(cell)

# Similar to https://github.com/JuliaCollections/IterTools.jl/blob/0ecaa88/src/IterTools.jl#L1028-L1032
function Base.iterate(itr::EachAtom, state = 1)
    if state > length(itr)
        return nothing
    else
        return (itr.atoms[state], itr.positions[state]), state + 1
    end
end

Base.eltype(::EachAtom{A,B}) where {A,B} = Tuple{A,B}

Base.length(itr::EachAtom) = length(itr.atoms)

Base.IteratorSize(::Type{<:EachAtom}) = Base.HasLength()

function Base.show(io::IO, cell::Cell)
    if get(io, :compact, false) || get(io, :typeinfo, nothing) == typeof(cell)
        Base.show_default(IOContext(io, :limit => true), cell)  # From https://github.com/mauro3/Parameters.jl/blob/ecbf8df/src/Parameters.jl#L556
    else
        println(io, string(typeof(cell)))
        println(io, " lattice:")
        for row in eachrow(cell.lattice.data)
            println(io, "  ", join(row, "  "))
        end
        N = natoms(cell)
        println(io, " $N atomic positions:")
        for pos in cell.positions
            println(io, "  ", pos)
        end
        println(io, " $N atoms:")
        println(io, "  ", cell.atoms)
    end
end
