function Base.show(io::IO, x::AbstractLattice)
    if get(io, :compact, false) || get(io, :typeinfo, nothing) == typeof(x)
        Base.show_default(IOContext(io, :limit => true), x)  # From https://github.com/mauro3/Parameters.jl/blob/ecbf8df/src/Parameters.jl#L556
    else
        println(io, string(typeof(x)))
        for row in eachrow(x.data)
            println(io, ' ', join(row, "  "))
        end
        println(io, " (a, b, c, α, β, γ) = ", latticeconstants(x))
    end
end
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
function Base.show(io::IO, x::ReciprocalPoint)
    if get(io, :compact, false) || get(io, :typeinfo, nothing) == typeof(x)
        Base.show_default(IOContext(io, :limit => true), x)  # From https://github.com/mauro3/Parameters.jl/blob/ecbf8df/src/Parameters.jl#L556
    else
        println(io, string(typeof(x)))
        print(io, " coord = ", x.coord, ", weight = ", x.weight)
    end
end
function Base.show(io::IO, x::ChangeOfBasis)
    if get(io, :compact, false) || get(io, :typeinfo, nothing) == typeof(x)
        Base.show_default(IOContext(io, :limit => true), x)  # From https://github.com/mauro3/Parameters.jl/blob/ecbf8df/src/Parameters.jl#L556
    else
        println(io, string(typeof(x)))
        for row in eachrow(x.tf)
            println(io, ' ', join(row, "  "))
        end
    end
end
