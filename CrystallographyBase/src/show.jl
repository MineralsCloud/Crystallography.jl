function Base.show(io::IO, lattice::AbstractLattice)
    if get(io, :compact, false)
        print(io, typeof(lattice), latticeconstants(lattice))
    else
        print(io, typeof(lattice), '(')
        join(
            io,
            (
                "$x=$y" for
                (x, y) in zip(('a', 'b', 'c', 'α', 'β', 'γ'), latticeconstants(lattice))
            ),
            ", ",
        )
        print(io, ')')
    end
end
function Base.show(io::IO, ::MIME"text/plain", lattice::AbstractLattice)
    summary(io, lattice)
    println(io)
    join(io, ' ' * join(row, "  ") * '\n' for row in eachrow(lattice.data))
    return nothing
end
function Base.show(io::IO, ::MIME"text/plain", cell::Cell)
    summary(io, cell)
    println(io)
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
    return nothing
end
function Base.show(io::IO, ::MIME"text/plain", point::ReciprocalPoint)
    println(io, string(typeof(point)))
    println(io, " coord=", point.coord)
    print(io, " weight=", point.weight)
    return nothing
end
function Base.show(io::IO, ::MIME"text/plain", A::ChangeOfBasis)
    println(io, string(typeof(A)))
    for row in eachrow(A.tf)
        println(io, ' ', join(row, "  "))
    end
end
