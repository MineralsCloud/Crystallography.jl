using Brillouin: wignerseitz

import Brillouin: vertices, faces

export WignerSeitzCell, vertices, faces

struct WignerSeitzCell{V,F,T<:AbstractLattice}
    vertices::SVector{V,SVector{3,Float64}}
    faces::SVector{F,Vector{Int}}
    lattice::T
end
"""
    WignerSeitzCell(lattice::AbstractLattice)

Construct a [Wigner–Seitz cell](https://en.wikipedia.org/wiki/Wigner%E2%80%93Seitz_cell) from an `AbstractLattice`.
"""
function WignerSeitzCell(lattice::T) where {T<:AbstractLattice}
    ws = wignerseitz(collect(primitivevectors(lattice)))
    return WignerSeitzCell{length(ws.verts),length(ws.faces),T}(ws.verts, ws.faces, lattice)
end

"""
    vertices(ws::WignerSeitzCell, cartesian=false)

Get the coordinates of the vertices of a Wigner–Seitz cell.

If `cartesian` is `true`, return the coordinates in the Cartesian coordinate system.
"""
vertices(ws::WignerSeitzCell, cartesian = false) =
    cartesian ? map(CartesianFromFractional(ws.lattice), ws.vertices) : ws.vertices

"""
    faces(ws::WignerSeitzCell, cartesian=false)

Get the coordinates of the vertices of each face of a Wigner–Seitz cell.

If `cartesian` is `true`, return the coordinates in the Cartesian coordinate system.
"""
faces(ws::WignerSeitzCell, cartesian = false) =
    map(I -> vertices(ws, cartesian)[I], ws.faces)

# Referenced from https://github.com/thchr/Brillouin.jl/blob/f32a826/src/WignerSeitz.jl#L59-L78
function Base.show(io::IO, ::MIME"text/plain", x::WignerSeitzCell)
    if get(io, :compact, false) || get(io, :typeinfo, nothing) == typeof(x)
        Base.show_default(IOContext(io, :limit => true), x)  # From https://github.com/mauro3/Parameters.jl/blob/ecbf8df/src/Parameters.jl#L556
    else
        summary(io, x)
        println(io, ":")
        println(io, "  vertices: ")
        foreach(v -> println(io, "    ", v), x.vertices)
        println(io, "  faces: ")
        foreach(f -> println(io, "    ", f), x.faces)
        println(io, "  base lattice: ")
        for row in eachrow(x.lattice.data)
            println(io, "    ", join(row, "  "))
        end
    end
end
