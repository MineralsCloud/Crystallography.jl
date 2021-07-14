using Brillouin

@testset "Test `WignerSeitzCell` against `Brillouin`" begin
    # Data from https://thchr.github.io/Brillouin.jl/dev/wignerseitz/
    lattice = Lattice([0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0])
    ws = WignerSeitzCell(lattice)
    cᴿ = wignerseitz([[0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]])
    @test ws.verts == [
        [0.25, -0.75, 0.25],
        [-0.5, -0.5, 0.5],
        [-0.75, 0.25, 0.25],
        [0.75, -0.25, -0.25],
        [-0.25, -0.25, -0.25],
        [0.25, 0.25, -0.75],
        [-0.5, 0.5, -0.5],
        [0.5, -0.5, -0.5],
        [-0.25, -0.25, 0.75],
        [0.5, -0.5, 0.5],
        [-0.25, 0.75, -0.25],
        [0.25, 0.25, 0.25],
        [0.5, 0.5, -0.5],
        [-0.5, 0.5, 0.5],
    ]
    @test ws.faces == [
        [3, 7, 5, 2],
        [6, 8, 5, 7],
        [7, 3, 14, 11],
        [11, 13, 6, 7],
        [14, 3, 2, 9],
        [1, 2, 5, 8],
        [1, 10, 9, 2],
        [13, 11, 14, 12],
        [14, 9, 10, 12],
        [4, 8, 6, 13],
        [4, 10, 1, 8],
        [13, 12, 10, 4],
    ]
    @test cartesianize(cᴿ).verts == WignerSeitzCell(lattice, true).verts
end
