# Example from https://github.com/spglib/spglib/blob/v2.1.0-rc2/example/python_api/example_full.py#L111-L117
@testset "Test MgB₂ structure" begin
    a = 3.07
    c = 3.52
    lattice = Lattice([[a, 0, 0], [-a / 2, a / 2 * sqrt(3), 0], [0, 0, c]])
    positions = [[0, 0, 0], [1 / 3, 2 / 3, 1 / 2], [2 / 3, 1 / 3, 1 / 2]]
    atoms = [12, 5, 5]
    cell = Cell(lattice, positions, atoms)
    mesh = [7, 7, 7]
    result = reciprocal_mesh(cell, mesh, 1e-5; is_time_reversal=true)
    allk = eachpoint(result, true)
    symops = getsymmetry(cell, 1e-5)
    for k in allk
        stabilizer, orbit = getstabilizer(symops, k), getorbit(symops, k)
        # Lagrange's theorem: | G_s | | O_s | = | G |
        @test length(stabilizer) * length(orbit) == length(symops)
    end
end
