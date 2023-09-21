# Example from https://github.com/spglib/spglib/blob/v2.1.0-rc2/example/python_api/example_full.py#L61-L77
@testset "Test distorted silicon structure" begin
    lattice = [(4.01, 0, 0), (0, 4, 0), (0, 0, 3.99)]
    positions = [
        (0.001, 0, 0),
        (0, 0.5, 0.5),
        (0.5, 0, 0.5),
        (0.5, 0.5, 0),
        (0.25, 0.25, 0.251),
        (0.25, 0.75, 0.75),
        (0.75, 0.25, 0.75),
        (0.75, 0.75, 0.25),
    ]
    atoms = [14, 14, 14, 14, 14, 14, 14, 14]
    cell = Cell(lattice, positions, atoms)
    mesh = [6, 5, 6]
    result = reciprocal_mesh(
        cell, mesh, 1e-5; is_shift=[true, false, true], is_time_reversal=true
    )
    allk = eachpoint(result, true)
    symops = getsymmetry(cell, 1e-5)
    for k in allk
        stabilizer, orbit = getstabilizer(symops, k), getorbit(symops, k)
        # Lagrange's theorem: | G_s | | O_s | = | G |
        @test length(stabilizer) * length(orbit) == length(symops)
    end
end

# From https://github.com/spglib/spglib/blob/d8c39f6/example/python_api/example_full.py#L259-L266
@testset "Test the primitive cell of silicon" begin
    lattice = Lattice([[0, 2, 2], [2, 0, 2], [2, 2, 0]])
    positions = [[0, 0, 0], [0.25, 0.25, 0.25]]
    atoms = [14, 14]
    cell = Cell(lattice, positions, atoms)
    mesh = [11, 11, 11]
    result = reciprocal_mesh(cell, mesh, 1e-5; is_time_reversal=true)
    allk = eachpoint(result, true)
    symops = getsymmetry(cell, 1e-5)
    for k in allk
        stabilizer, orbit = getstabilizer(symops, k), getorbit(symops, k)
        # Lagrange's theorem: | G_s | | O_s | = | G |
        @test length(stabilizer) * length(orbit) == length(symops)
    end
end

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
