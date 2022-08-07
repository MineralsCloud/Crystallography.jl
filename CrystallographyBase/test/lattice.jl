@testset "Test creating `Lattice` with units" begin
    a = 4u"nm"
    b = 180u"bohr"
    c = 3u"angstrom"
    lattice = Lattice(a, b, c, 90, 90, 90)
    @test lattice == Lattice(
        [
            4u"nm" 0u"m" 0.0u"cm"
            0u"cm" 180.0u"bohr" 0u"m"
            0u"bohr" 0u"nm" (3//1)*u"angstrom"
        ],
    )
    @test latticesystem(lattice; lengthtol = 1e-5u"bohr") == LatticeSystem.Orthorhombic
end

# See https://github.com/LaurentRDC/crystals/blob/7c544fe/crystals/tests/test_lattice.py#L156-L209
@testset "Test `latticesystem` from 6 lattice constants" begin
    @test latticesystem(2, 2, 3, 90, 90, 90) == LatticeSystem.Tetragonal
    @test latticesystem(1, 1, 1, 87, 87, 87) == LatticeSystem.Rhombohedral
    @test latticesystem(1, 2, 3, 90, 115, 90) == LatticeSystem.Monoclinic
    @test latticesystem(2, 3, 1, 115, 90, 90) == LatticeSystem.Monoclinic
    @test latticesystem(3, 1, 2, 90, 90, 115) == LatticeSystem.Monoclinic
    @test latticesystem(2, 2, 3, 90, 90, 120) == LatticeSystem.Hexagonal
    @test latticesystem(3, 2, 2, 120, 90, 90) == LatticeSystem.Hexagonal
    @test latticesystem(2, 3, 2, 90, 120, 90) == LatticeSystem.Hexagonal
    @test latticesystem(2, 2, 2, 90, 120, 90) == LatticeSystem.Hexagonal
    @test latticesystem(1, 2, 3, 75, 40, 81) == LatticeSystem.Triclinic
end

@testset "Test `latticeconstants`" begin
    # From https://github.com/LaurentRDC/crystals/blob/7c544fe/crystals/tests/test_lattice.py#L96
    @test collect(latticeconstants(Lattice(2, 1, 5, 90, 90, 90))) ≈ [2, 1, 5, 90, 90, 90]  # Orthorombic
    # From https://github.com/LaurentRDC/crystals/blob/7c544fe/crystals/tests/test_lattice.py#L104
    @test collect(latticeconstants(Lattice(1, 2, 3, 90, 120, 90))) ≈ [1, 2, 3, 90, 120, 90]  # Monoclinic
    # From https://github.com/LaurentRDC/crystals/blob/7c544fe/crystals/tests/test_lattice.py#L117
    @test collect(latticeconstants(Lattice(1, 2, 3, 75, 40, 81))) ≈ [1, 2, 3, 75, 40, 81]  # Triclinic
end
