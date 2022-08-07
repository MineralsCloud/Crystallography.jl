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
    @test latticesystem(lattice) == LatticeSystem.Orthorhombic
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
