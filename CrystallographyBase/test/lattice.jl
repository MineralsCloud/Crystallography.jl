@testset "Test creating `Lattice` with units" begin
    a = 4u"nm"
    b = 180u"bohr"
    c = 3u"angstrom"
    @test Lattice(a, b, c, 90, 90, 90) == Lattice(
        [
            4u"nm" 0u"m" 0.0u"cm"
            0u"cm" 180.0u"bohr" 0u"m"
            0u"bohr" 0u"nm" (3//1)*u"angstrom"
        ],
    )
end
