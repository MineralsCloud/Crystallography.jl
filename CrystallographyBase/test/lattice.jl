@testset "Test creating `Lattice`s with units" begin
    a = 4u"nm"
    b = 180u"bohr"
    c = 3u"angstrom"
    lattice = Lattice(a, b, c, 90, 90, 90)
    @test lattice == Lattice(
        [
            4u"nm" 0u"m" 0.0u"cm"
            0u"cm" 180.0u"bohr" 0u"m"
            0u"bohr" 0u"nm" (3//1)*u"angstrom"
        ]
    )
    @test latticesystem(lattice; lengthtol=1e-5u"bohr") == LatticeSystem.Orthorhombic
end

@testset "Test creating `Lattice`s from 6 lattice constants" begin
    # From https://github.com/LaurentRDC/crystals/blob/7c544fe/crystals/tests/test_lattice.py#L96
    @test basisvectors(Lattice(2, 1, 5, 90, 90, 90; axis=:c)) ==
        ([2, 0, 0], [0, 1, 0], [0, 0, 5])  # Orthorombic
    # From https://github.com/LaurentRDC/crystals/blob/7c544fe/crystals/tests/test_lattice.py#L104
    @test basisvectors(Lattice(1, 2, 3, 90, 120, 90; axis=:c)) ==
        ([0.8660254037844387, 0, -0.5], [0, 2, 0], [0, 0, 3])  # Monoclinic
    # From https://github.com/LaurentRDC/crystals/blob/7c544fe/crystals/tests/test_lattice.py#L117
    # Compared with Python's results
    @test all(
        basisvectors(Lattice(1, 2, 3, 75, 40, 81; axis=:c)) .≈
        ([0.641327, -0.04330811, 0.76604444], [0, 1.93185165, 0.51763809], [0, 0, 3]),
    )  # Triclinic
    # From https://github.com/LaurentRDC/crystals/blob/7c544fe/crystals/tests/test_lattice.py#L131
    # Compared with Python's results
    @test all(
        basisvectors(Lattice(3, 4, 20, 45, 90, 126; axis=:c)) .≈ (
            [1.66767891, -2.49376163, 1.83697020e-16],
            [0, 2.82842712, 2.82842712],
            [0, 0, 20],
        ),
    )  # Triclinic
    # From https://github.com/LaurentRDC/crystals/blob/7c544fe/crystals/tests/test_lattice.py#L158
    @test all(
        basisvectors(Lattice(2, 2, 3, 90, 90, 90; axis=:c)) .≈
        ([2, 0, 0], [0, 2, 0], [0, 0, 3]),
    )  # Tetragonal
    # From https://github.com/LaurentRDC/crystals/blob/7c544fe/crystals/tests/test_lattice.py#L165
    # Compared with Python's results
    @test all(
        basisvectors(Lattice(1, 1, 1, 87, 87, 87; axis=:c)) .≈
        ([0.99739377, 0.04966497, 0.05233596], [0, 0.99862953, 0.05233596], [0, 0, 1]),
    )  # Rhombohedral
    # From https://github.com/LaurentRDC/crystals/blob/7c544fe/crystals/tests/test_lattice.py#L173
    # Compared with Python's results
    @test all(
        basisvectors(Lattice(1, 2, 3, 90, 115, 90; axis=:c)) .≈
        ([0.906307787, 8.71102450e-17, -0.422618262], [0, 2, 1.2246468e-16], [0, 0, 3]),
    )  # Monoclinic
    # From https://github.com/LaurentRDC/crystals/blob/7c544fe/crystals/tests/test_lattice.py#L177
    # Compared with Python's results
    @test all(
        basisvectors(Lattice(2, 3, 1, 115, 90, 90; axis=:c)) .≈
        ([2, 1.92231042e-16, 1.22464680e-16], [0, 2.71892336, -1.26785479], [0, 0, 1]),
    )  # Monoclinic
    # From https://github.com/LaurentRDC/crystals/blob/7c544fe/crystals/tests/test_lattice.py#L181
    # Compared with Python's results
    @test all(
        basisvectors(Lattice(3, 1, 2, 90, 90, 115; axis=:c)) .≈
        ([2.71892336, -1.26785479, 1.83697020e-16], [0, 1, 6.123234e-17], [0, 0, 2]),
    )  # Monoclinic
    # From https://github.com/LaurentRDC/crystals/blob/7c544fe/crystals/tests/test_lattice.py#L189
    # Compared with Python's results
    @test all(
        basisvectors(Lattice(2, 2, 3, 90, 90, 120; axis=:c)) .≈
        ([1.73205081, -1, 1.22464680e-16], [0, 2, 1.2246468e-16], [0, 0, 3]),
    )  # Hexagonal
    # From https://github.com/LaurentRDC/crystals/blob/7c544fe/crystals/tests/test_lattice.py#L193
    # Compared with Python's results
    @test all(
        basisvectors(Lattice(3, 2, 2, 120, 90, 90; axis=:c)) .≈
        ([3, 3.18172572e-16, 1.83697020e-16], [0, 1.73205081, -1], [0, 0, 2]),
    )  # Hexagonal
    # From https://github.com/LaurentRDC/crystals/blob/7c544fe/crystals/tests/test_lattice.py#L197
    # Compared with Python's results
    @test all(
        basisvectors(Lattice(2, 3, 2, 90, 120, 90; axis=:c)) .≈
        ([1.73205081, 1.83697020e-16, -1], [0, 3, 1.8369702e-16], [0, 0, 2]),
    )  # Hexagonal
    # From https://github.com/LaurentRDC/crystals/blob/7c544fe/crystals/tests/test_lattice.py#L201
    # Compared with Python's results
    @test all(
        basisvectors(Lattice(2, 2, 2, 90, 120, 90; axis=:c)) .≈
        ([1.73205081, 1.83697020e-16, -1], [0, 2, 1.2246468e-16], [0, 0, 2]),
    )  # Hexagonal
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

@testset "Test `periodicity`" begin
    # From https://github.com/LaurentRDC/crystals/blob/7c544fe/crystals/tests/test_lattice.py#L117
    # Compared with Python's results
    @test collect(
        periodicity(Lattice(1u"angstrom", 2u"angstrom", 3u"angstrom", 75, 40, 81; axis=:c))
    ) ≈ [0.6413269980200662, 1.975159767026871, 4.283682533324019] * u"angstrom"
    # From https://github.com/LaurentRDC/crystals/blob/7c544fe/crystals/tests/test_lattice.py#L131
    # Compared with Python's results
    @test collect(periodicity(Lattice(3, 4, 20, 45, 90, 126; axis=:c))) ≈
        [1.6676789107542636, 5.322188751410911, 22.82842712474619]
    # From https://github.com/LaurentRDC/crystals/blob/7c544fe/crystals/tests/test_lattice.py#L181
    # Compared with Python's results
    @test collect(periodicity(Lattice(3, 1, 2, 90, 90, 115; axis=:c))) ≈
        [2.71892336110995, 2.267854785222098, 2.0000000000000004]
end

@testset "Test `latticeconstants`" begin
    # From https://github.com/LaurentRDC/crystals/blob/7c544fe/crystals/tests/test_lattice.py#L96
    @test collect(latticeconstants(Lattice(2, 1, 5, 90, 90, 90))) ≈ [2, 1, 5, 90, 90, 90]  # Orthorombic
    # From https://github.com/LaurentRDC/crystals/blob/7c544fe/crystals/tests/test_lattice.py#L104
    @test collect(latticeconstants(Lattice(1, 2, 3, 90, 120, 90))) ≈ [1, 2, 3, 90, 120, 90]  # Monoclinic
    # From https://github.com/LaurentRDC/crystals/blob/7c544fe/crystals/tests/test_lattice.py#L117
    @test collect(latticeconstants(Lattice(1, 2, 3, 75, 40, 81))) ≈ [1, 2, 3, 75, 40, 81]  # Triclinic
end
