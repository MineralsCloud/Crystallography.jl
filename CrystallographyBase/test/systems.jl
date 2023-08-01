@testset "Test construction from symbols and numbers" begin
    @test LatticeSystem(1) == LatticeSystem(:Triclinic)
    @test LatticeSystem(:Monoclinic) == LatticeSystem(2)
    @test LatticeSystem(:Rhombohedral) == LatticeSystem(5)
    @test CrystalSystem(:Triclinic) == CrystalSystem(1)
    @test CrystalSystem(4) == CrystalSystem(:Tetragonal)
    @test CrystalSystem(:Trigonal) == CrystalSystem(5)
    @test Bravais(:PrimitiveTetragonal) == Bravais(8)
    @test Bravais(:BodyCenteredTetragonal) == Bravais(9)
    @test Bravais(14) == Bravais(:FaceCenteredCubic)
end

@testset "Test instances of different types are different" begin
    @test LatticeSystem(4) != CrystalSystem(3)
    @test LatticeSystem(:Tetragonal) != Bravais(:PrimitiveTetragonal)
    @test CrystalSystem(:Monoclinic) != Bravais(:BaseCenteredMonoclinic)
end

@testset "Test instances of different types are different even if they have the same name" begin
    @test LatticeSystem(:Orthorhombic) != CrystalSystem(:Orthorhombic)
    @test LatticeSystem(:Cubic) != CrystalSystem(:Cubic)
    @test LatticeSystem(4) != CrystalSystem(4)
end

@testset "Test `instances`" begin
    @test length(instances(LatticeSystem)) == 7
    @test length(instances(CrystalSystem)) == 7
    @test length(instances(Bravais)) == 14
end
