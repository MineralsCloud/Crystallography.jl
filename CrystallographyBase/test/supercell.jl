@testset "Test the simplest `supercell` for a `Lattice`" begin
    lattice = Lattice(2, 1, 5, 90, 90, 90; axis = :c)
    a = supercell(lattice, 3)
    b = Lattice(6, 3, 15, 90, 90, 90; axis = :c)
    @test a == b
    @test reciprocal(a) == reciprocal(b)
end

@testset "Test `supercell` from primitive to body-centered for a `Lattice`" begin
    lattice = Lattice(2, 1, 5, 90, 90, 90; axis = :c)
    P = [
        0 1 1
        1 0 1
        1 1 0
    ]
    a = supercell(lattice, P)
    b = StandardizedFromPrimitive(P)(a)
    @test a == b
    @test reciprocal(a) == reciprocal(b)
end

@testset "Test `supercell` from orthogonal to hexagonal for a `Lattice`" begin
    lattice = Lattice(2, 1, 5, 90, 90, 90; axis = :c)
    P = [
        1 1 0
        0 2 0
        0 0 1
    ]
    a = supercell(lattice, P)
    b = StandardizedFromPrimitive(P)(a)
    @test a == b
    @test reciprocal(a) == reciprocal(b)
end
