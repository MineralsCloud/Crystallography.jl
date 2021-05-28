using LinearAlgebra
using Test

using SymPy: symbols

using Crystallography

@testset "Test `bravaislattices`" begin end

@testset "Test constructors" begin
    @test_throws TypeError MillerIndices{Int,Int}
    @test_throws TypeError MillerBravaisIndices{Int,Int}
    @test_throws TypeError MillerIndices{Int}(1, 2, 3)
    @test_throws TypeError MillerBravaisIndices{Int}(1, 2, 3, 4)
    @test_throws TypeError MillerIndices{Int}([1, 2, 3])
    @test_throws TypeError MillerBravaisIndices{Int}([1, 2, 3, 4])
    @test_throws TypeError MillerIndices{Int}((1, 2, 3))
    @test_throws TypeError MillerBravaisIndices{Int}((1, 2, 3, 4))
end

@testset "Test conversion between real `MillerIndices` and `MillerBravaisIndices`" begin
    millerreal = [[1, 1, 1], [-1, 0, 1], [0, -1, 1], [-1, -1, 1], [1, 0, 1], [0, 1, 1]]
    millerbravaisreal = [
        [1, 1, -2, 3],
        [-2, 1, 1, 3],
        [1, -2, 1, 3],
        [-1, -1, 2, 3],
        [2, -1, -1, 3],
        [-1, 2, -1, 3],
    ]
    for (x, y) in zip(millerreal, millerbravaisreal)
        m = MillerIndices{RealSpace}(x)
        mb = MillerBravaisIndices{RealSpace}(y)
        @test convert(MillerBravaisIndices{RealSpace}, m) == mb
        @test convert(MillerIndices{RealSpace}, mb) == m
        @test convert(MillerIndices{RealSpace}, m) == m
        @test convert(MillerBravaisIndices{RealSpace}, mb) == mb
        @test_throws MethodError convert(MillerIndices{ReciprocalSpace}, m)
        @test_throws MethodError convert(MillerBravaisIndices{ReciprocalSpace}, mb)
        @test_throws MethodError convert(MillerBravaisIndices{ReciprocalSpace}, m)
        @test_throws MethodError convert(MillerIndices{ReciprocalSpace}, mb)
        @test_throws MethodError convert(MillerBravaisIndices, m)
        @test_throws MethodError convert(MillerIndices, mb)
    end
end

@testset "Test conversion between reciprocal `MillerIndices` and `MillerBravaisIndices`" begin
    millerreciprocal =
        [[1, 0, 0], [0, 1, 0], [1, -1, 0], [-1, 0, 0], [0, -1, 0], [-1, 1, 0]]
    millerbravaisreciprocal = [
        [1, 0, -1, 0],
        [0, 1, -1, 0],
        [1, -1, 0, 0],
        [-1, 0, 1, 0],
        [0, -1, 1, 0],
        [-1, 1, 0, 0],
    ]
    for (x, y) in zip(millerreciprocal, millerbravaisreciprocal)
        m = MillerIndices{ReciprocalSpace}(x)
        mb = MillerBravaisIndices{ReciprocalSpace}(y)
        @test convert(MillerBravaisIndices{ReciprocalSpace}, m) == mb
        @test convert(MillerIndices{ReciprocalSpace}, mb) == m
        @test convert(MillerIndices{ReciprocalSpace}, m) == m
        @test convert(MillerBravaisIndices{ReciprocalSpace}, mb) == mb
        @test_throws MethodError convert(MillerIndices{RealSpace}, m)
        @test_throws MethodError convert(MillerBravaisIndices{RealSpace}, mb)
        @test_throws MethodError convert(MillerBravaisIndices{RealSpace}, m)
        @test_throws MethodError convert(MillerIndices{RealSpace}, mb)
        @test_throws MethodError convert(MillerBravaisIndices, m)
        @test_throws MethodError convert(MillerIndices, mb)
    end
end

@testset "Test length in a hexagonal lattice" begin
    g = MetricTensor(CellParameters(Hexagonal() * Primitive(), 1, 1, 2, 0, 0, 0))
    a = CrystalCoordinates(1, 2, 1)
    @test dot(a, g, a) == 7
    @test norm(CrystalCoordinates(1, 2, 1), g)^2 == 7
end

@testset "Test distance between atoms in a hexagonal lattice" begin
    g = MetricTensor(CellParameters(Hexagonal() * Primitive(), 1, 1, 2, 0, 0, 0))
    a = CrystalCoordinates(1, 1, 1)
    b = CrystalCoordinates(1 / 3, 1 / 3, 1 / 2)
    @test distance(a, g, b)^2 == 13 / 9
end

@testset "Test direction cosine in a tetragonal lattice" begin
    g = MetricTensor(CellParameters(Tetragonal() * Primitive(), 2, 2, 3, 0, 0, 0))
    a = CrystalCoordinates(1, 2, 1)
    b = CrystalCoordinates(0, 0, 1)
    @test directioncosine(a, g, b)^2 == 9 // 29
end

@testset "Symbolic calculation" begin
    a, b, c = symbols("a, b, c", positive = true)
    @test MetricTensor(CellParameters(Hexagonal() * Primitive(), a, b, c, 0, 0, 0)) ==
          MetricTensor([a^2 -a^2/2 0; -a^2/2 a^2 0; 0 0 c^2])
end # testset

@testset "Fractional coordinates to Cartesian coordinates" begin
    lattice = Lattice([
        1/2 0 0
        0 1/2 0
        0 0 1
    ])
    @test CartesianFromFractional(lattice)([2, 3, 1]) == [1, 3 / 2, 1]
end
