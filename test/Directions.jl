module Directions

using LinearAlgebra
using Test

using Crystallography
using Crystallography.Directions

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
        @test_throws ErrorException convert(MillerIndices{ReciprocalSpace}, m)
        @test_throws ErrorException convert(MillerBravaisIndices{ReciprocalSpace}, mb)
        @test_throws ErrorException convert(MillerBravaisIndices{ReciprocalSpace}, m)
        @test_throws ErrorException convert(MillerIndices{ReciprocalSpace}, mb)
        @test_throws ErrorException convert(MillerBravaisIndices, m)
        @test_throws ErrorException convert(MillerIndices, mb)
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
        @test_throws ErrorException convert(MillerIndices{RealSpace}, m)
        @test_throws ErrorException convert(MillerBravaisIndices{RealSpace}, mb)
        @test_throws ErrorException convert(MillerBravaisIndices{RealSpace}, m)
        @test_throws ErrorException convert(MillerIndices{RealSpace}, mb)
        @test_throws ErrorException convert(MillerBravaisIndices, m)
        @test_throws ErrorException convert(MillerIndices, mb)
    end
end

@testset "Test length in a hexagonal lattice" begin
    g = MetricTensor(BravaisLattice(Hexagonal(), Primitive()), 1, 2)
    a = CrystalCoordinates(1, 2, 1)
    @test dot(a, g, a) ≈ 7
    @test norm(CrystalCoordinates(1, 2, 1), g) == sqrt(7)
end

@testset "Test distance between atoms in a hexagonal lattice" begin
    g = MetricTensor(BravaisLattice(Hexagonal(), Primitive()), 1, 2)
    a = CrystalCoordinates(1, 1, 1)
    b = CrystalCoordinates(1 / 3, 1 / 3, 1 / 2)
    @test distance(a, g, b) ≈ sqrt(13) / 3
end

@testset "Test direction cosine in a tetragonal lattice" begin
    g = MetricTensor(BravaisLattice(Tetragonal(), Primitive()), 2, 3)
    a = CrystalCoordinates(1, 2, 1)
    b = CrystalCoordinates(0, 0, 1)
    @test directioncosine(a, g, b) ≈ 3 / sqrt(29)
end

end
