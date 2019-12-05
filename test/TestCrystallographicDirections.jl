#=
TestCrystallographicDirections.jl:
- Julia version: 1.0
- Author: qz
- Date: Jun 2, 2019
=#
module TestCrystallographicDirections

using Test

using Crystallography
using Crystallography.CrystallographicDirections

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

end
