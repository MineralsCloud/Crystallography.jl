#=
TestMetric.jl:
- Julia version: 1.0
- Author: qz
- Date: Jun 3, 2019
=#
module TestMetric

using Test

using LinearAlgebra

using Crystallography
using Crystallography.Metric

@testset "Test length in a hexagonal lattice" begin
    g = MetricTensor{RealSpace}(Hexagonal, 1, 2)
    a = CrystalCoordinates{RealSpace}(1, 2, 1)
    @test dot(a, g, a) ≈ 7
    @test length(CrystalCoordinates{RealSpace}(1, 2, 1), g) == sqrt(7)
end

@testset "Test distance between atoms in a hexagonal lattice" begin
    g = MetricTensor{RealSpace}(Hexagonal, 1, 2)
    a = CrystalCoordinates{RealSpace}(1, 1, 1)
    b = CrystalCoordinates{RealSpace}(1 / 3, 1 / 3, 1 / 2)
    @test distance(a, g, b) ≈ sqrt(13) / 3
end

@testset "Test direction cosine in a tetragonal lattice" begin
    g = MetricTensor{RealSpace}(Tetragonal, 2, 3)
    a = CrystalCoordinates{RealSpace}(1, 2, 1)
    b = CrystalCoordinates{RealSpace}(0, 0, 1)
    @test directioncosine(a, g, b) ≈ 3 / sqrt(29)
end

end
