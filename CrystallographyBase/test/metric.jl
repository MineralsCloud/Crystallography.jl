using LinearAlgebra: dot, norm, diagm
# using SymPy: symbols

@testset "Test length in a hexagonal lattice" begin
    g = MetricTensor(1, 1, 2, 90, 90, 120)  # Primitive hexagonal
    @test g ≈ MetricTensor([1, 0, 0], [-1 // 2, sqrt(3) / 2, 0], [0, 0, 2])
    @test g ≈ MetricTensor(Lattice([
        1.0 -1/2 0.0
        0.0 sqrt(3)/2 0.0
        0.0 0.0 2.0
    ]))
    a = [1, 2, 1]
    @test dot(a, g, a) ≈ 7
    @test norm([1, 2, 1], g)^2 ≈ 7
end

@testset "Test distance between atoms in a hexagonal lattice" begin
    g = MetricTensor(1, 1, 2, 90, 90, 120)  # Primitive hexagonal
    a = [1, 1, 1]
    b = [1 // 3, 1 // 3, 1 // 2]
    @test distance(a, g, b)^2 == 13 / 9
end

@testset "Test direction cosine in a tetragonal lattice" begin
    g = MetricTensor(2, 2, 3, 90, 90, 90)  # Primitive tetragonal
    @test g == MetricTensor(Lattice(diagm([2, 2, 3])))
    @test g == MetricTensor([2, 0, 0], [0, 2, 0], [0, 0, 3])
    a = [1, 2, 1]
    b = [0, 0, 1]
    @test directioncosine(a, g, b)^2 ≈ 9 / 29
end

# @testset "Symbolic calculation" begin
#     a, c = symbols("a, c", positive = true)
#     @test MetricTensor(a, a, c, 90, 90, 120) ==
#           MetricTensor([a^2 -0.5*a^2 0; -0.5*a^2 a^2 0; 0 0 c^2])  # Primitive hexagonal
# end

# From https://ssd.phys.strath.ac.uk/wp-content/uploads/Crystallographic_maths.pdf
@testset "Find angle θ between [100] and [111] directions in GaN" begin
    a, c = 3.19u"angstrom", 5.19u"angstrom"
    g = MetricTensor(a, a, c, 90, 90, 120)
    @test directioncosine([1, 0, 0], g, [1, 1, 1]) == a / 2 / sqrt(a^2 + c^2)
end
