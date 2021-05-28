@testset "Fractional coordinates to Cartesian coordinates" begin
    lattice = Lattice([
        1/2 0 0
        0 1/2 0
        0 0 1
    ])
    @test CartesianFromFractional(lattice)([2, 3, 1]) == [1, 3 / 2, 1]
end
