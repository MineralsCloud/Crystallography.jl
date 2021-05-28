@testset "Fractional coordinates to Cartesian coordinates" begin
    lattice = Lattice([
        1/2 0 0
        0 1/2 0
        0 0 1
    ])
    @test CartesianFromFractional(lattice)([2, 3, 1]) == [1, 3 / 2, 1]
end

# Example from http://ww1.iucr.org/iucr-top/comm/cteach/pamphlets/22/node34.html#SECTION00073300000000000000
@testset "Refer ZrSiO₄ primitive cell with origin at 2/m" begin
    𝐚ₚ, 𝐛ₚ, 𝐜ₚ = [1, 0, 0], [0, 1, 0], [1, 1, 1] // 2
    P = PrimitiveFromStandardized(hcat(𝐚ₚ, 𝐛ₚ, 𝐜ₚ))
    @test P([0, 1 // 4, 7 // 8]) == [-7, -5, 14] // 8  # The new coordinates of the first Zr atom
    @test P([0, 1 // 4, 3 // 8]) == [-3, -1, 6] // 8  # The new coordinates of the first Si atom
    @test P([0, 0.45, 0.215]) ≈ [-0.215, 0.235, 0.43]  # The new coordinates of the first O atom
end
