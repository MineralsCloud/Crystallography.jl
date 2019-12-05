#=
TestBravaisLattices.jl:
- Julia version: 1.0
- Author: qz
- Date: Jun 2, 2019
=#
module TestBravaisLattices

using Test

using Crystallography

@testset "Test `allbravaislattices`" begin
    @test allbravaislattices(
        if_nomenclature = true,
    ) == (
        "aP",
        "mP",
        "mB",
        "oP",
        "oC",
        "oI",
        "oF",
        "tP",
        "tI",
        "cP",
        "cI",
        "cF",
        "hP",
        "hR",
    )
end

end
