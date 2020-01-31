module BravaisLattices

using Test

using Crystallography

@testset "Test `bravaislattices`" begin
    @test bravaislattices(
        symbol = true,
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
