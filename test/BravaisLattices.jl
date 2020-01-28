module BravaisLattices

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
